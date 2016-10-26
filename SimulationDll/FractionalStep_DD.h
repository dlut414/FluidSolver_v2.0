/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Released under CC BY-NC
*/
//FractionalStep_DD.h
///defination of class FractionalStep_DD (KK)

#pragma once
#include "Simulator.h"
#include "Particle_x.h"
#include "Shifter.h"

namespace SIM {

	template <typename R, int D, int P>
	class FractionalStep_DD : public Simulator<R,D,FractionalStep_DD<R,D,P>> {};

	template <typename R, int P>
	class FractionalStep_DD<R,1,P> : public Simulator<R,1,FractionalStep_DD<R,1,P>> {};

	template <typename R, int P>
	class FractionalStep_DD<R,2,P> : public Simulator<R,2,FractionalStep_DD<R,2,P>> {
		typedef mMath::Polynomial_A<R,2,P> PN;
		typedef mMath::Derivative_A<R,2,P> DR;
		typedef Eigen::Matrix<R,PN::value,1> VecP;
		typedef Eigen::Matrix<int,2,1> iVec;
		typedef Eigen::Matrix<R,2,1> Vec;
		typedef Eigen::Triplet<R> Tpl;
		typedef Eigen::Matrix<R,PN::value,PN::value> MatPP;
	public:
		FractionalStep_DD() {}
		~FractionalStep_DD() {}

		void init_() {
			part = new Particle_x<R, 2, P>();
			part->clean();
			*part << "Geo.in";
			part->init(para.k);
			part->buildCell();
			part->makeBdc();
			part->b2b();
			part->b2normal();
			part->b2neumann();
			part->b2dirichlet();
			part->init_x();
			sen = new Sensor<R, 2, Particle_x<R, 2, P>>(part);
			*sen << "Sensor.in";
		}

		void step() {
			calInvMat();

			visTerm_i_q1();
			presTerm_i_q1();
			temperatureTerm_i_q1();

			syncPos();
			updateVelocity_q1();
			updatePressure_q1();
			updatePosition_s1();

			calCell();
			calInvMat();
			calForVis();
			check();

			Redistribute();
			sync();
		}

		void Redistribute() {
			shi.StaticUpwindModel(part);
			//shi.StaticWENOModel(part);
			//shi.SpringUpwindModel(part, para);
		}

		void visTerm_i_q1() {
			makeLhs_v_q1();
			makeRhs_v_q1();
			solvMat_v();
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == FLUID) {
					part->vel_p1[0][p] += part->vel[0][p];
					part->vel_p1[1][p] += part->vel[1][p];
				}
			}
		}

		void presTerm_i_q1() {
			makeLhs_p();
			makeRhs_p_q1();
			solvMat_phi();
		}

		void temperatureTerm_i_q1() {
			makeLhs_t();
			makeRhs_t_q1();
			solvMat_t();
		}

		void updateVelocity_q1() {
			const R coefL = para.dt;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == FLUID) {
					const Vec Gd = part->Grad(part->phi.data(), p);
					part->vel_p1[0][p] += -coefL* Gd[0];
					part->vel_p1[1][p] += -coefL* Gd[1];
				}
			}
		}

		void updatePosition_s1() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == FLUID) {
					part->pos[0][p] += para.dt * part->vel[0][p];
					part->pos[1][p] += para.dt * part->vel[1][p];
				}
			}
		}

		void updatePressure_q1() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == FLUID || part->type[p] == BD1) part->pres[p] += part->phi[p];
			}
		}

	public:
		Particle_x<R,2,P>* part;
		Sensor<R,2,Particle_x<R,2,P>>* sen;

	private:
		void makeLhs_v_q1() {
			coef.clear();
			const R coefL = R(0.5)* para.dt* para.Pr;
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD1 || part->type[p] == BD2) {
					coef.push_back(Tpl(2 * p, 2 * p, 1.0));
					coef.push_back(Tpl(2 * p + 1, 2 * p + 1, 1.0));
					continue;
				}
				R pp = R(0);
				const auto& mm = part->invMat[p];
				const auto& cell = part->cell;
				const int cx = cell->pos2cell(part->pos[0][p]);
				const int cy = cell->pos2cell(part->pos[1][p]);
				for (int i = 0; i < cell->blockSize::value; i++) {
					const int key = cell->hash(cx, cy, i);
					for (int m = 0; m < cell->linkList[key].size(); m++) {
						const int q = cell->linkList[key][m];
						if (part->type[q] == BD2) continue;
						const R dr[2] = { part->pos[0][q] - part->pos[0][p], part->pos[1][q] - part->pos[1][p] };
						const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
						if (dr1 > part->r0) continue;
						const R w = part->ww(dr1);
						VecP npq;
						part->poly(dr, npq.data());
						const VecP aa = mm * (w* npq);
						const R pq = -coefL* (part->pn_lap_o.dot(aa));
						pp -= pq;
						if (q == p) continue;
						coef.push_back(Tpl(2 * p, 2 * q, pq));
						coef.push_back(Tpl(2 * p + 1, 2 * q + 1, pq));
					}
				}
				pp += R(1.0);
				coef.push_back(Tpl(2 * p, 2 * p, pp));
				coef.push_back(Tpl(2 * p + 1, 2 * p + 1, pp));
			}
			mSol->au.setFromTriplets(coef.begin(), coef.end());
		}

		void makeRhs_v_q1() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD1 || part->type[p] == BD2) {
					mSol->rhs[2 * p + 0] = R(0);
					mSol->rhs[2 * p + 1] = R(0);
					continue;
				}
				const Vec Gd = part->Grad(part->pres.data(), p);
				const Vec Lp = part->Lap(part->vel[0].data(), part->vel[1].data(), p);
				const R rhsx = para.dt*(-Gd[0] + para.Pr* Lp[0]);
				const R rhsy = para.dt*(-Gd[1] + para.Pr* Lp[1] + para.Ra* part->temp[p]);
				mSol->rhs[2 * p + 0] = rhsx;
				mSol->rhs[2 * p + 1] = rhsy;
			}
		}

		void makeLhs_p() {
			coef.clear();
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD2) {
					coef.push_back(Tpl(p, p, 1.));
					continue;
				}
				R pqsum = R(0);
				R pp = R(0);
				MatPP* mm;
				if (IS(part->bdc[p], P_NEUMANN))	mm = &(part->invNeu.at(p));
				else								mm = &(part->invMat[p]);
				const auto& cell = part->cell;
				const int cx = cell->pos2cell(part->pos[0][p]);
				const int cy = cell->pos2cell(part->pos[1][p]);
				for (int i = 0; i < cell->blockSize::value; i++) {
					const int key = cell->hash(cx, cy, i);
					for (int m = 0; m < cell->linkList[key].size(); m++) {
						const int q = cell->linkList[key][m];
						if (part->type[q] == BD2) continue;
						const R dr[2] = { part->pos[0][q] - part->pos[0][p], part->pos[1][q] - part->pos[1][p] };
						const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
						if (dr1 > part->r0) continue;
						const R w = part->ww(dr1);
						VecP npq;
						part->poly(dr, npq.data());
						const VecP aa = (*mm) * (w* npq);
						const R pq = part->pn_lap_o* aa;
						pp -= pq;
						if (q == p) continue;
						coef.push_back(Tpl(p, q, pq));
						pqsum += abs(pq);
					}
				}
				coef.push_back(Tpl(p, p, pp));
				if (pqsum < para.eps) coef.push_back(Tpl(p, p, 1.0));
			}
			mSol->a.setFromTriplets(coef.begin(), coef.end());
		}

		void makeRhs_p_q1() {
			const R coefL = R(1) / para.dt;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD2) {
					mSol->b[p] = R(0);
					continue;
				}
				mSol->b[p] = coefL * part->Div(part->vel_p1[0].data(), part->vel_p1[1].data(), p);
				if (IS(part->bdc[p], P_NEUMANN)) {
					VecP inner = VecP::Zero();
					inner.block<2, 1>(0, 0) = part->bdnorm.at(p);
					const VecP aa = part->invNeu.at(p)* inner;
					const R cst = part->p_neumann.at(p)*part->ww(0.0)* (1.0 / part->varrho) * (part->pn_lap_o.dot(aa));
					mSol->b[p] -= cst;
				}
			}
		}

		void makeLhs_t() {
			coef.clear();
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD2) {
					coef.push_back(Tpl(p, p, 1.0));
					continue;
				}
				if (IS(part->bdc[p], T_DIRICHLET)) {
					coef.push_back(Tpl(p, p, 1.0));
					continue;
				}
				R pqsum = R(0.0);
				R pp = R(0.0);
				MatPP* mm;
				if (IS(part->bdc[p], T_NEUMANN))	mm = &(part->invNeu.at(p));
				else								mm = &(part->invMat[p]);
				const R coefL = -R(0.5)* para.dt;
				const auto& cell = part->cell;
				const int cx = cell->pos2cell(part->pos[0][p]);
				const int cy = cell->pos2cell(part->pos[1][p]);
				for (int i = 0; i < cell->blockSize::value; i++) {
					const int key = cell->hash(cx, cy, i);
					for (int m = 0; m < cell->linkList[key].size(); m++) {
						const int q = cell->linkList[key][m];
						if (part->type[q] == BD2) continue;
						const R dr[2] = { part->pos[0][q] - part->pos[0][p], part->pos[1][q] - part->pos[1][p] };
						const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
						if (dr1 > part->r0) continue;
						const R w = part->ww(dr1);
						VecP npq;
						part->poly(dr, npq.data());
						const VecP aa = (*mm) * (w* npq);
						const R pq = coefL* (part->pn_lap_o.dot(aa));
						pp -= pq;
						if (q == p) continue;
						coef.push_back(Tpl(p, q, pq));
					}
				}
				pp += 1.0;
				coef.push_back(Tpl(p, p, pp));
			}
			mSol->a.setFromTriplets(coef.begin(), coef.end());
		}

		void makeRhs_t_q1() {
			const R coefL = 0.5* para.dt;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD2) {
					mSol->b[p] = 0.0;
					continue;
				}
				if (IS(part->bdc[p], T_DIRICHLET)) {
					mSol->b[p] = part->t_dirichlet.at(p);
					continue;
				}
				mSol->b[p] = part->temp[p] + coefL* part->Lap(part->temp.data(), p);
				if (IS(part->bdc[p], T_NEUMANN)) {
					VecP inner = VecP::Zero();
					inner.block<2, 1>(0, 0) = part->bdnorm.at(p);
					const VecP aa = part->invNeu.at(p)* inner;
					const R cst = part->p_neumann.at(p)*part->ww(0.0)* (1.0 / part->varrho) * (part->pn_lap_o.dot(aa));
					mSol->b[p] -= cst;
				}
			}
		}

		void sync() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				part->vel_m1[0][p] = part->vel[0][p];
				part->vel_m1[1][p] = part->vel[1][p];
				part->vel[0][p] = part->vel_p1[0][p];
				part->vel[1][p] = part->vel_p1[1][p];
			}
		}
		void syncPos() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				part->pos_m1[0][p] = part->pos[0][p];
				part->pos_m1[1][p] = part->pos[1][p];
			}
		}

	private:
		Shifter<R,2> shi;
		std::vector<Tpl> coef;
	};

	template <typename R, int P>
	class FractionalStep_DD<R,3,P> : public Simulator<R,3,FractionalStep_DD<R,3,P>> {};

}