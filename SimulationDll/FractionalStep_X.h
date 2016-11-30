/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Released under CC BY-NC
*/
//FractionalStep_X.h
///defination of class FractionalStep_X
/// solve pressure term first then viscousity term; 
#pragma once
#include "Simulator.h"
#include "Particle_x.h"
#include "Shifter.h"

namespace SIM {
	
	template <typename R, int D, int P>
	class FractionalStep_X : public Simulator<R,D,FractionalStep_X<R,D,P>> {};

	template <typename R, int P>
	class FractionalStep_X<R,1,P> : public Simulator<R,1,FractionalStep_X<R,1,P>> {};

	template <typename R, int P>
	class FractionalStep_X<R,2,P> : public Simulator<R,2,FractionalStep_X<R,2,P>> {
		typedef mMath::Polynomial_A<R,2,P> PN;
		typedef mMath::Derivative_A<R,2,P> DR;
		typedef Eigen::Matrix<R,PN::value,1> VecP;
		typedef Eigen::Matrix<R,2,1> Vec;
		typedef Eigen::Matrix<R,PN::value,PN::value> MatPP;
		typedef Eigen::Triplet<R> Tpl;
	public:
		FractionalStep_X() {}
		~FractionalStep_X() {}

		void init_() {
			part = new Particle_x<R,2,P>();
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
			sen = new Sensor<R,2,Particle_x<R,2,P>>(part);
			*sen << "Sensor.in";
			Div_old.resize(part->np);
		}

		void step() {
			calInvMat();

			TPE_CN1();
			PPE_q2();
			VPE_q2r0();

			syncPos();
			adVel_q2();
			adPos_s2();

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

		void VPE_q2r1() {
			LHS_v_q2();
			RHS_v_q2r1();
			solveMat_v();
		}

		void VPE_q1r0() {
			LHS_v_q1();
			RHS_v_q1r0();
			solveMat_v();
		}

		void VPE_q2r0() {
			LHS_v_q2();
			RHS_v_q2r0();
			solveMat_v();
		}

		void PPE_q2() {
			LHS_p();
			RHS_p_q2();
			solveMat_phi();
			//for (int p = 0; p < part->np; p++) {
			//	//if (p == 165) PRINT( (R(3.0) / (R(2.0)* para.dt))* part->Div(part->vel_p1[0].data(), part->vel_p1[1].data(), p) - part->Lap(part->phi.data(), p) - mSol->x[part->np] );
			//	//if (p == 165) PRINT( (mSol->a.row(p).dot( mSol->x)) - mSol->b[p] );
			//}
		}

		void PPE_q1() {
			LHS_p();
			RHS_p_q1();
			solveMat_phi();
		}

		void TPE_CN1() {
			LHS_t_CN1();
			RHS_t_CN1();
			solveMat_t();
		}

		void adVel_q1() {
			const R coefL = para.dt;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == FLUID || part->type[p] == BD1) {
					part->pres[p] = part->phi[p];
				}
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == FLUID) {
					const Vec du = -coefL * part->Grad(part->phi.data(), p);
					part->vel_p1[0][p] += du[0];
					part->vel_p1[1][p] += du[1];
				}
			}
		}

		void adVel_q2() {
			const R coefL = (R(2)* para.dt) / (R(3));
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				part->pres[p] = part->phi[p];
			}
		}

		void adPos_s1() {
			const R coefL = R(0.5)* para.dt;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == FLUID) {
					part->pos[0][p] += coefL * (part->vel[0][p] + part->vel_p1[0][p]);
					part->pos[1][p] += coefL * (part->vel[1][p] + part->vel_p1[1][p]);
				}
			}
		}

		void adPos_s2() {
			const R coefL = 0.5* para.dt;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == FLUID) {
					part->pos[0][p] += coefL * (3.0* part->vel[0][p] - 1.0* part->vel_m1[0][p]);
					part->pos[1][p] += coefL * (3.0* part->vel[1][p] - 1.0* part->vel_m1[1][p]);
				}
			}
		}

	public:
		Particle_x<R,2,P>* part;
		Sensor<R,2,Particle_x<R,2,P>>* sen;

	private:
		void LHS_v_q2() {
			coef.clear();
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD1 || part->type[p] == BD2) {
					coef.push_back(Tpl(2 * p, 2 * p, R(1)));
					coef.push_back(Tpl(2 * p + 1, 2 * p + 1, R(1)));
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
						const R pq = -(part->pn_lap_o.dot(aa));
						pp -= pq;
						if (q == p) continue;
						coef.push_back(Tpl(2 * p, 2 * q, pq));
						coef.push_back(Tpl(2 * p + 1, 2 * q + 1, pq));
					}
				}
				pp += 3.0 / (2.0 * para.dt * para.Pr);
				coef.push_back(Tpl(2 * p, 2 * p, pp));
				coef.push_back(Tpl(2 * p + 1, 2 * p + 1, pp));
			}
			mSol->au.setFromTriplets(coef.begin(), coef.end());
		}

		void LHS_v_q1() {
			coef.clear();
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
						const R pq = -(part->pn_lap_o.dot(aa));
						pp -= pq;
						if (q == p) continue;
						coef.push_back(Tpl(2 * p, 2 * q, pq));
						coef.push_back(Tpl(2 * p + 1, 2 * q + 1, pq));
					}
				}
				pp += 1.0 / (para.dt * para.Pr);
				coef.push_back(Tpl(2 * p, 2 * p, pp));
				coef.push_back(Tpl(2 * p + 1, 2 * p + 1, pp));
			}
			mSol->au.setFromTriplets(coef.begin(), coef.end());
		}

		void RHS_v_q2r1() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD1 || part->type[p] == BD2) {
					mSol->rhs[2 * p + 0] = part->vel[0][p];
					mSol->rhs[2 * p + 1] = part->vel[1][p];
					continue;
				}
				const Vec Gp = part->Grad(part->pres, p);
				const R coefL = 1.0 / (2.0* para.dt * para.Pr);
				const R rhsx = coefL* (4.0* part->vel[0][p] - part->vel_m1[0][p]) - (1.0 / para.Pr)* Gp[0];
				const R rhsy = coefL* (4.0* part->vel[1][p] - part->vel_m1[1][p]) - (1.0 / para.Pr)* Gp[1] + para.Ra* part->temp[p];
				mSol->rhs[2 * p + 0] = rhsx;
				mSol->rhs[2 * p + 1] = rhsy;
			}
		}

		void RHS_v_q1r0() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD1 || part->type[p] == BD2) {
					mSol->rhs[2 * p + 0] = part->vel[0][p];
					mSol->rhs[2 * p + 1] = part->vel[1][p];
					continue;
				}
				const R coefL = (R(1.0) / para.dt * para.Pr);
				const R rhsx = coefL* part->vel[0][p];
				const R rhsy = coefL* part->vel[1][p] + para.Ra* part->temp[p];
				mSol->rhs[2 * p + 0] = rhsx;
				mSol->rhs[2 * p + 1] = rhsy;
			}
		}

		void RHS_v_q2r0() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD1 || part->type[p] == BD2) {
					mSol->rhs[2 * p + 0] = part->vel[0][p];
					mSol->rhs[2 * p + 1] = part->vel[1][p];
					continue;
				}
				const R coef_vel_local = R(1) / (R(2)* para.dt * para.Pr);
				const R coef_pres_local = R(1) / para.Pr;
				const Vec p_grad_local = part->Grad(part->phi.data(), p);
				const R rhsx = coef_vel_local* (R(4)* part->vel[0][p] - part->vel_m1[0][p]) - coef_pres_local* p_grad_local[0];
				const R rhsy = coef_vel_local* (R(4)* part->vel[1][p] - part->vel_m1[1][p]) + para.Ra* part->temp[p] - coef_pres_local* p_grad_local[1];
				mSol->rhs[2 * p + 0] = rhsx;
				mSol->rhs[2 * p + 1] = rhsy;
			}
		}


		void LHS_p() {
			coef.clear();
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD2) {
					coef.push_back(Tpl(p, p, R(1)));
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
				if (pqsum < para.eps) coef.push_back(Tpl(p, p, R(1.0)));
			}
			mSol->a.setFromTriplets(coef.begin(), coef.end());
		}

		void RHS_p_q2() {
			const R RaPr = para.Ra* para.Pr;
			const R one_over_dt = R(1) / para.dt;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				const R pt_py_local = part->DerY(part->temp.data(), p);
				const R div_u_local = part->Div(part->vel[0].data(), part->vel[1].data(), p);
				const R div_um1_local = part->Div(part->vel_m1[0].data(), part->vel_m1[1].data(), p);
				mSol->b[p] = RaPr * pt_py_local + one_over_dt* (R(2)* div_u_local - R(0.5)* div_um1_local);
				if (IS(part->bdc[p], P_NEUMANN)) {
					Vec& normal = part->bdnorm.at(p);
					VecP inner = VecP::Zero();
					inner.block<2, 1>(0, 0) = normal;
					const Vec lap_ustar_local = part->Lap(part->vel_p1[0].data(), part->vel_p1[1].data(), p);
					const R neumannX = para.Pr* lap_ustar_local[0];
					const R neumannY = para.Pr* lap_ustar_local[1] + para.Ra* para.Pr* part->temp[p];
					const R neumann = neumannX* normal[0] + neumannY* normal[1];
					const VecP aa = part->invNeu.at(p)* inner;
					const R cst = neumann *part->ww(R(0))* (R(1) / part->varrho) * (part->pn_lap_o.dot(aa));
					mSol->b[p] -= cst;
				}
			}
		}

		void RHS_p_q1() {
			const R coefL = R(1.0) / para.dt;
			part->DivGrad(part->phi.data(), DGP_old.data());
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD2) {
					mSol->b[p] = 0.0;
					continue;
				}
				const R div_local = part->Div(part->vel_p1[0].data(), part->vel_p1[1].data(), p);
				const R LP_old_local = part->Lap(part->phi.data(), p);
				mSol->b[p] = coefL * div_local + (LP_old_local - DGP_old[p]);
				if (IS(part->bdc[p], P_NEUMANN)) {
					VecP inner = VecP::Zero();
					inner.block<2, 1>(0, 0) = part->bdnorm.at(p);
					const VecP aa = part->invNeu.at(p)* inner;
					const R cst = part->p_neumann.at(p)*part->ww(0.0)* (1.0 / part->varrho) * (part->pn_lap_o.dot(aa));
					mSol->b[p] -= cst;
				}
			}
		}

		void LHS_t_CN1() {
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

		void RHS_t_CN1() {
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
					inner.block<2,1>(0, 0) = part->bdnorm.at(p);
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
		std::vector<R> Div_old;
	};

	template <typename R, int P>
	class FractionalStep_X<R,3,P> : public Simulator<R,3,FractionalStep_X<R,3,P>>  {};
}