/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Released under CC BY-NC
*/
//FractionalStep_KM_A_BFSF.h
///defination of class FractionalStep_KM_A_ (Kim & Moin)
///Backward-facing step flow
/// FractionalStep scheme with modified Dirichlet condition for velocity PE and Neumann condition for Pressure PE
/// inlet/outlet B.C.: velocity -> Dirichlet condition    pressure -> homogeneous Neumann condition
/// inlet: velocity profile specified    outlet: fully developed flow
#pragma once
#include "Simulator.h"
#include "Particle_x.h"
#include "Shifter.h"

namespace SIM {
	
	template <typename R, int D, int P>
	class FractionalStep_KM_A_BFSF : public Simulator<R,D,FractionalStep_KM_A_BFSF<R,D,P>> {};

	template <typename R, int P>
	class FractionalStep_KM_A_BFSF<R,1,P> : public Simulator<R,1,FractionalStep_KM_A_BFSF<R,1,P>> {};

	template <typename R, int P>
	class FractionalStep_KM_A_BFSF<R,2,P> : public Simulator<R,2,FractionalStep_KM_A_BFSF<R,2,P>> {
		typedef mMath::Polynomial_A<R,2,P> PN;
		typedef mMath::Derivative_A<R,2,P> DR;
		typedef Eigen::Matrix<R,PN::value,1> VecP;
		typedef Eigen::Matrix<R,2,1> Vec;
		typedef Eigen::Matrix<R,PN::value,PN::value> MatPP;
		typedef Eigen::Triplet<R> Tpl;
	public:
		FractionalStep_KM_A_BFSF() {}
		~FractionalStep_KM_A_BFSF() {}

		void init_() {
			part = new Particle_x<R,2,P>();
			part->clean();
			*part << "Geo.in";
			part->init(para.k);
			part->buildCell();
			part->init_x();
			makeBC();
			b2neumann();
			b2dirichlet();
			calCell();
			calInvMat();
			sen = new Sensor<R, 2, Particle_x<R, 2, P>>(part);
			*sen << "Sensor.in";
			mSol = new MatSolver<R, 2, 1>(int(derived().part->np), para.eps);
		}

		void makeBC() {
			for (int p = 0; p < part->np; p++) {
				part->bdc[p] = 0;
				if (part->type[p] == BD1 || part->type[p] == INLET || part->type[p] == OUTLET) part->bdc[p] = ON(part->bdc[p], P_NEUMANN);
			}
		}
		void b2neumann() {
			for (int p = 0; p < part->np; p++) {
				if (IS(part->bdc[p], P_NEUMANN)) part->p_neumann[p] = R(0);
			}
		}
		void b2dirichlet() {
			for (int p = 0; p < part->np; p++) {
			}
		}

		void step() {
			VPE_q1r0();
			PPE_q1();
			adPres();

			syncPos();
			adVel_q1();
			adPos_s1();

			calCell();
			calInvMat();
			calForVis();
			check();

			Redistribute();
			InletOutletPart();

			sync();
			calCell();
			calInvMat();
		}

		void Redistribute() {
			//SpringLSIModel();
			SpringULSIModel();
			//StaticLSIModel();
			//StaticULSIModel();
		}

		void VPE_q1r0() {
			LHS_v_q1();
			RHS_v_q1r0_nml();
			solveMat_v();
		}

		void PPE_q1() {
			LHS_p();
			RHS_p_q1();
			solveMat_phi();
		}

		void adPres() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == FLUID || part->type[p] == BD1 || part->type[p] == INLET || part->type[p] == OUTLET) {
					part->pres[p] = part->phi[p];
				}
			}
		}

		void adVel_q1() {
			const R coef_local = para.dt;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == FLUID) {
					const Vec du = -coef_local * part->Grad(part->phi.data(), p);
					part->vel_p1[0][p] += du[0];
					part->vel_p1[1][p] += du[1];
				}
				else if (part->type[p] == BD1) {
					part->vel_p1[0][p] = part->vel[0][p];
					part->vel_p1[1][p] = part->vel[1][p];
				}
				else if (part->type[p] == INLET) {
					part->vel_p1[0][p] = part->vel[0][p];
					part->vel_p1[1][p] = part->vel[1][p];
				}
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == OUTLET) {
					int q = part->NearestFluid(p);
					const Vec gradX_q = part->Grad(part->vel[0].data(), q);
					const Vec gradY_q = part->Grad(part->vel[1].data(), q);
					const Vec norm = part->bdnorm.at(p);
					const Vec gradX_p = gradX_q - gradX_q.dot(norm) * norm;
					const Vec gradY_p = gradY_q - gradY_q.dot(norm) * norm;
					Vec Dqp;
					Dqp[0] = part->pos[0][p] - part->pos[0][q];
					Dqp[1] = part->pos[1][p] - part->pos[1][q];
					part->vel_p1[0][p] = part->vel_p1[0][q] + gradX_p.dot(Dqp);
					part->vel_p1[1][p] = part->vel_p1[1][q] + gradY_p.dot(Dqp);
					R flowRate = part->vel_p1[0][p] * norm[0] + part->vel_p1[1][p] * norm[1];
					if (flowRate < R(0)) {
						std::cout << " Flow in at OUTLET ! " << std::endl;
						std::cout << " Clamp to Zero ! " << std::endl;
						part->vel_p1[0][p] = part->vel_p1[0][p] - flowRate * norm[0];
						part->vel_p1[1][p] = part->vel_p1[1][p] - flowRate * norm[1];
						flowRate = 0;
					}
				}
			}
		}
		
		void adPos_s1() {
//			const R coef_local = R(0.5)* para.dt;
//#if OMP
//#pragma omp parallel for
//#endif
//			for (int p = 0; p < part->np; p++) {
//				if (part->type[p] == FLUID) {
//					part->pos[0][p] += coef_local * (part->vel[0][p] + part->vel_p1[0][p]);
//					part->pos[1][p] += coef_local * (part->vel[1][p] + part->vel_p1[1][p]);
//				}
//			}
			const R coef_local = para.dt;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == FLUID || part->type[p] == INLET || part->type[p] == OUTLET) {
					part->pos[0][p] += coef_local * (part->vel[0][p]);
					part->pos[1][p] += coef_local * (part->vel[1][p]);
				}
			}
		}

	public:
		MatSolver<R, 2, 1>* mSol;
		Particle_x<R, 2, P>* part;
		Sensor<R,2,Particle_x<R,2,P>>* sen;

	private:
		void LHS_v_q1() {
			coef.clear();
			R coefI_local = R(1) / (para.dt * para.Pr);
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD1 || part->type[p] == BD2 || part->type[p] == INLET || part->type[p] == OUTLET) {
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
				pp += coefI_local;
				coef.push_back(Tpl(2 * p, 2 * p, pp));
				coef.push_back(Tpl(2 * p + 1, 2 * p + 1, pp));
			}
			mSol->au.setFromTriplets(coef.begin(), coef.end());
		}

		void RHS_v_q1r0_nml() {
///@ concept 1
			const R coef_local = R(1) / (para.dt * para.Pr);
			const R epsilon_local = 1.e-6;
			R inletSum = 0;
			R outletSum = 0;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD1) {
					const Vec lap_local = part->Lap(part->vel[0].data(), part->vel[1].data(), p);
					mSol->rhs[2 * p + 0] = part->vel[0][p] + para.dt*(para.Pr* lap_local[0]);
					mSol->rhs[2 * p + 1] = part->vel[1][p] + para.dt*(para.Pr* lap_local[1]);
				}
				else if (part->type[p] == INLET) {
					const Vec lap_local = part->Lap(part->vel[0].data(), part->vel[1].data(), p);
					mSol->rhs[2 * p + 0] = part->vel[0][p] + para.dt*(para.Pr* lap_local[0]);
					mSol->rhs[2 * p + 1] = part->vel[1][p] + para.dt*(para.Pr* lap_local[1]);
					inletSum += -(part->vel[0][p] * part->bdnorm[p][0] + part->vel[1][p] * part->bdnorm[p][1]) * part->dp;
				}
				else if (part->type[p] == OUTLET) {
					int q = part->NearestFluid(p);
					const Vec gradX_q = part->Grad(part->vel[0].data(), q);
					const Vec gradY_q = part->Grad(part->vel[1].data(), q);
					const Vec norm = part->bdnorm.at(p);
					const Vec gradX_p = gradX_q - gradX_q.dot(norm) * norm;
					const Vec gradY_p = gradY_q - gradY_q.dot(norm) * norm;
					Vec Dqp;
					Dqp[0] = part->pos[0][p] - part->pos[0][q];
					Dqp[1] = part->pos[1][p] - part->pos[1][q];
					part->vel_p1[0][p] = part->vel[0][q] + gradX_p.dot(Dqp);
					part->vel_p1[1][p] = part->vel[1][q] + gradY_p.dot(Dqp);
					R flowRate = part->vel_p1[0][p] * norm[0] + part->vel_p1[1][p] * norm[1];
					outletSum += flowRate * part->dp;
				}
				else if (part->type[p] == FLUID) {
					const R rhsx = coef_local* part->vel[0][p];
					const R rhsy = coef_local* part->vel[1][p];
					mSol->rhs[2 * p + 0] = rhsx;
					mSol->rhs[2 * p + 1] = rhsy;
				}
				else {
					mSol->rhs[2 * p + 0] = 0;
					mSol->rhs[2 * p + 1] = 0;
				}
			}
			if (abs(outletSum) < epsilon_local) outletSum = epsilon_local;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == OUTLET) {
					//const Vec lap_local = part->Lap(part->vel[0].data(), part->vel[1].data(), p);
					mSol->rhs[2 * p + 0] = part->vel_p1[0][p] * abs(inletSum / outletSum);
					mSol->rhs[2 * p + 1] = part->vel_p1[1][p] * abs(inletSum / outletSum);
				}
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

		void RHS_p_q1() {
			const R coef_local = R(1) / para.dt;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD2) {
					mSol->b[p] = R(0);
					continue;
				}
				const R div_local = part->Div(part->vel_p1[0].data(), part->vel_p1[1].data(), p);
				mSol->b[p] = coef_local * div_local;
				if (IS(part->bdc[p], P_NEUMANN)) {
					if (part->type[p] == BD1) {
						Vec& normal = part->bdnorm.at(p);
						VecP inner = VecP::Zero();
						inner.block<2, 1>(0, 0) = normal;
						//const Vec lap_ustar_local = part->Lap(part->vel_p1[0].data(), part->vel_p1[1].data(), p);
						//const R neumannX = para.Pr* lap_ustar_local[0];
						//const R neumannY = para.Pr* lap_ustar_local[1];
						//const R neumann = neumannX* normal[0] + neumannY* normal[1];
						const R neumann = (normal[0] * part->vel_p1[0][p] + normal[1] * part->vel_p1[1][p]) / para.dt;
						const VecP aa = part->invNeu.at(p)* inner;
						const R cst = neumann *part->ww(R(0))* (R(1) / part->varrho) * (part->pn_lap_o.dot(aa));
						mSol->b[p] -= cst;
					}
					else if (part->type[p] == INLET) {
						///???
						///neumann != 0;
						Vec& normal = part->bdnorm.at(p);
						VecP inner = VecP::Zero();
						inner.block<2, 1>(0, 0) = normal;
						const Vec lap_ustar_local = part->Lap(part->vel_p1[0].data(), part->vel_p1[1].data(), p);
						const R neumannX = para.Pr* lap_ustar_local[0];
						const R neumannY = para.Pr* lap_ustar_local[1];
						const R neumann = neumannX* normal[0] + neumannY* normal[1];
						const VecP aa = part->invNeu.at(p)* inner;
						const R cst = neumann *part->ww(R(0))* (R(1) / part->varrho) * (part->pn_lap_o.dot(aa));
						mSol->b[p] -= cst;
					}
					else if (part->type[p] == OUTLET) {
						///neumann = 0;
					}
				}
			}
		}

		void LHS_t_q1CN() {
			coef.clear();
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD2) {
					coef.push_back(Tpl(p, p, R(1)));
					continue;
				}
				if (IS(part->bdc[p], T_DIRICHLET)) {
					coef.push_back(Tpl(p, p, R(1)));
					continue;
				}
				R pqsum = R(0);
				R pp = R(0);
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
				pp += R(1);
				coef.push_back(Tpl(p, p, pp));
			}
			mSol->a.setFromTriplets(coef.begin(), coef.end());
		}

		void RHS_t_q1CN() {
			const R coef_local = R(0.5)* para.dt;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD2) {
					mSol->b[p] = R(0);
					continue;
				}
				if (IS(part->bdc[p], T_DIRICHLET)) {
					mSol->b[p] = part->t_dirichlet.at(p);
					continue;
				}
				mSol->b[p] = part->temp[p] + coef_local* part->Lap(part->temp.data(), p);
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
				if (part->type[p] != INLET && part->type[p] != OUTLET) {
					part->pos_m1[0][p] = part->pos[0][p];
					part->pos_m1[1][p] = part->pos[1][p];
				}
			}
		}

		void InletOutletPart() {
///@ concept 3: define virtual line for OUTLET, remove distant particle, do OUTLET particle detection

///@ concept 1: fix OUTLET
			std::vector<int> inId;
			std::vector<int> rmId;
			const R dp2 = part->dp* part->dp;
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == INLET) {
					const R dr[2] = { part->pos[0][p] - part->pos_m1[0][p], part->pos[1][p] - part->pos_m1[1][p] };
					const R dr2 = dr[0] * dr[0] + dr[1] * dr[1];
					if (dr2 > dp2) {
						inId.push_back(p);
					}
				}
				else if (part->type[p] == OUTLET) {
					part->pos[0][p] = part->pos_m1[0][p];
					part->pos[1][p] = part->pos_m1[1][p];
				}
			}
			for (auto it = 0; it < inId.size(); it++) {
				const int id = inId[it];
				Vec pos, vel;
				pos[0] = part->pos[0][id];
				pos[1] = part->pos[1][id];
				vel[0] = part->vel_p1[0][id];
				vel[1] = part->vel_p1[1][id];
				part->addPart(FLUID, pos, vel, 0);
				part->pos[0][id] = part->pos_m1[0][id];
				part->pos[1][id] = part->pos_m1[1][id];
			}
			//for (int p = 0; p < part->np; p++) {
			//	if (part->type[p] == OUTLET) {
			//		const auto& cell = part->cell;
			//		const int cx = cell->pos2cell(part->pos[0][p]);
			//		const int cy = cell->pos2cell(part->pos[1][p]);
			//		for (int i = 0; i < cell->blockSize::value; i++) {
			//			const int key = cell->hash(cx, cy, i);
			//			for (int m = 0; m < cell->linkList[key].size(); m++) {
			//				const int q = cell->linkList[key][m];
			//				if (part->type[q] == FLUID) {
			//					const R dr[2] = { part->pos[0][p] - part->pos[0][q], part->pos[1][p] - part->pos[1][q] };
			//					const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
			//					if (dr1 < 0.6 * part->dp) {
			//						rmId.push_back(q);
			//					}
			//				}
			//			}
			//		}
			//	}
			//}
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == FLUID) {
					int flag = 0;
					const auto& cell = part->cell;
					const int cx = cell->pos2cell(part->pos[0][p]);
					const int cy = cell->pos2cell(part->pos[1][p]);
					for (int i = 0; i < cell->blockSize::value; i++) {
						const int key = cell->hash(cx, cy, i);
						for (int m = 0; m < cell->linkList[key].size(); m++) {
							const int q = cell->linkList[key][m];
							if (part->type[q] == OUTLET) flag = 1;
						}
					}
					if (!flag) continue;
					for (int it = 0; it < int(part->outlet.size()); it++) {
						Vec p0;
						p0[0] = part->pos[0][p];
						p0[1] = part->pos[1][p];
						const R dis = part->distance(part->outlet[it], p0);
						if (dis <= 0.3 * part->dp) rmId.push_back(p);
					}
				}
			}
			std::sort(rmId.begin(), rmId.end());
			for (int it = 0; it < (int(rmId.size()) - 1); it++) {
				while ((it+1) < rmId.size() && rmId[it + 1] == rmId[it]) {
					rmId.erase(rmId.begin() + it + 1);
				}
			}
			for (int it = (int(rmId.size()) - 1); it >= 0; it--) {
				part->erasePart(rmId[it]);
			}

///@ concept 2: move OUTLET
			//std::vector<int> inId;
			//std::vector<int> outId;
			//std::vector<int> rmId;
			//const R dp2 = part->dp* part->dp;
			//for (int p = 0; p < part->np; p++) {
			//	if (part->type[p] == INLET) {
			//		const R dr[2] = { part->pos[0][p] - part->pos_m1[0][p], part->pos[1][p] - part->pos_m1[1][p] };
			//		const R dr2 = dr[0] * dr[0] + dr[1] * dr[1];
			//		if (dr2 > dp2) {
			//			inId.push_back(p);
			//		}
			//	}
			//	else if (part->type[p] == OUTLET) {
			//		const R dr[2] = { part->pos[0][p] - part->pos_m1[0][p], part->pos[1][p] - part->pos_m1[1][p] };
			//		const R dr2 = dr[0] * dr[0] + dr[1] * dr[1];
			//		if (dr2 > dp2) {
			//			outId.push_back(p);
			//		}
			//	}
			//}
			//for (auto it = 0; it < inId.size(); it++) {
			//	const int id = inId[it];
			//	Vec pos, vel;
			//	pos[0] = part->pos[0][id];
			//	pos[1] = part->pos[1][id];
			//	vel[0] = part->vel_p1[0][id];
			//	vel[1] = part->vel_p1[1][id];
			//	part->addPart(FLUID, pos, vel, 0);
			//	part->pos[0][id] = part->pos_m1[0][id];
			//	part->pos[1][id] = part->pos_m1[1][id];
			//}
			//for (auto it = 0; it < outId.size(); it++) {
			//	const int id = outId[it];
			//	int q = part->NearestFluid(id);
			//	part->copyPart_phisicalValue(id, q);
			//	rmId.push_back(q);
			//	part->type[q] = BD2;
			//}
			//std::sort(rmId.begin(), rmId.end());
			//for (int it = int(rmId.size())-1; it >= 0; it--) {
			//	const int id = rmId[it];
			//	part->erasePart(id);
			//}
			mSol->resize(part->np);
		}

		template <int LOOP = 3>
		void SpringULSIModel() {
			std::vector<R> Dposx(part->np, R(0));
			std::vector<R> Dposy(part->np, R(0));
			std::vector<R> Du1x(part->np, R(0));
			std::vector<R> Du1y(part->np, R(0));
			std::vector<R> Du2x(part->np, R(0));
			std::vector<R> Du2y(part->np, R(0));
			const R coef = para.umax* para.dt;
			R dis_min = std::numeric_limits<R>::max();
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == OUTLET) continue;
				const auto& cell = part->cell;
				const int cx = cell->pos2cell(part->pos[0][p]);
				const int cy = cell->pos2cell(part->pos[1][p]);
				for (int i = 0; i < cell->blockSize::value; i++) {
					const int key = cell->hash(cx, cy, i);
					for (int m = 0; m < cell->linkList[key].size(); m++) {
						const int q = cell->linkList[key][m];
						if (q == p || part->type[q] == OUTLET) continue;
						const R dr[2] = { part->pos[0][q] - part->pos[0][p], part->pos[1][q] - part->pos[1][p] };
						const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
						dis_min = dr1 < dis_min ? dr1 : dis_min;
					}
				}
			}
			if (dis_min > 0.6* part->dp) return;
			std::cout << " Redistribute " << std::endl;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				Dposx[p] = part->pos[0][p];
				Dposy[p] = part->pos[1][p];
			}
			for (int iter = 0; iter < LOOP; iter++) {
#if OMP
#pragma omp parallel for
#endif
				for (int p = 0; p < part->np; p++) {
					if (part->type[p] != FLUID) continue;
					int flag = 0;
					R Dpq[2] = { 0.0, 0.0 };
					const auto& cell = part->cell;
					const int cx = cell->pos2cell(part->pos[0][p]);
					const int cy = cell->pos2cell(part->pos[1][p]);
					for (int i = 0; i < cell->blockSize::value; i++) {
						const int key = cell->hash(cx, cy, i);
						for (int m = 0; m < cell->linkList[key].size(); m++) {
							const int q = cell->linkList[key][m];
							if (q == p) continue;
							const R dr[2] = { Dposx[q] - Dposx[p], Dposy[q] - Dposy[p] };
							const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
							if (dr1 > part->r0) continue;
							if (part->type[q] == INLET || part->type[q] == OUTLET) flag = 1;
							const R w = part->ww(dr1);
							const R coeff = w / dr1;
							Dpq[0] -= coeff * dr[0];
							Dpq[1] -= coeff * dr[1];
						}
					}
					if (flag) continue;
					Dposx[p] += coef* Dpq[0];
					Dposy[p] += coef* Dpq[1];
				}
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] != FLUID) continue;
				const Vec u1 = part->interpolateLSAU(part->vel[0].data(), part->vel[1].data(), p, Dposx[p], Dposy[p]);
				Du1x[p] = u1[0];
				Du1y[p] = u1[1];
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] != FLUID) continue;
				const Vec u2 = part->interpolateLSAU(part->vel_p1[0].data(), part->vel_p1[1].data(), p, Dposx[p], Dposy[p]);
				Du2x[p] = u2[0];
				Du2y[p] = u2[1];
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] != FLUID) continue;
				part->pos[0][p] = Dposx[p];
				part->pos[1][p] = Dposy[p];
				part->vel[0][p] = Du1x[p];
				part->vel[1][p] = Du1y[p];
				part->vel_p1[0][p] = Du2x[p];
				part->vel_p1[1][p] = Du2y[p];
			}
		}

		template <int LOOP = 1>
		void SpringLSIModel() {
			std::vector<R> Dposx(part->np, R(0));
			std::vector<R> Dposy(part->np, R(0));
			std::vector<R> Du1x(part->np, R(0));
			std::vector<R> Du1y(part->np, R(0));
			std::vector<R> Du2x(part->np, R(0));
			std::vector<R> Du2y(part->np, R(0));
			const R coef = para.umax* para.dt;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				Dposx[p] = part->pos[0][p];
				Dposy[p] = part->pos[1][p];
			}
			for (int iter = 0; iter < LOOP; iter++) {
#if OMP
#pragma omp parallel for
#endif
				for (int p = 0; p < part->np; p++) {
					if (part->type[p] != FLUID) continue;
					R Dpq[2] = { 0.0, 0.0 };
					const auto& cell = part->cell;
					const int cx = cell->pos2cell(part->pos[0][p]);
					const int cy = cell->pos2cell(part->pos[1][p]);
					for (int i = 0; i < cell->blockSize::value; i++) {
						const int key = cell->hash(cx, cy, i);
						for (int m = 0; m < cell->linkList[key].size(); m++) {
							const int q = cell->linkList[key][m];
							if (q == p) continue;
							const R dr[2] = { Dposx[q] - Dposx[p], Dposy[q] - Dposy[p] };
							const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
							if (dr1 > part->r0) continue;
							const R w = part->ww(dr1);
							const R coeff = w / dr1;
							Dpq[0] -= coeff * dr[0];
							Dpq[1] -= coeff * dr[1];
						}
					}
					Dposx[p] += coef* Dpq[0];
					Dposy[p] += coef* Dpq[1];
				}
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] != FLUID) continue;
				const Vec u1 = part->interpolateLSA(part->vel[0].data(), part->vel[1].data(), p, Dposx[p], Dposy[p]);
				Du1x[p] = u1[0];
				Du1y[p] = u1[1];
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] != FLUID) continue;
				const Vec u2 = part->interpolateLSA(part->vel_p1[0].data(), part->vel_p1[1].data(), p, Dposx[p], Dposy[p]);
				Du2x[p] = u2[0];
				Du2y[p] = u2[1];
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] != FLUID) continue;
				part->pos[0][p] = Dposx[p];
				part->pos[1][p] = Dposy[p];
				part->vel[0][p] = Du1x[p];
				part->vel[1][p] = Du1y[p];
				part->vel_p1[0][p] = Du2x[p];
				part->vel_p1[1][p] = Du2y[p];
			}
		}

	private:
		Shifter<R,2> shi;
		std::vector<Tpl> coef;
	};

	template <typename R, int P>
	class FractionalStep_KM_A_BFSF<R,3,P> : public Simulator<R,3,FractionalStep_KM_A_BFSF<R,3,P>>  {};
}