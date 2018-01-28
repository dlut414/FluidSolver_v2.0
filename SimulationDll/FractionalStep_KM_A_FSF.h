/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Released under CC BY-NC
*/
//FractionalStep_KM_A_FSF.h
///defination of class FractionalStep_KM_A_ (Kim & Moin)
///FSF: Free-Surface Flow
/// FractionalStep scheme
/// free-surface B.C.: pressure -> zero Dirichlet condition
/// wall boundary B.C.: pressure -> homogeneous Neumann condition
#pragma once
#include "Simulator.h"
#include "Particle_x.h"
#include "Shifter.h"
#include <algorithm>
#include <queue>
#define BOOST_PYTHON_DYNAMIC_LIB
#define BOOST_NUMPY_DYNAMIC_LIB
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <Python.h>

namespace SIM {

	template <typename R, int D, int P>
	class FractionalStep_KM_A_FSF : public Simulator<R, D, FractionalStep_KM_A_FSF<R, D, P>> {};

	template <typename R, int P>
	class FractionalStep_KM_A_FSF<R, 1, P> : public Simulator<R, 1, FractionalStep_KM_A_FSF<R, 1, P>>{};

	template <typename R, int P>
	class FractionalStep_KM_A_FSF<R, 2, P> : public Simulator<R, 2, FractionalStep_KM_A_FSF<R, 2, P>>{
		typedef mMath::Polynomial_A<R, 2, P> PN;
		typedef mMath::Derivative_A<R, 2, P> DR;
		typedef Eigen::Matrix<R, PN::value, 1> VecP;
		typedef Eigen::Matrix<R, 2, 1> Vec;
		typedef Eigen::Matrix<R, PN::value, PN::value> MatPP;
		typedef Eigen::Triplet<R> Tpl;
		typedef MatSolver<R, 2, 0> Solver;
		typedef Particle_x<R, 2, P> PartX;
	public:
		FractionalStep_KM_A_FSF() {}
		~FractionalStep_KM_A_FSF() {}

		void init_() {
			part = new PartX();
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
			fs = std::vector<R>(part->np, 0);
			sen = new Sensor<R, 2, PartX>(part);
			*sen << "Sensor.in";
			mSol = new Solver(int(derived().part->np), para.eps);
		}

		R* Surface() {
			return fs.data();
		}

		void makeBC() {
			for (int p = 0; p < part->np; p++) {
				part->bdc[p] = 0;
				if (part->type[p] == BD1 || part->type[p] == INLET) part->bdc[p] = ON(part->bdc[p], P_NEUMANN);
			}
		}
		void b2neumann() {
		}
		void b2dirichlet() {
		}

		void step() {
			makeFs();
			VPE_q1r0();
			PPE_q1();

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
			nn_ULSIModel();
		}

		void VPE_q1r0() {
			LHS_v_q1();
			RHS_v_q1r0();
			solveMat_v();
		}

		void PPE_q1() {
			LHS_p();
			RHS_p_q1();
			solveMat_p();
		}

		void adVel_q1() {
			const R coef_local = para.dt;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == FLUID) {
					const Vec du = -coef_local * part->Grad(part->pres.data(), p);
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
				else if (part->type[p] == OUTLET) {
				}
			}
		}

		void adPos_s1() {
			const R coef_local = R(0.5)* para.dt;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == FLUID) {
					part->pos[0][p] += coef_local * (part->vel[0][p] + part->vel_p1[0][p]);
					part->pos[1][p] += coef_local * (part->vel[1][p] + part->vel_p1[1][p]);
				}
			}
//			const R coef_local = para.dt;
//#if OMP
//#pragma omp parallel for
//#endif
//			for (int p = 0; p < part->np; p++) {
//				if (part->type[p] == FLUID || part->type[p] == INLET || part->type[p] == OUTLET) {
//					part->pos[0][p] += coef_local * (part->vel[0][p]);
//					part->pos[1][p] += coef_local * (part->vel[1][p]);
//				}
//			}
		}

	public:
		Solver* mSol;
		PartX* part;
		Sensor<R, 2, Particle_x<R, 2, P>>* sen;

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

		void RHS_v_q1r0() {
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
					const Vec gradX_q = part->Grad(part->vel[0].data(), q, FLUID | BD1 | OUTLET);
					const Vec gradY_q = part->Grad(part->vel[1].data(), q, FLUID | BD1 | OUTLET);
					const Vec norm = part->bdnorm.at(p);
					const Vec gradX_p = gradX_q - gradX_q.dot(norm) * norm;
					const Vec gradY_p = gradY_q - gradY_q.dot(norm) * norm;
					Vec Dqp;
					Dqp[0] = part->pos[0][p] - part->pos[0][q];
					Dqp[1] = part->pos[1][p] - part->pos[1][q];
					part->vel_p1[0][p] = part->vel[0][q] + gradX_p.dot(Dqp);
					part->vel_p1[1][p] = part->vel[1][q] + gradY_p.dot(Dqp);
					R flowRate = part->vel_p1[0][p] * norm[0] + part->vel_p1[1][p] * norm[1];
					if (flowRate < 0) {
						std::cout << " minus flowRate " << std::endl;
						PRINT(p);
						PRINT(flowRate);
					}
					outletSum += flowRate * part->dp;
				}
				else if (part->type[p] == FLUID) {
					/// gravity is adjusted here
					const R rhsx = coef_local* part->vel[0][p];
					const R rhsy = coef_local* part->vel[1][p] + (1.0 / para.Pr);
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
					//makeup for pressure gradient
					//const Vec Du = part->Grad(part->pres.data(), p, FLUID | BD1);
					mSol->rhs[2 * p + 0] = part->vel_p1[0][p];// +para.dt* Du[0];// *abs(inletSum / outletSum);
					mSol->rhs[2 * p + 1] = part->vel_p1[1][p];// +para.dt* Du[1];// *abs(inletSum / outletSum);
				}
			}
		}

		void LHS_p() {
			coef.clear();
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD2 || part->type[p] == OUTLET || fs[p] > 0.5) {
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
				if (part->type[p] == BD2 || fs[p] > 0.5) {
					mSol->b[p] = R(0);
					continue;
				}
				else if (part->type[p] == OUTLET) {
					///Dirichlet condition by extrapolation
					int q = part->NearestFluid(p);
					const Vec grad = part->Grad(part->pres.data(), q, FLUID|BD1);
					Vec Dqp;
					Dqp[0] = part->pos[0][p] - part->pos[0][q];
					Dqp[1] = part->pos[1][p] - part->pos[1][q];
					mSol->b[p] = part->pres[q] + grad.dot(Dqp);
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
					//else if (part->type[p] == OUTLET) {
					//	///neumann = 0;
					//}
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
					if (dr2 > dp2) inId.push_back(p);
				}
				else if (part->type[p] == OUTLET) {
					///make outlet able to inject particles
					//const R dr[2] = { part->pos[0][p] - part->pos_m1[0][p], part->pos[1][p] - part->pos_m1[1][p] };
					//const R dr2 = dr[0] * dr[0] + dr[1] * dr[1];
					//const Vec norm = part->bdnorm.at(p);
					//R flowRate = part->vel_p1[0][p] * norm[0] + part->vel_p1[1][p] * norm[1];
					//if (flowRate > 0) {
					//	part->pos[0][p] = part->pos_m1[0][p];
					//	part->pos[1][p] = part->pos_m1[1][p];
					//}
					//if (dr2 > dp2) inId.push_back(p);
					///previous code
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
						//remove particle threshold
						if (dis <= 0.45 * part->dp) rmId.push_back(p);
					}
				}
			}
			std::sort(rmId.begin(), rmId.end());
			for (int it = 0; it < (int(rmId.size()) - 1); it++) {
				while ((it + 1) < rmId.size() && rmId[it + 1] == rmId[it]) {
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
		void nn_ULSIModel() {
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
			if (dis_min > 0.7* part->dp) return;
			///distribute when dis_min <= 0.7* dp
			///Dposx hold new pos
			std::cout << " Redistribute " << std::endl;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				Dposx[p] = part->pos[0][p];
				Dposy[p] = part->pos[1][p];
			}

			using namespace boost::python;
			namespace NPY = boost::python::numpy;
			if (!python_initialized) {
				Py_Initialize();
				python_initialized = true;
			}
			if (!numpy_initialized) {
				NPY::initialize();
				main_module = import("__main__");
				global = main_module.attr("__dict__");
				numpy_initialized = true;
			}
			try {
				exec("import numpy as np", global, global);
				exec("import tensorflow as tf", global, global);
				exec("sess = tf.Session()", global, global);
				exec("saver = tf.train.import_meta_graph('./python/tf_model/model.meta')", global, global);
				exec("saver.restore(sess, tf.train.latest_checkpoint('./python/tf_model/'))", global, global);
				exec("graph = tf.get_default_graph()", global, global);
				exec("x = graph.get_tensor_by_name('x:0')", global, global);
				exec("y = graph.get_tensor_by_name('y:0')", global, global);
			}
			catch (const error_already_set&) {
				PyErr_Print();
			}
			const auto& pos = part->pos;
			const auto& dp = part->dp;
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] != FLUID || fs[p] > 0.5) continue;
				const int N = 24;
				std::vector<int> nbr;
				nNearestNeighbor<N>(nbr, p);
				if (!numpy_initialized) {
					NPY::initialize();
					numpy_initialized = true;
				}
				NPY::ndarray xx = NPY::zeros(make_tuple(1, 2 * N), NPY::dtype::get_builtin<float>());
				for (size_t i = 0; i < nbr.size(); i++) {
					xx[0][i * 2] = (pos[0][nbr[i]] - pos[0][p]) / dp;
					xx[0][i * 2 + 1] = (pos[1][nbr[i]] - pos[1][p]) / dp;
				}
				for (size_t i = nbr.size(); i < N; i++) {
					xx[0][i * 2] = 0;
					xx[0][i * 2 + 1] = 0;
				}
				///predict here
				try {
					object x = global["x"];
					object y = global["y"];
					dict feed_dict;
					feed_dict[x] = xx;
					///sess argument to y.attr("eval")() must be defined as object first
					object sess = global["sess"];
					object res = y.attr("eval")(feed_dict, sess);
					object item = res.attr("item");
					R vx = extract<R>(item(0));
					R vy = extract<R>(item(1));
					Dposx[p] += dp* vx;
					Dposy[p] += dp* vy;
				}
				catch (const error_already_set&) {
					PyErr_Print();
				}
			}
			exec("sess.close()", global, global);

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

		void makeFs() {
			using namespace boost::python;
			namespace NPY = boost::python::numpy;
			if (!python_initialized) {
				Py_Initialize();
				main_module = import("__main__");
				global = main_module.attr("__dict__");
				python_initialized = true;
			}
			exec("import sys", global, global);
			exec("sys.path.append('./python')", global, global);
			exec("import numpy", global, global);
			exec("from NN import NN as NN", global, global);
			exec("Layers = (17, 8, 8, 1)", global, global);
			exec("nn = NN(Layers = Layers)", global, global);
			exec("nn.load('./python/config')", global, global);
			const auto& pos = part->pos;
			for (int p = 0; p < part->np; p++) {
				static const int N = 8;
				std::vector<int> nbr;
				nNearestNeighbor<N>(nbr, p);
				if (!numpy_initialized) {
					NPY::initialize();
					numpy_initialized = true;
				}
				NPY::ndarray x = NPY::zeros(make_tuple(2 * N + 1, 1), NPY::dtype::get_builtin<float>());
				x[0][0] = part->type[p];
				for (size_t i = 0; i < nbr.size(); i++) {
					x[i * 2 + 1][0] = (pos[0][nbr[i]] - pos[0][p]) / part->dp;
					x[i * 2 + 2][0] = (pos[1][nbr[i]] - pos[1][p]) / part->dp;
				}
				for (size_t i = nbr.size(); i < N; i++) {
					x[i * 2 + 1][0] = 0;
					x[i * 2 + 2][0] = 0;
				}
				object nn = global["nn"];
				object predict01 = nn.attr("predict01");
				object ret = predict01(x);
				object np = global["numpy"];
				object asscalar = np.attr("asscalar");
				object iret = asscalar(ret);
				fs[p] = R(extract<int>(iret));
			}
		}

		template <int N = 8>
		void nNearestNeighbor(std::vector<int>& nbr, const int p) const {
			typedef std::pair<R, int> R2i;
			std::priority_queue<R2i, std::vector<R2i>, std::greater<R2i>> que;
			const auto& cell = part->cell;
			const auto& pos = part->pos;
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (part->type[q] == BD2) continue;
					const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
					const R dr2 = dr[0] * dr[0] + dr[1] * dr[1];
					que.push({ dr2, q });
				}
			}
			for (int i = 0; i < N; i++) {
				if (que.empty()) break;
				nbr.push_back(que.top().second);
				que.pop();
			}
		}

	private:
		Shifter<R, 2> shi;
		std::vector<Tpl> coef;
		std::vector<R> fs;
		///for python
		bool python_initialized;
		bool numpy_initialized;
		boost::python::object main_module;
		boost::python::object global;
	};

	template <typename R, int P>
	class FractionalStep_KM_A_FSF<R, 3, P> : public Simulator<R, 3, FractionalStep_KM_A_FSF<R, 3, P>>{};
}