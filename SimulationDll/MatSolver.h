/*
*/
#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <Eigen/OrderingMethods>
#include "Header.h"

#define AUGMENT (1)
#define AG AUGMENT

namespace SIM {

	template <typename R, int D>
	class MatSolver {
		typedef Eigen::Triplet<R> Tpl;
		typedef Eigen::Matrix<R, Eigen::Dynamic, 1> dVec;
		typedef Eigen::SparseMatrix<R, Eigen::RowMajor> sMat;

		typedef Eigen::IncompleteLUT<R> preconditionerILU;
		typedef Eigen::DiagonalPreconditioner<R> preconditionerDia;
	public:
		MatSolver(const int& _n, const R& e) : n(_n), Dn(D*_n), eps(e) {
			init();
		}
		~MatSolver() {}

		void resize(const int& nsize) {
			n = nsize, Dn = D* nsize;
			a.resize(n + AG, n + AG), x.resize(n + AG), b.resize(n + AG);
			au.resize(D*n, D*n), u.resize(D*n), rhs.resize(D*n);
			vr.resize(n + AG), vr_hat.resize(n + AG), vv.resize(n + AG), vp.resize(n + AG), vt.resize(n + AG), vh.resize(n + AG), vs.resize(n + AG);
			Dvr.resize(D*n), Dvr_hat.resize(D*n), Dvv.resize(D*n), Dvp.resize(D*n), Dvt.resize(D*n), Dvh.resize(D*n), Dvs.resize(D*n);
		}

		void biCg() {
			solverBiCgDia.compute(a);
			//x = solverBiCgDia.solveWithGuess(b, x);
			x = solverBiCgDia.solve(b);
			std::cout << " iterations ----------> " << solverBiCgDia.iterations() << std::endl;
			std::cout << " error ---------------> " << solverBiCgDia.error() << std::endl;
//			R rho;
//			R alpha;
//			R beta;
//			R omega;
//			R residual;
//			int loop = 0;
//			sMat pre(n + AG, n + AG);
//			std::vector<Tpl> coefList;
//			for (int p = 0; p < n; p++) {
//				coefList.push_back(Tpl(p, p, R(1.0) / a.coeff(p, p)));
//				x[p] = R(0.0);
//			}
//#if AUGMENT
//			coefList.push_back(Tpl(n,n,R(1.0)));
//			x[n] = b[n] = R(0.0);
//#endif
//			pre.setFromTriplets(coefList.begin(), coefList.end());
//			a = pre * a;
//			b = pre * b;
//			a.coeffRef(n, n) = R(1.0);
//			vr = b - a*x;
//			residual = sqrt(vr.dot(vr));
//			if (residual < eps) {
//				std::cout << " iterations ----------> " << loop << std::endl;
//				std::cout << " error ---------------> " << residual << std::endl;
//				return;
//			}
//#if OMP
//#pragma omp parallel for
//#endif
//			for (int p = 0; p < n + AG; p++) {
//				vr_hat[p] = vr[p];
//				vv[p] = vp[p] = R(0.0);
//			}
//			rho = alpha = omega = R(1.0);
//			for (loop = 0; loop < maxIter; loop++) {
//				const R rho_m1 = rho;
//				rho = vr_hat.dot(vr);
//				if (abs(rho) < eps) rho = eps;
//				beta = (rho / rho_m1)*(alpha / omega);
//#if OMP
//#pragma omp parallel for
//#endif
//				for (int p = 0; p < n + AG; p++) {
//					vp[p] = vr[p] + beta*(vp[p] - omega* vv[p]);
//				}
//				vv = a* vp;
//				R vr_hat_vv = vr_hat.dot(vv);
//				if (abs(vr_hat_vv) < eps) vr_hat_vv = eps;
//				alpha = rho / vr_hat_vv;
//#if OMP
//#pragma omp parallel for
//#endif
//				for (int p = 0; p < n + AG; p++) {
//					vh[p] = x[p] + alpha*vp[p];
//					vs[p] = vr[p] - alpha* vv[p];
//				}
//				vt = b - a * vh;
//				residual = sqrt(vt.dot(vt));
//				if (residual < eps) {
//#if OMP
//#pragma omp parallel for
//#endif
//					for (int p = 0; p < n + AG; p++) {
//						x[p] = vh[p];
//					}
//					std::cout << " iterations ----------> " << loop << std::endl;
//					std::cout << " error ---------------> " << residual << std::endl;
//					return;
//				}
//				vt = a*vs;
//				omega = (vt.dot(vs)) / (vt.dot(vt));
//#if OMP
//#pragma omp parallel for
//#endif
//				for (int p = 0; p < n + AG; p++) {
//					x[p] = vh[p] + omega*vs[p];
//				}
//				vr = b - a*x;
//				residual = sqrt(vr.dot(vr));
//				if (residual < eps) break;
//#if OMP
//#pragma omp parallel for
//#endif
//				for (int p = 0; p < n + AG; p++) {
//					vr[p] = vs[p] - omega*vt[p];
//				}
//			}
//			std::cout << " iterations ----------> " << loop << std::endl;
//			std::cout << " error ---------------> " << residual << std::endl;
		}
		void biCg_v() {
			solverBiCgDia.compute(au);
			u = solverBiCgDia.solveWithGuess(rhs, u);
			u = solverBiCgDia.solve(rhs);
			std::cout << " iterations ----------> " << solverBiCgDia.iterations() << std::endl;
			std::cout << " error ---------------> " << solverBiCgDia.error() << std::endl;
//			R rho;
//			R alpha;
//			R beta;
//			R omega;
//			R residual;
//			int loop = 0;
//			sMat pre(Dn, Dn);
//			std::vector<Tpl> coefList;
//			for (int p = 0; p < Dn; p++) {
//				coefList.push_back(Tpl(p, p, R(1.0) / au.coeff(p, p)));
//				u[p] = R(0.0);
//			}
//			pre.setFromTriplets(coefList.begin(), coefList.end());
//			au = pre * au;
//			rhs = pre * rhs;
//			Dvr = rhs - au*u;
//			residual = sqrt(Dvr.dot(Dvr));
//			if (residual < eps) {
//				std::cout << " iterations ----------> " << loop << std::endl;
//				std::cout << " error ---------------> " << residual << std::endl;
//				return;
//			}
//#if OMP
//#pragma omp parallel for
//#endif
//			for (int p = 0; p < Dn; p++) {
//				Dvr_hat[p] = Dvr[p];
//				Dvv[p] = Dvp[p] = R(0.0);
//			}
//			rho = alpha = omega = R(1.0);
//			for (loop = 0; loop < maxIter; loop++) {
//				const R rho_m1 = rho;
//				rho = Dvr_hat.dot(Dvr);
//				beta = (rho / rho_m1)*(alpha / omega);
//#if OMP
//#pragma omp parallel for
//#endif
//				for (int p = 0; p < Dn; p++) {
//					Dvp[p] = Dvr[p] + beta*(Dvp[p] - omega* Dvv[p]);
//				}
//				Dvv = au* Dvp;
//				alpha = rho / (Dvr_hat.dot(Dvv));
//#if OMP
//#pragma omp parallel for
//#endif
//				for (int p = 0; p < Dn; p++) {
//					Dvh[p] = u[p] + alpha*Dvp[p];
//					Dvs[p] = Dvr[p] - alpha* Dvv[p];
//				}
//				Dvt = rhs - au * Dvh;
//				residual = sqrt(Dvt.dot(Dvt));
//				if (residual < eps) {
//#if OMP
//#pragma omp parallel for
//#endif
//					for (int p = 0; p < Dn; p++) {
//						u[p] = Dvh[p];
//					}
//					std::cout << " iterations ----------> " << loop << std::endl;
//					std::cout << " error ---------------> " << residual << std::endl;
//					return;
//				}
//				Dvt = au*Dvs;
//				omega = (Dvt.dot(Dvs)) / (Dvt.dot(Dvt));
//#if OMP
//#pragma omp parallel for
//#endif
//				for (int p = 0; p < Dn; p++) {
//					u[p] = Dvh[p] + omega*Dvs[p];
//				}
//				Dvr = rhs - au*u;
//				residual = sqrt(Dvr.dot(Dvr));
//				if (residual < eps) break;
//#if OMP
//#pragma omp parallel for
//#endif
//				for (int p = 0; p < Dn; p++) {
//					Dvr[p] = Dvs[p] - omega*Dvt[p];
//				}
//			}
//			std::cout << " iterations ----------> " << loop << std::endl;
//			std::cout << " error ---------------> " << residual << std::endl;
		}
		void qr() {
			a.makeCompressed();
			solverQR.compute(a);
			x = solverQR.solve(b);
			std::cout << " rank ----------------> " << solverQR.rank() << std::endl;
			if (solverQR.info() != Eigen::Success) {
				std::cout << " info ----------------> " << solverQR.lastErrorMessage() << std::endl;
			}
		}
		void lnBiCg() {
			sMat at(n, n);
			at = a.transpose();
			a = (a*at);
			solverBiCgDia.compute(a);
			//x = solver1.solveWithGuess(b, 0.5*x);
			x = solverBiCgDia.solve(b);
			x = (at*x);
			std::cout << " iterations ----------> " << solverBiCgDia.iterations() << std::endl;
			std::cout << " error ---------------> " << solverBiCgDia.error() << std::endl;
		}
		void ccBiCg_augment(const std::vector<int>& type) {
			sMat d(n + AG, n + AG);
			std::vector<Tpl> coef;
			for (int p = 0; p<n; p++) {
				if (type[p] == BD2) continue;
				coef.push_back(Tpl(n, p, R(1)));
				coef.push_back(Tpl(p, n, R(1)));
			}
			d.setFromTriplets(coef.begin(), coef.end());
			a = a + d;
			b[n] = R(0.0);
			solverBiCgILU.compute(a);
			//x = solverBiCgILU.solveWithGuess(b, x);
			x = solverBiCgILU.solve(b);
			std::cout << " iterations ----------> " << solverBiCgILU.iterations() << std::endl;
			std::cout << " error ---------------> " << solverBiCgILU.error() << std::endl;
//			R rho;
//			R alpha;
//			R beta;
//			R omega;
//			R residual;
//			int loop = 0;
//			sMat pre(n + AG, n + AG);
//			std::vector<Tpl> coefList;
//			for (int p = 0; p < n; p++) {
//				coefList.push_back(Tpl(p, p, R(1.0) / a.coeff(p, p)));
//				x[p] = R(0.0);
//			}
//			x[n] = R(0.0);
//			pre.setFromTriplets(coefList.begin(), coefList.end());
//			a = pre * a;
//			b = pre * b;
//			sMat d(n + AG, n + AG);
//			std::vector<Tpl> coef;
//			for (int p = 0; p < n; p++) {
//				if (type[p] == BD2) continue;
//				coef.push_back(Tpl(n, p, 1.0));
//				coef.push_back(Tpl(p, n, 1.0));
//			}
//			d.setFromTriplets(coef.begin(), coef.end());
//			a = a + d;
//			b[n] = R(0.0);
//			vr = b - a*x;
			//residual = sqrt(vr.dot(vr));
			//if (residual < eps) {
			//	std::cout << " iterations ----------> " << loop << std::endl;
			//	std::cout << " error ---------------> " << residual << std::endl;
			//	return;
			//}
//#if OMP
//#pragma omp parallel for
//#endif
//			for (int p = 0; p < n + AG; p++) {
//				vr_hat[p] = vr[p];
//				vv[p] = vp[p] = R(0.0);
//			}
//			rho = alpha = omega = R(1.0);
//			for (loop = 0; loop < maxIter; loop++) {
//				const R rho_m1 = rho;
//				rho = vr_hat.dot(vr);
//				beta = (rho / rho_m1)*(alpha / omega);
//#if OMP
//#pragma omp parallel for
//#endif
//				for (int p = 0; p < n + AG; p++) {
//					vp[p] = vr[p] + beta*(vp[p] - omega* vv[p]);
//				}
//				vv = a* vp;
//				alpha = rho / (vr_hat.dot(vv));
//#if OMP
//#pragma omp parallel for
//#endif
//				for (int p = 0; p < n + AG; p++) {
//					vh[p] = x[p] + alpha*vp[p];
//					vs[p] = vr[p] - alpha* vv[p];
//				}
//			vt = b - a * vh;
//			residual = sqrt(vt.dot(vt));
//			if (residual < eps) {
//#if OMP
//#pragma omp parallel for
//#endif
//				for (int p = 0; p < n + AG; p++) {
//					x[p] = vh[p];
//				}
//				std::cout << " iterations ----------> " << loop << std::endl;
//				std::cout << " error ---------------> " << residual << std::endl;
//				return;
//			}
//				vt = a*vs;
//				omega = (vt.dot(vs)) / (vt.dot(vt));
//#if OMP
//#pragma omp parallel for
//#endif
//				for (int p = 0; p < n + AG; p++) {
//					x[p] = vh[p] + omega*vs[p];
//				}
//				vr = b - a*x;
//				residual = sqrt(vr.dot(vr));
//				if (residual < eps) break;
//#if OMP
//#pragma omp parallel for
//#endif
//				for (int p = 0; p < n + AG; p++) {
//					vr[p] = vs[p] - omega*vt[p];
//				}
//			}
//			std::cout << " iterations ----------> " << loop << std::endl;
//			std::cout << " error ---------------> " << residual << std::endl;
		}
		void ccCgs_horibata() {
			sMat at(n, n);
			at = a.transpose();
			dVec e(n), zero(n);
			solverBiCgDia.compute(at);
			e = solverBiCgDia.solve(0.*zero);
			const R e2 = e.dot(e);
			if (e2 > eps) {
				b = b - (b.dot(e) / e2)*e;
			}
			const dVec r0 = b - a*x;
			dVec r1 = r0;
			dVec r2 = r1;
			dVec p = r0;
			dVec e1 = b;
			dVec h = b;
			R res = r1.dot(r1);
			int k;
			for (k = 0; k < maxIter; k++) {
				res = r1.dot(r1);
				if (res < eps) break;
				const R alpha = res / r0.dot(a*p);
				h = e1 - alpha*(a*p);
				x += alpha*(e1 + h);
				r2 = r1 - alpha*(e1 + h);
				const R beta = r0.dot(r2) / r0.dot(r1);
				e1 = r2 + beta* h;
				p = e1 + beta*(h + beta*p);
				r1 = r2;
			}
			std::cout << " iterations ----------> " << k << std::endl;
			std::cout << " error ---------------> " << res << std::endl;
		}
		void lsqr() {}

	public:
		int n;
		int Dn;
		int maxIter;
		R eps;
		sMat a, au;
		dVec x, b, u, rhs;
		Eigen::BiCGSTAB< sMat, preconditionerILU > solverBiCgILU;
		Eigen::BiCGSTAB< sMat, preconditionerDia > solverBiCgDia;
		Eigen::SparseQR< sMat, Eigen::NaturalOrdering<int> > solverQR;

	private:
		void init() {
			a.resize(n + AG, n + AG), x.resize(n + AG), b.resize(n + AG);
			au.resize(D*n, D*n), u.resize(D*n), rhs.resize(D*n);
			vr.resize(n + AG), vr_hat.resize(n + AG), vv.resize(n + AG), vp.resize(n + AG), vt.resize(n + AG), vh.resize(n + AG), vs.resize(n + AG);
			Dvr.resize(D*n), Dvr_hat.resize(D*n), Dvv.resize(D*n), Dvp.resize(D*n), Dvt.resize(D*n), Dvh.resize(D*n), Dvs.resize(D*n);
			for (int i = 0; i < n; i++) {
				x[i] = b[i] = R(0.0);
			}
			for (int i = 0; i < D*n; i++) {
				u[i] = rhs[i] = (0.0);
			}
			maxIter = 1000;
			solverBiCgILU.preconditioner().setDroptol(eps);
			solverBiCgILU.preconditioner().setFillfactor(1);
			solverBiCgILU.setMaxIterations(maxIter);
			solverBiCgILU.setTolerance(eps);
			solverBiCgDia.setMaxIterations(maxIter);
			solverBiCgDia.setTolerance(eps);
			solverQR.setPivotThreshold(1.0 / n);
		}
		void fina() {}

	private:
		dVec vr, vr_hat, vv, vp, vt, vh, vs;
		dVec Dvr, Dvr_hat, Dvv, Dvp, Dvt, Dvh, Dvs;
	};

}