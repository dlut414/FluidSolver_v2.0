/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Released under CC BY-NC
*/
//Simulator.h
///defination of class Simulator
#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include "Header.h"
#include "Parameter.h"
#include "Particle.h"
#include "MatSolver.h"
#include "Sensor.h"
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif

namespace SIM {

	template <typename R, int D, typename Derived>
	class Simulator {};

	template <typename R, typename Derived>
	class Simulator<R,1,Derived> {};

	template <typename R, typename Derived>
	class Simulator<R,2,Derived> {
		typedef Eigen::Matrix<R, 2, 1> Vec;
		typedef Eigen::Matrix<R, 2, 2> Mat;
		typedef Eigen::Triplet<R> Tpl;
	public:
		Simulator() { numOfSteps = 0; }
		~Simulator() {}

		Derived& derived() { return *static_cast<Derived*>(this); }
		const Derived& derived() const { return *static_cast<const Derived*>(this); }

		void operator >> (const std::string& str) const {
			saveData(str);
		}
		void operator << (const std::string& str) {
			std::ifstream file(str);
			std::istringstream iss;
			std::string line;
			if (!file.is_open()) std::cout << " No file Para. found ! " << std::endl;
			std::getline(file, line); std::getline(file, line); iss.str(line);
			iss >> para.k; iss.clear();
			std::getline(file, line); std::getline(file, line); iss.str(line);
			iss >> para.Pr; iss.clear();
			std::getline(file, line); std::getline(file, line); iss.str(line);
			iss >> para.Ra; iss.clear();
			std::getline(file, line); std::getline(file, line); iss.str(line);
			iss >> para.cfl; iss.clear();
			std::getline(file, line); std::getline(file, line); iss.str(line);
			iss >> para.dtMax; iss.clear();
			std::getline(file, line); std::getline(file, line); iss.str(line);
			iss >> para.tt; iss.clear();
			std::getline(file, line); std::getline(file, line); iss.str(line);
			iss >> para.eps; iss.clear();
			std::getline(file, line); std::getline(file, line); iss.str(line);
			iss >> para.alpha; iss.clear();
			std::getline(file, line); std::getline(file, line); iss.str(line);
			iss >> para.beta; iss.clear();
			std::cout << " Effective radius (times of dp)   : " << para.k << std::endl;
			std::cout << " Prandtl number                   : " << para.Pr << std::endl;
			std::cout << " Rayleigh number                  : " << para.Ra << std::endl;
			std::cout << " CFL number                       : " << para.cfl << std::endl;
			std::cout << " Maximum time step (1)            : " << para.dtMax << std::endl;
			std::cout << " Total time (1)                   : " << para.tt << std::endl;
			std::cout << " EPS                              : " << para.eps << std::endl;
			std::cout << " Parameter Alpha                  : " << para.alpha << std::endl;
			std::cout << " Parameter Beta                   : " << para.beta << std::endl;
			std::cout << " Reading Para.txt done " << std::endl;
			file.close();
		}

		void init() {
			*this << "Para.txt";
			derived().init_();
			mSol = new MatSolver<R,2>(int(derived().part->np), para.eps);
			std::cout << " Particle number : " << derived().part->np << std::endl;
			R tmp = timeStep();
			para.dt = tmp < para.dtMax ? tmp : para.dtMax;
			numOfSteps = int(derived().part->ct / para.dt);
		}

		void mainLoop() {
			auto* const part = derived().part;
			while (part->ct <= para.tt) {
				std::cout << " step ----------------------------------> " << numOfSteps << std::endl;
				R tmp = timeStep();
				para.dt = tmp < para.dtMax ? tmp : para.dtMax;
				part->updateCell();
				derived().step();
				part->ct += para.dt;	numOfSteps++;
				std::cout << " time --------> " << part->ct << std::endl;
				std::cout << " dt ----------> " << para.dt << std::endl;
			}
			saveData();
		}

		R stepGL() {
			auto* const part = derived().part;
			if (part->ct > para.tt) {
				saveData();
			}
			std::cout << " step ----------------------------------> " << numOfSteps << std::endl;
			R tmp = timeStep();
			para.dt = tmp < para.dtMax ? tmp : para.dtMax;
			part->updateCell();
			derived().step();
			part->ct += para.dt;	numOfSteps++;
			std::cout << " time --------> " << part->ct << std::endl;
			std::cout << " dt ----------> " << para.dt << std::endl;
			return part->ct;
		}

		void sensorOut() {
			static int i = 0;
			std::ostringstream convert;
			convert << i++;
			*(derived().sen) >> convert.str();
		}

		void profileOut() {
			auto* const part = derived().part;
			static std::string pf = "profile";
			derived().sen->profile(rt, pf);
		}

		void saveData() const {
			static int i = 0;
			std::ostringstream convert;
			convert << i++;
			*(derived().part) >> ("./out/" + convert.str() + ".out");
		}
		void saveData(const std::string& str) const {
			*(derived().part) >> ("./out/" + str + ".out");
		}

		__forceinline const R* PositionX() const {
			return derived().part->pos[0].data();
		}
		__forceinline const R* PositionY() const {
			return derived().part->pos[1].data();
		}
		__forceinline const R* VelocityX() const {
			return derived().part->vel[0].data();
		}
		__forceinline const R* VelocityY() const {
			return derived().part->vel[1].data();
		}
		__forceinline const R* Pressure() const {
			return derived().part->pres.data();
		}
		__forceinline const R* Temperature() const {
			return derived().part->temp.data();
		}
		__forceinline const R* Divergence() const {
			return derived().part->div.data();
		}
		__forceinline const int* type() const {
			return (int*)(derived().part->type.data());
		}

	public:
		Parameter<R, 2> para;
		MatSolver<R, 2>* mSol;

	protected:
		void step() {}
		void convect() {}
		void visTerm_e() {}
		void visTerm_i() {}
		void presTerm_e() {}
		void presTerm_i() {}
		void makeDirchlet_v() {}

		void makeNeumann_p() {
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] != BD1) continue;
				part->neumann[p] = R(0);
			}
		}

		void solvMat_p() {
			mSol->biCg();
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				part->pres[p] = mSol->x[p];
				//if (part->pres[p] < -1.e5) part->pres[p] = -1.e5;
				//if (part->pres[p] > 1.e5) part->pres[p] = 1.e5;
			}
		}

		void solvMat_phi() {
			auto* const part = derived().part;
			mSol->ccBiCg_augment(part->type);
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				part->phi[p] = mSol->x[p];
				//if (part->phi[p] < -1.e5) part->phi[p] = -1.e5;
				//if (part->phi[p] > 1.e5) part->phi[p] = 1.e5;
			}
		}

		void solvMat_t() {
			mSol->biCg();
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				part->temp_m1[p] = part->temp[p];
				part->temp[p] = mSol->x[p];
			}
		}

		void solvMat_v() {
			mSol->biCg_v();
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				part->vel_p1[0][p] = mSol->u[2*p];
				part->vel_p1[1][p] = mSol->u[2*p+1];
			}
		}

		R timeStep() {
			R umax = R(0);
			const auto* const part = derived().part;
			for (int p = 0; p < part->np; p++) {
				const R ux = part->vel[0][p];
				const R uy = part->vel[1][p];
				const R tmp = sqrt(ux*ux + uy*uy);
				if (tmp > umax) umax = tmp;
			}
			para.umax = umax;
			return para.cfl * part->dp / umax;
			//const R term1 = 4 * para.Pr / (umax* umax);
			//const R term2 = para.Pr*(part->dp*part->dp) / 2;
			//const R ret = para.cfl* (term1 < term2 ? term1 : term2);
			//return ret;
		}

		void calCell() {
			derived().part->updateCell();
		}

		void calInvMat() {
			derived().part->updateInvMat();
		}

		void makeFs() {
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) part->fs[p] = part->_isFs(p);
		}

		void calForVis() {
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				part->vort[p] = part->Rot(part->vel[0].data(), part->vel[1].data(), p);
				part->div[p] = part->Div(part->vel[0].data(), part->vel[1].data(), p);
			}
		}

		void check() const {
			const auto* const part = derived().part;
			R velMax = std::numeric_limits<R>::min();
			R phiMax = std::numeric_limits<R>::min();
			R divMax = std::numeric_limits<R>::min();
			int idv = 0, idp = 0, idd = 0;
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD2) continue;
				const R ux = part->vel[0][p];
				const R uy = part->vel[1][p];
				const R vel = sqrt(ux*ux + uy*uy);
				const R phi = part->phi[p];
				const R div = part->div[p];
				if (vel > velMax) {
					velMax = vel;
					idv = p;
				}
				if (abs(phi) > abs(phiMax)) {
					phiMax = phi;
					idp = p;
				}
				if (abs(div) > abs(divMax)) {
					divMax = div;
					idd = p;
				}
			}
			std::cout << " max vel: " << velMax << " --- id: " << idv << std::endl;
			std::cout << " max phi: " << phiMax << " --- id: " << idp << std::endl;
			std::cout << " max Div: " << divMax << " --- id: " << idd << std::endl;
		}

		void insertRand() {
			auto* const part = derived().part;
			R coef = 0.25;
			std::default_random_engine gen;
			std::normal_distribution<R> dis(0., 0.5);
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				const R dr = coef* part->dp* dis(gen);
				const R theta = 2.* M_PI * (R(rand()) / RAND_MAX);
				const R dx = cos(theta)*dr;
				const R dy = sin(theta)*dr;
				part->pos[0][p] += dx;
				part->pos[1][p] += dy;
				part->pos_m1[0][p] += dx;
				part->pos_m1[1][p] += dy;
			}
		}

	protected:
		int numOfSteps;
	};

	template <typename R, typename Derived>
	class Simulator<R,3,Derived> {};
}