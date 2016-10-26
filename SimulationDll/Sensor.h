/*
*/
#pragma once
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include "Particle.h"

namespace SIM {

	template <typename R, unsigned D, typename Der>
	class Sensor {};

	template <typename R, typename Der>
	class Sensor<R,1,Der> {};

	template <typename R, typename Der>
	class Sensor<R,2,Der> {
		typedef Eigen::Matrix<R,2,1> Vec;
	public:
		Sensor(Particle<R,2,Der>* _ptr) : ptr(_ptr) {
			pos.clear();	vel.clear();	pres.clear(); gd.clear();
			const auto& pt = ptr->derived();
			for (int p = 0; p < pt.np; p++) {
				if (IS(pt.bdc[p], T_DIRICHLET1)) {
					id.push_back(p); gd.push_back(Vec::Zero());
				}
			}
			for (int p = 0; p < pt.np; p++) {
				if (IS(pt.bdc[p], T_DIRICHLET0)) {
					id.push_back(p); gd.push_back(Vec::Zero());
				}
			}
		}
		~Sensor() {}

		void operator >> (const std::string& str) {
			interpolateVelocity();
			std::ofstream file("./out/s" + str + ".out", std::ofstream::out);
			for (auto s = 0; s < pos.size(); s++) {
				for (auto d = 0; d < 2; d++) {
					file << std::setprecision(6) << std::scientific << pos[s][d] << " ";
				}
				for (auto d = 0; d < 2; d++) {
					file << std::setprecision(6) << std::scientific << vel[s][d] << " ";
				}
				file << std::endl;
			}
			file.close();
			computeNusselt();
			const auto& pt = ptr->derived();
			file.open("./out/nu" + str + ".out", std::ofstream::out);
			for (auto s = 0; s < id.size(); s++) {
				auto p = id[s];
				for (auto d = 0; d < 2; d++) {
					file << std::setprecision(6) << std::scientific << pt.pos[d][p] << " ";
				}
				for (auto d = 0; d < 2; d++) {
					file << std::setprecision(6) << std::scientific << gd[s][d] << " ";
				}
				file << std::endl;
			}
			file.close();
			std::cout << " Writing Sensor. done " << std::endl;
		}
		void operator << (const std::string& str) {
			std::ifstream file(str);
			if (!file.is_open()) std::cout << " No file Sensor. file ! " << std::endl;
			while (file.good()) {
				Vec censorPos;
				file >> censorPos[0] >> censorPos[1];
				pos.push_back(censorPos);
				vel.push_back(Vec::Zero());
				pres.push_back(0.0);
			}
			file.close();
			std::cout << " Reading Sensor. done " << std::endl;
		}

		void profile(const R& time, const std::string& str) const {
			interpolatePressure();
			std::ofstream file("./out/s" + str + ".out", std::ofstream::app);
			for (unsigned s = 0; s < pos.size(); s++) {
				file << std::setprecision(6) << std::scientific << time << " "
					<< std::setprecision(6) << std::scientific << pres[s]
					<< std::endl;
			}
			file.close();
			std::cout << " Writing profile. done " << std::endl;
		}

		void interpolateVelocity() {
			const auto& pt = ptr->derived();
			for (auto s = 0; s < pos.size(); s++) {
				vel[s] = pt.interpolateLSA(pt.vel[0].data(), pt.vel[1].data(), pos[s][0], pos[s][1]);
			}
		}
		void interpolatePressure() {
			const auto& pt = ptr->derived();
			for (unsigned s = 0; s < pos.size(); s++) {
				pres[s] = pt.interpolateLSA(pt.pres.data(), pos[s][0], pos[s][1]);
			}
		}
		void computeNusselt() {
			const auto& pt = ptr->derived();
			for (unsigned s = 0; s < id.size(); s++) {
				gd[s] = pt.Grad(pt.temp.data(), id[s]);
			}
		}

	public:
		std::vector<int> id;
		std::vector<Vec> pos;
		std::vector<Vec> vel;
		std::vector<Vec> gd;
		std::vector<R> pres;

	private:
		Particle<R,2,Der>* ptr;
	};

	template <typename R, typename Der>
	class Sensor<R,3,Der> {};
}