/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Released under CC BY-NC
*/
//Particle.h
///defination of class Particle
///2016.4.22 fixed bug in b2normal() --- 
/// normal calculation is different when initial particles position are random
#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include "LinkCell.h"
#include "Header.h"
#include "Base.h"
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif

namespace SIM {

	template <typename R, int D, typename Derived>
	class Particle : public Base<R,D> {};

	template <typename R, typename Derived>
	class Particle<R,1,Derived> : public Base<R,1> {};

	template <typename R, typename Derived>
	class Particle<R,2,Derived> : public Base<R,2> {
	public:
		typedef Eigen::Matrix<int,2,1> iVec;
		typedef Eigen::Matrix<R,2,1> Vec;
		typedef Eigen::Matrix<R,2,2> Mat;
	public:
		Particle() {}
		~Particle() {}

		Derived& derived() { return *static_cast<Derived*>(this); }
		const Derived& derived() const { return *static_cast<const Derived*>(this); }

		void clean() {
			type.clear();
			pos[0].clear(); pos[1].clear(); pos_m1[0].clear(); pos_m1[1].clear();
			vel[0].clear(); vel[1].clear(); vel_p1[0].clear(); vel_p1[1].clear(); vel_m1[0].clear(); vel_m1[1].clear();
			pres.clear();
			phi.clear(); vort.clear(); div.clear();
			bdc.clear();
			bdnorm.clear();
			p_dirichlet.clear();
			t_dirichlet.clear();
			p_neumann.clear();
			t_neumann.clear();
			np = 0;
		}
		void operator >> (const std::string str) const {
			std::ofstream file(str, std::ofstream::out);
			file << std::scientific << std::setprecision(6) << ct << std::endl;
			file << std::scientific << std::setprecision(6) << dp << std::endl;
			file << np << std::endl;
			for (int p = 0; p < np; p++) {
				file << std::scientific << std::setprecision(6);
				file << type[p] << " " << pos[0][p] << " " << pos[1][p] << " " << vel[0][p] << " " << vel[1][p] << " " << temp[p];
				if (type[p] == BD1 || type[p] == INLET || type[p] == OUTLET) file << " " << bdnorm.at(p)[0] << " " << bdnorm.at(p)[1];
				file << std::endl;
			}
			std::cout << " Writing Geo.in done. " << std::endl;
			file.close();
		}
		void operator << (const std::string str) {
			int n; int t; Vec p; Vec v; R tp; Vec norm;
			std::ifstream file(str);
			if (!file.is_open()) std::cout << " File Geo.in not found ! " << std::endl;
			file >> ct >> dp >> n;
			while (n-- > 0) {
				file >> t;
				file >> p[0] >> p[1];
				file >> v[0] >> v[1];
				file >> tp;
				if (t == BD1 || t == INLET || t == OUTLET) {
					file >> norm[0] >> norm[1];
					addPart(pType(t), p, v, tp, norm);
				}
				else addPart(pType(t), p, v, tp);
			}
			file.close();
			std::cout << " Reading Geo.in done " << std::endl;
		}

		void addPart(const pType& t, const Vec& p, const Vec& v, const R& tp) {
			type.push_back(t);
			pos[0].push_back(p[0]);	pos[1].push_back(p[1]); pos_m1[0].push_back(p[0]); pos_m1[1].push_back(p[1]);
			vel[0].push_back(v[0]); vel[1].push_back(v[1]); vel_p1[0].push_back(v[0]); vel_p1[1].push_back(v[1]); vel_m1[0].push_back(v[0]); vel_m1[1].push_back(v[1]);
			temp.push_back(tp); temp_m1.push_back(tp); pres.push_back(R(0)); phi.push_back(R(0)); vort.push_back(R(0)); div.push_back(R(0));
			bdc.push_back(0);
			bdnorm.push_back(Vec());
			p_dirichlet.push_back(0);
			t_dirichlet.push_back(0);
			p_neumann.push_back(0);
			t_neumann.push_back(0);
			np++;
		}
		void addPart(const pType& t, const Vec& p, const Vec& v, const R& tp, const Vec& norm) {
			type.push_back(t);
			pos[0].push_back(p[0]);	pos[1].push_back(p[1]); pos_m1[0].push_back(p[0]); pos_m1[1].push_back(p[1]);
			vel[0].push_back(v[0]); vel[1].push_back(v[1]); vel_p1[0].push_back(v[0]); vel_p1[1].push_back(v[1]); vel_m1[0].push_back(v[0]); vel_m1[1].push_back(v[1]);
			temp.push_back(tp); temp_m1.push_back(tp); pres.push_back(R(0)); phi.push_back(R(0)); vort.push_back(R(0)); div.push_back(R(0));
			bdc.push_back(0);
			bdnorm.push_back(norm);
			p_dirichlet.push_back(0);
			t_dirichlet.push_back(0);
			p_neumann.push_back(0);
			t_neumann.push_back(0);
			np++;
		}
		void erasePart(const int& offset) {
			type.erase(type.begin() + offset);
			pos[0].erase(pos[0].begin() + offset); pos[1].erase(pos[1].begin() + offset); pos_m1[0].erase(pos_m1[0].begin() + offset); pos_m1[1].erase(pos_m1[1].begin() + offset);
			vel[0].erase(vel[0].begin() + offset); vel[1].erase(vel[1].begin() + offset); vel_p1[0].erase(vel_p1[0].begin() + offset); vel_p1[1].erase(vel_p1[1].begin() + offset); vel_m1[0].erase(vel_m1[0].begin() + offset); vel_m1[1].erase(vel_m1[1].begin() + offset);
			temp.erase(temp.begin() + offset); temp_m1.erase(temp_m1.begin() + offset); pres.erase(pres.begin() + offset); phi.erase(phi.begin() + offset); vort.erase(vort.begin() + offset); div.erase(div.begin() + offset);
			bdc.erase(bdc.begin() + offset);
			bdnorm.erase(bdnorm.begin() + offset);
			p_dirichlet.erase(p_dirichlet.begin() + offset);
			t_dirichlet.erase(t_dirichlet.begin() + offset);
			p_neumann.erase(p_neumann.begin() + offset);
			t_neumann.erase(t_neumann.begin() + offset);
			np--;
		}

		void buildCell() {
			BBox<R> b = BBox<R>();
			for (int p = 0; p < np; p++) {
				Vec pp;
				pp[0] = pos[0][p];
				pp[1] = pos[1][p];
				b += pp;
			}
			b.Expand(0.1);
			cell = new LinkCell<R,2>(b, r0);
			updateCell();
		}
		void updateCell() {
			cell->update(pos[0].data(), pos[1].data(), np);
		}
		void getBBox(R& left, R& right, R& bottom, R& top) const {
			cell->getBBox(left, right, bottom, top);
		}

		void makeBdc() {
			for (int p = 0; p < np; p++) {
				bdc[p] = 0;
				if (type[p] == BD1 || type[p] == INLET || type[p] == OUTLET) {
					bdc[p] = ON(bdc[p], P_NEUMANN);
				}
			}
		}

		void b2neumann() {
			for (int p = 0; p < np; p++) {
				if (IS(bdc[p], P_NEUMANN)) p_neumann[p] = R(0);
				if (IS(bdc[p], T_NEUMANN)) t_neumann[p] = R(0);
			}
		}
		void b2dirichlet() {
			for (int p = 0; p < np; p++) {
				if (IS(bdc[p], P_DIRICHLET)) p_dirichlet[p] = R(0);
				if (IS(bdc[p], T_DIRICHLET0)) t_dirichlet[p] = R(0);
				if (IS(bdc[p], T_DIRICHLET1)) t_dirichlet[p] = R(1);
			}
		}

	public:
		R ct;
		int np;
		std::vector<R> pos[2];
		std::vector<R> pos_m1[2];
		std::vector<R> vel[2];
		std::vector<R> vel_p1[2];
		std::vector<R> vel_m1[2];

		std::vector<R> temp;
		std::vector<R> temp_m1;
		std::vector<R> pres;
		std::vector<R> div;
		std::vector<pType> type;
		std::vector<int> bdc;
		std::vector<R> phi;
		std::vector<R> vort;
		std::vector<Vec> bdnorm;
		std::vector<R> p_dirichlet;
		std::vector<R> t_dirichlet;
		std::vector<R> p_neumann;
		std::vector<R> t_neumann;

		LinkCell<R,2>* cell;
	};

	template <typename R, typename Derived>
	class Particle<R,3,Derived> : public Base<R,3> {};

}