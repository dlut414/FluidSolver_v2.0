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
				if (type[p] == INLET) file << " " << pos_m1[0][p] << " " << pos_m1[1][p];
				file << std::endl;
			}
			std::cout << " Writing Geo.in done. " << std::endl;
			file.close();
		}
		void operator << (const std::string str) {
			int n; int t; Vec p; Vec v; R tp; Vec norm; Vec p_m1;
			std::ifstream file(str);
			if (!file.is_open()) std::cout << " File Geo.in not found ! " << std::endl;
			file >> ct >> dp >> n;
			while (n-- > 0) {
				file >> t;
				file >> p[0] >> p[1];
				file >> v[0] >> v[1];
				file >> tp;
				if (t == BD1 || t == OUTLET) {
					file >> norm[0] >> norm[1];
					addPart(int(t), p, v, tp, norm);
				}
				else if (t == INLET) {
					file >> norm[0] >> norm[1] >> p_m1[0] >> p_m1[1];
					addPart(int(t), p, v, tp, norm, p_m1);
				}
				else addPart(int(t), p, v, tp);
			}
			file.close();
			std::cout << " Reading Geo.in done " << std::endl;
		}

		void addPart(const int& t, const Vec& p, const Vec& v, const R& tp) {
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
		void addPart(const int& t, const Vec& p, const Vec& v, const R& tp, const Vec& norm) {
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
		void addPart(const int& t, const Vec& p, const Vec& v, const R& tp, const Vec& norm, const Vec& p_m1) {
			type.push_back(t);
			pos[0].push_back(p[0]);	pos[1].push_back(p[1]); pos_m1[0].push_back(p_m1[0]); pos_m1[1].push_back(p_m1[1]);
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
		void swapPart(const int& p, const int& q) {
			int t = type[p]; type[p] = type[q], type[q] = t;
			R posx = pos[0][p]; pos[0][p] = pos[0][q], pos[0][q] = posx;
			R posy = pos[1][p]; pos[1][p] = pos[1][q], pos[1][q] = posy;
			R posx_m1 = pos_m1[0][p]; pos_m1[0][p] = pos_m1[0][q], pos_m1[0][q] = posx_m1;
			R posy_m1 = pos_m1[1][p]; pos_m1[1][p] = pos_m1[1][q], pos_m1[1][q] = posy_m1;
			R velx = vel[0][p]; vel[0][p] = vel[0][q], vel[0][q] = velx;
			R vely = vel[1][p]; vel[1][p] = vel[1][q], vel[1][q] = vely;
			R velx_p1 = vel_p1[0][p]; vel_p1[0][p] = vel_p1[0][q], vel_p1[0][q] = velx_p1;
			R vely_p1 = vel_p1[1][p]; vel_p1[1][p] = vel_p1[1][q], vel_p1[1][q] = vely_p1;
			R velx_m1 = vel_m1[0][p]; vel_m1[0][p] = vel_m1[0][q], vel_m1[0][q] = velx_m1;
			R vely_m1 = vel_m1[1][p]; vel_m1[1][p] = vel_m1[1][q], vel_m1[1][q] = vely_m1;
			R tmp = temp[p]; temp[p] = temp[q], temp[q] = tmp;
			R tmp_m1 = temp_m1[p]; temp_m1[p] = temp_m1[q], temp_m1[q] = tmp_m1;
			R ps = pres[p]; pres[p] = pres[q], pres[q] = ps;
			R ph = phi[p]; phi[p] = phi[q], phi[q] = ph;
			R vt = vort[p]; vort[p] = vort[q], vort[q] = vt;
			R dv = div[p]; div[p] = div[q], div[q] = dv;
			int bc = bdc[p]; bdc[p] = bdc[q], bdc[q] = bc;
			Vec bm = bdnorm[p]; bdnorm[p] = bdnorm[q], bdnorm[q] = bm;
			R pd = p_dirichlet[p]; p_dirichlet[p] = p_dirichlet[q], p_dirichlet[q] = pd;
			R td = t_dirichlet[p]; t_dirichlet[p] = t_dirichlet[q], t_dirichlet[q] = td;
			R pn = p_neumann[p]; p_neumann[p] = p_neumann[q], p_neumann[q] = pn;
			R tn = t_neumann[p]; t_neumann[p] = t_neumann[q], t_neumann[q] = tn;
		}
		void copyPart_phisicalValue(const int& p, const int& q) {
			pos[0][p] = pos[0][q];
			pos[1][p] = pos[1][q];
			pos_m1[0][p] = pos_m1[0][q];
			pos_m1[1][p] = pos_m1[1][q];
			vel[0][p] = vel[0][q];
			vel[1][p] = vel[1][q];
			vel_p1[0][p] = vel_p1[0][q];
			vel_p1[1][p] = vel_p1[1][q];
			vel_m1[0][p] = vel_m1[0][q];
			vel_m1[1][p] = vel_m1[1][q];
			temp[p] = temp[q];
			temp_m1[p] = temp_m1[q];
			pres[p] = pres[q];
			phi[p] = phi[q];
			vort[p] = vort[q];
			div[p] = div[q];
		}

		void buildCell() {
			BBox<R> b = BBox<R>();
			for (int p = 0; p < np; p++) {
				Vec pp;
				pp[0] = pos[0][p];
				pp[1] = pos[1][p];
				b += pp;
			}
			b.Expand(0.05);
			cell = new LinkCell<R,2>(b, r0);
			updateCell();
		}
		void updateCell() {
			cell->update(pos[0].data(), pos[1].data(), np);
		}
		void getBBox(R& left, R& right, R& bottom, R& top) const {
			cell->getBBox(left, right, bottom, top);
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
		std::vector<int> type;
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