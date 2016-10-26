/*
*/
#pragma once
#include <iostream>
#include <vector>
#include <unordered_map>
#include <cassert>
#include <Eigen/Dense>
#include "mMath.h"
#include <BBox.h>

namespace SIM {

	template <typename R, unsigned D>
	class LinkCell {
		typedef Eigen::Matrix<R,D,1>	VecD;
		typedef Eigen::Matrix<int,D,1>	iVecD;
	public:
		template <typename T>
		LinkCell(const BBox<T,D>& b, const R& s) : box(b), cSize(s) { init(); }
		~LinkCell() { fina(); }

		template <unsigned D_ = 0, bool Over = (D_==(D-1))> struct Convert {
			template <typename U, typename V> static __forceinline void Gen(const U* const in, V* const out)				{ out[D_] = static_cast<V>(in[D_]); Convert<D_+1>::Gen(in,out); } 
			template <typename U, typename V> static __forceinline void Gen(const R& s, const U* const in, V* const out)	{ out[D_] = static_cast<V>(s*in[D_]); Convert<D_+1>::Gen(s, in, out); }
		};
		template <unsigned D_>								struct Convert<D_,true> { 
			template <typename U, typename V> static __forceinline void Gen(const U* const in, V* const out)				{ out[D_] = static_cast<V>(in[D_]); } 
			template <typename U, typename V> static __forceinline void Gen(const R& s, const U* const in, V* const out)	{ out[D_] = static_cast<V>(s*in[D_]); }
		};

		template <unsigned D_ = D>	__forceinline const unsigned hash		(const iVecD& c)	const {}
		template <>					__forceinline const unsigned hash<1>	(const iVecD& c)	const { return unsigned((c[0]) % cNum); }
		template <>					__forceinline const unsigned hash<2>	(const iVecD& c)	const { return unsigned((c[0] + dv[0]* c[1]) % cNum); }
		template <>					__forceinline const unsigned hash<3>	(const iVecD& c)	const { return unsigned((c[0] + dv[0]* c[1] + sheet* c[2] + 100000) % cNum); }
		
		__forceinline const iVecD iCoord(const VecD& p) const { static R cSizeInv = 1. / cSize; iVecD ret; Convert<>::Gen(cSizeInv, p.data(), ret.data()); return ret; }
		__forceinline const unsigned hash(const iVecD& c) const { return hash<>(c); }
		__forceinline const unsigned hash_(const VecD& p) const { return hash<>(iCoord(p)); }
		
		struct blockSize { enum { value = mMath::Power<3,D>::value, }; };
		
		__forceinline const unsigned hash(const iVecD& c, const unsigned& i) const { return hash(c+loopTable.at(i)); }

		void update(const std::vector<VecD>& pos) {
			for (unsigned i = 0; i < linkList.size(); i++) {
				linkList[i].clear();
			}
			for (unsigned i = 0; i < pos.size(); i++) {
				linkList[hash_(pos[i])].push_back(i);
			}
		}

	public:
		BBox<R,D> box;
		std::vector< std::vector<unsigned> > linkList;

	private:
		void init() {
			VecD v = VecD(box.pMax - box.pMin);
			for (int d = 0; d < int(D); d++) { dv[d] = int(v[d] / cSize) + 1; dv[d] = dv[d] > 3 ? dv[d] : 3; }
			cNum = 1;
			for (int d = 0; d < int(D); d++) { cNum *= unsigned(dv[d]); }
			if(D==3) sheet = unsigned(dv[0] * dv[1]);
			else sheet = 1;

			if (cNum >= linkList.max_size()) {
				std::cout << " linkList overflow ! " << std::endl;
				exit(0);
			}
			linkList.clear();
			linkList = std::vector< std::vector<unsigned> >(cNum);
			std::cout << " Cell num: " << cNum << std::endl;

			loopTable.clear();
			switch (D) {
			case 1:
				for (auto i = 0; i < 3; i++) {
					loopTable[i] << (i - 1);
				}
				break;
			case 2:
				for (auto i = 0; i < 9; i++) {
					loopTable[i] << (i % 3 - 1), (i / 3 - 1);
				}
				break;
			case 3:
				for (auto i = 0; i < 27; i++) {
					const unsigned z = i / 9;
					const unsigned y = (i % 9) / 3;
					const unsigned x = (i % 9) % 3;
					loopTable[i] << (x - 1), (y - 1), (z - 1);
				}
				break;
			default:
				break;
			}
		}
		void fina() { linkList.clear(); }

	private:
		iVecD dv;
		unsigned sheet;
		unsigned cNum;
		R cSize;
		std::unordered_map<unsigned, iVecD> loopTable;
	};

}