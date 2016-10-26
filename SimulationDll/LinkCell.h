/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Released under CC BY-NC
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

	template <typename R, int D>
	class LinkCell {};

	template <typename R>
	class LinkCell<R,2> {
		typedef Eigen::Matrix<int,2,1>	iVecD;
	public:
		template <typename T>
		LinkCell(const BBox<T,2>& b, const R& s) : box(b), cSize(s), cSizeInv(R(1)/s) { initialize(); }
		~LinkCell() { finalize(); }

		struct blockSize { enum { value = mMath::Power<3,2>::value, }; };
		
		__forceinline const int pos2cell(const R& x) const { return int(cSizeInv*x); }
		__forceinline const int hash_(const int& cx, const int& cy) const { return int((cx + nocX * cy + cNum) % cNum); }
		__forceinline const int hash(const R& px, const R& py) const { return hash_(pos2cell(px), pos2cell(py)); }
		__forceinline const int hash(const int& cx, const int& cy, const int& i) const { return hash_(cx + loopTable.at(i)[0], cy + loopTable.at(i)[1]); }

		void update(const R* const px, const R* const py, const int& np) {
			for (int i = 0; i < linkList.size(); i++) {
				linkList[i].clear();
			}
			for (int i = 0; i < np; i++) {
				linkList[hash(px[i], py[i])].push_back(i);
			}
		}

	public:
		BBox<R,2> box;
		std::vector<std::vector<int>> linkList;

	private:
		void initialize() {
			R vx = R(box.pMax[0] - box.pMin[0]);
			R vy = R(box.pMax[1] - box.pMin[1]);
			nocX = int(vx / cSize) + 1; nocX = nocX > 3 ? nocX : 3;
			nocY = int(vy / cSize) + 1; nocY = nocY > 3 ? nocY : 3;
			cNum = nocX * nocY;

			if (cNum >= linkList.max_size()) {
				std::cout << " linkList overflow ! " << std::endl;
				exit(0);
			}
			linkList.clear();
			linkList = std::vector< std::vector<int> >(cNum);
			std::cout << " Cell num : " << cNum << std::endl;

			loopTable.clear();
			for (int i = 0; i < blockSize::value; i++) {
				loopTable[i] << (i % 3 - 1), (i / 3 - 1);
			}
		}
		void finalize() { linkList.clear(); }

	private:
		int nocX;
		int nocY;
		int cNum;
		R cSize;
		R cSizeInv;
		std::unordered_map<int, iVecD> loopTable;
	};

	template <typename R>
	class LinkCell<R,3> {};

}