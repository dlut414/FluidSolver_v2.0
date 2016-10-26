/*
*/
#pragma once
#include <vector>
#include <Eigen/Dense>
#include "Header.h"
#define _USE_MATH_DEFINES
#include <cmath>

namespace SIM {

	template <typename R, int D>
	class Base {
		typedef Eigen::Matrix<R, D, 1> vec;
		typedef Eigen::Matrix<R, 2, 1> vec2;
		typedef Eigen::Matrix<R, 3, 1> vec3;
		typedef Eigen::Matrix<R, 5, 1> vec5;
		typedef Eigen::Matrix<R, 6, 1> vec6;
		typedef Eigen::Matrix<R, 7, 1> vec7;
		typedef Eigen::Matrix<R, 8, 1> vec8;
		typedef Eigen::Matrix<R, 9, 1> vec9;
	public:
		Base() {}
		~Base() {}

		inline const R ww(const R& r) const {
			if (r >= r0) {
				return 0.;
			}
			else {
				return pow((1 - r / r0), 2);
			}
		}

	public:
		R dp;
		R k;
		R r0;
		R eps;
		R eps_mat;

	public:
		void init(const R& _k) {
			k = _k; r0 = k* dp;
			eps = static_cast<R>(1.e-6);
			eps_mat = static_cast<R>(1.e-6);
		}

	};

}