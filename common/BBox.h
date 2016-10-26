/*
 * LICENCE
 * copyright 2014 ~ ****
 * Some rights reserved.
 * Author: HUFANGYUAN
 * Released under CC BY-NC
*/
#ifndef BBOX_H
#define BBOX_H

#include <limits>
#include <algorithm>
#include <Eigen/Dense>
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif

template <typename R, int D = 2>
class BBox {
	typedef Eigen::Matrix<R,D,1> vec;
public:
    BBox() {
        R _infinity = std::numeric_limits<R>::infinity();
		for (int i = 0; i < D; i++) {
			pMin[i] = _infinity;
			pMax[i] = -_infinity;
		}
    }

	template <typename T>
	BBox(const BBox<T, D>& b) {
		for (int i = 0; i < D; i++) {
			pMin[i] = static_cast<R>(b.pMin[i]);
			pMax[i] = static_cast<R>(b.pMax[i]);
		}
	}

    BBox(const vec& p) {
		for (int i = 0; i < D; i++) {
			pMin[i] = p[i];
			pMax[i] = p[i];
		}
	}

    BBox(const vec& p1, const vec& p2) {
		for (int i = 0; i < D; i++) {
			pMin[i] = std::min(p1[i], p2[i]);
			pMax[i] = std::max(p1[i], p2[i]);
		}
    }

	void operator+=(const vec& p) {
		for (int i = 0; i < D; i++) {
			pMin[i] = std::min(pMin[i], p[i]);
			pMax[i] = std::max(pMax[i], p[i]);
		}
	}

    BBox Union(const BBox& b1, const BBox &b2) const {
		BBox ret = b1;
		for (int i = 0; i < D; i++) {
			ret.pMin[i] = std::min(b1.pMin[i], b2.pMin[i]);
			ret.pMax[i] = std::max(b1.pMax[i], b2.pMax[i]);
		}
        return ret;
    }

    bool bOverlaps(const BBox &b) const {
		bool ret = true;
		for (int i = 0; i < D; i++) {
			ret = ret && (pMax[i] >= b.pMin[i]) && (pMin[i] <= b.pMax[i]);
		}
        return ret;
    }

    bool bInside(const vec &p) const {
		bool ret = true;
		for (int i = 0; i < D; i++) {
			ret = ret && (p[i] >= pMin[i] && p[i] <= pMax[i]);
		}
		return ret;
    }

    vec pCenter() const {
		vec ret;;
		for (int i = 0; i < D; i++) {
			ret[i] = 0.5* (pMin[i] + pMax[i]);
		}
		return ret;
    }

    void Expand(const R& s) {
        vec v = s * (pMax - pMin);
        pMax += v;
        pMin -= v;
    }

public:
    vec pMin, pMax;

};

#endif // BBOX_H
