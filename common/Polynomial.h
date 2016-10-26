/*
* LICENSE
*/
///Polynomial.h
#pragma once
#include <vector>
#include "mMath.h"

namespace mMath {

	template <typename R, unsigned D, unsigned P>
	class Polynomial_A_ {
	public:
		enum { value = H<D, P>::value + Polynomial_A_<R,D,P-1>::value, };
		template <unsigned N>	static __forceinline R Power(const R& x) { return x*Power<N - 1>(x); }
		template <>				static __forceinline R Power<0>(const R& x) { return 1; }
	};
	template <typename R, unsigned D>	class Polynomial_A_<R,D,0> { public: enum { value = 0, }; };

	template <typename R, unsigned D, unsigned P>	class Polynomial_A : public Polynomial_A_<R,D,P>{};

	template <typename R, unsigned P>				class Polynomial_A<R,1,P> : public Polynomial_A_<R,1,P> {
	private:
		template <unsigned P_ = 1>	static __forceinline void Run_		(const R& x, R* const out) { out[P_-1] = Power<P_>(x); Run_<P_+1>(x, out); }
		template <>					static __forceinline void Run_<P>	(const R& x, R* const out) { out[P-1] = Power<P>(x); }
		template <unsigned P_ = 1>	static __forceinline void Run_		(const R& s, const R& x, R* const out) { out[P_-1] = Power<P_>(s*x); Run_<P_+1>(s, x, out); }
		template <>					static __forceinline void Run_<P>	(const R& s, const R& x, R* const out) { out[P-1] = Power<P>(s*x); }
	public:
		static __forceinline void Run(const R& x, R* const out)				{ Run_(x, out); }
		static __forceinline void Run(const R& s, const R& x, R* const out) { Run_(s, x, out); }
	};

	template <typename R, unsigned P>				class Polynomial_A<R,2,P> : public Polynomial_A_<R,2,P> {
	private:
		template <unsigned P_ = 1, unsigned PX = P_, unsigned I = 0>
		struct Run__			{ 
			static __forceinline void Run(const R* const x, R* const out)				{ out[I] = Power<PX>(x[0])*Power<P_-PX>(x[1]); Run__<P_,PX-1,I+1>::Run(x, out); } 
			static __forceinline void Run(const R& s, const R* const x, R* const out)	{ out[I] = Power<PX>(s*x[0])*Power<P_-PX>(s*x[1]); Run__<P_,PX-1,I+1>::Run(s, x, out); }
		};
		template <unsigned P_, unsigned I>						
		struct Run__<P_,0,I>	{ 
			static __forceinline void Run(const R* const x, R* const out)				{ out[I] = Power<P_>(x[1]); }
			static __forceinline void Run(const R& s, const R* const x, R* const out)	{ out[I] = Power<P_>(s*x[1]); }
		};
		template <unsigned P_ = 1, unsigned I = 0, bool Over = (P==P_)>
		struct Run_				{ 
			static __forceinline void Run(const R* const x, R* const out)				{ Run__<P_,P_,I>::Run(x, out); Run_<P_+1,I+H<2,P_>::value>::Run(x, out); } 
			static __forceinline void Run(const R& s, const R* const x, R* const out)	{ Run__<P_,P_,I>::Run(s, x, out); Run_<P_+1,I+H<2,P_>::value>::Run(s, x, out); }
		};
		template <unsigned P_, unsigned I>
		struct Run_<P_,I,true>	{ 
			static __forceinline void Run(const R* const x, R* const out)				{ Run__<P_,P_,I>::Run(x, out); } 
			static __forceinline void Run(const R& s, const R* const x, R* const out)	{ Run__<P_,P_,I>::Run(s, x, out); }
		};
	public:
		static __forceinline void Run(const R* const x, R* const out)				{ Run_<>::Run(x, out); }
		static __forceinline void Run(const R& s, const R* const x, R* const out)	{ Run_<>::Run(s, x, out); }
	};
	


	template <typename R, unsigned D, unsigned P>
	class Polynomial_B_ {
	public:
		enum { value = H<D, P>::value + Polynomial_B_<R, D, P-1>::value, };
	};
	template <typename R, unsigned D>
	class Polynomial_B_ < R, D, 0 > {
	public:
		enum { value = 1, };
	};

}