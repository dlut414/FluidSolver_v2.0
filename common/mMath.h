/*
*/
///math functions using metaprogramming

#pragma once

namespace mMath {

	template <unsigned M, unsigned N>	struct Power		{ enum { value = M*Power<M,N-1>::value, }; };
	template <unsigned M>				struct Power<M,0>	{ enum { value = 1, }; };

	template <typename R, unsigned N>	struct PowerR		{ static __forceinline R Run(const R& x) { return x*PowerR<R,N-1>::Run(x); } };
	template <typename R>				struct PowerR<R,0>	{ static __forceinline R Run(const R& x) { return 1; } };

	template <unsigned N>	struct Factorial	{ enum { value = N* Factorial<N-1>::value, }; };
	template <>				struct Factorial<0> { enum { value = 1, }; };

	template <unsigned N, unsigned M>
	struct H { enum { value = Factorial<N + M - 1>::value / (Factorial<M>::value * Factorial<N - 1>::value), };	};

}