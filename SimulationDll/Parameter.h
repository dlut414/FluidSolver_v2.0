/*
*/
#pragma once
#include "Header.h"

namespace SIM {

	template <typename R, unsigned D>
	class Parameter {
	public:
		Parameter() { init(); }
		~Parameter() {}

	public:
		R k;		//radius k*dp
		R eps;		//eps

		R Pr;		//Prandtl number
		R Ra;		//Rayleigh number

		R cfl;		//cfl num
		R dtMax;	//max time step
		R tt;		//total time
		R dt;		//current time step
		R umax;		//current umax

		R alpha;	//arbitrary parameter
		R beta;		//arbitrary parameter
	private:
		void init() {
			dt = 1.0;
		}
	};

}