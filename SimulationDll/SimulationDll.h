/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Published under CC BY-NC
*/
//SimulationDll.h
///defination of class SimulationDll
#pragma once
#ifdef SIMULATIONDLL_EXPORTS
#define SIMULATIONDLL_API __declspec(dllexport)
#else
#define SIMULATIONDLL_API __declspec(dllimport)
#endif

namespace SIM {

	typedef void* NPtr;

	class SimulationDll2D {
	public:
		static SIMULATIONDLL_API void Initialize();
		static SIMULATIONDLL_API void Run();
		static SIMULATIONDLL_API int Number();
		static SIMULATIONDLL_API NPtr Type();
		static SIMULATIONDLL_API NPtr PositionX();
		static SIMULATIONDLL_API NPtr PositionY();
		static SIMULATIONDLL_API NPtr VelocityX();
		static SIMULATIONDLL_API NPtr VelocityY();
		static SIMULATIONDLL_API NPtr Pressure();
		static SIMULATIONDLL_API NPtr Temperature();
		static SIMULATIONDLL_API NPtr Divergence();
		static SIMULATIONDLL_API NPtr Vorticity();
		static SIMULATIONDLL_API void SaveData();
		static SIMULATIONDLL_API void SensorOut();
		static SIMULATIONDLL_API void BBox(double&, double&, double&, double&);
	};

}