/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Published under CC BY-NC
*/
//VisualizationDll.h
///defination of class VisualizationDll
#pragma once
#ifdef VISUALIZATIONDLL_EXPORTS
#define VISUALIZATIONDLL_API __declspec(dllexport)
#else
#define VISUALIZATIONDLL_API __declspec(dllimport)
#endif
#include <Controller.h>

namespace VIS {

	typedef void* NPtr;

	class VisualizationDll {
	public:
		static VISUALIZATIONDLL_API void Initialize();
		static VISUALIZATIONDLL_API void Run(const Controller* const controlPtr, const int& dim, const int& num, NPtr tp, NPtr pos, NPtr s);
		static VISUALIZATIONDLL_API void Run(const Controller* const controlPtr, const int& num, NPtr tp, NPtr posX, NPtr posY, NPtr s);
		static VISUALIZATIONDLL_API int IntersectColorPick(const Controller* const controlPtr, const int& num, const GLuint& mouseX, const GLuint& mouseY);
		static VISUALIZATIONDLL_API void Finalize();
	};

}