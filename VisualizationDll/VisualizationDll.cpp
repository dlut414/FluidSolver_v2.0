// Visualization.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "VisualizationDll.h"
#include "Header.h"
#include "DrawParticle.h"
#include "DrawStreamline.h"
#include <PreInformation.h>

namespace VIS {
	
	typedef Parameters::DataType Type;
	typedef DrawParticle<Type> DP;
	typedef DrawParticle<Type>* DPPtr;

	DrawParticle<Parameters::DataType>* drawer_particle;
	DrawStreamline<Type>* drawer_stream;

	void VisualizationDll::Initialize() {
		glEnable(GL_TEXTURE_1D);
		glEnable(GL_TEXTURE_2D);
		glEnable(GL_TEXTURE_3D);
		glEnable(GL_CULL_FACE);
		//glDisable(GL_CULL_FACE);
		glFrontFace(GL_CCW);
		glEnable(GL_POINT_SPRITE_ARB);
		glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
		glEnable(GL_DEPTH_TEST);
		glDepthFunc(GL_LESS);
		glEnable(GL_ALPHA_TEST);
		glAlphaFunc(GL_GREATER, 0.f);
		//glEnable(GL_BLEND);
		//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		//glEnable(GL_LINE_SMOOTH);
		glLineWidth(1.f);
		glClearColor(1.f, 1.f, 1.f, 0.f);
		glewInit();
		drawer_particle = new DP();
		drawer_stream = new DrawStreamline<Type>;
	}

	void VisualizationDll::Run(const Controller* const controlPtr, const int& dim, const int& num, NPtr tp, NPtr pos, NPtr s) {
		drawer_particle->Draw(controlPtr, dim, num, tp, pos, s);
	}
	void VisualizationDll::Run(const Controller* const controlPtr, const int& num, NPtr tp, NPtr posX, NPtr posY, NPtr s) {
		drawer_particle->Draw(controlPtr, num, tp, posX, posY, s);
	}

	int VisualizationDll::IntersectColorPick(const Controller* const controlPtr, const int& num, const GLuint& mouseX, const GLuint& mouseY) {
		return drawer_particle->IntersectColorPick(controlPtr, num, mouseX, mouseY);
	}
	void VisualizationDll::Run_stream(const Controller* const controlPtr, const int refreshBuffer, void(*interp)(const Type, const Type, Type&, Type&)) {
		if (refreshBuffer) drawer_stream->makeStreamline(controlPtr, interp);
		drawer_stream->Draw(controlPtr);
	}

	void VisualizationDll::Finalize() {

	}

}


