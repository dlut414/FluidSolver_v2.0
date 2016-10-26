// Manipulator.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Released under CC BY-NC
*/
///Manipulator.cpp main loop
#pragma once
#include <Controller.h>
#include "Header.h"
#include "Bitmap.h"
#include "../VisualizationDll/VisualizationDll.h"
#include "../SimulationDll/SimulationDll.h"
#include <PreInformation.h>

typedef VIS::VisualizationDll Visualization;
typedef SIM::SimulationDll2D Simulation;

static VIS::Controller control;
static TwBar* GUIBar;

static void Render() {
	//Visualization::Run(&control, Parameters::Dimension, Simulation::Number(), Simulation::Type(), Simulation::Position(), Simulation::Scalar());
	switch (control.m_mode) {
	case VIS::DMODE_ONE:
		Visualization::Run(&control, Simulation::Number(), Simulation::Type(), Simulation::PositionX(), Simulation::PositionY(), Simulation::Temperature());
		break;
	case VIS::DMODE_TWO:
		Visualization::Run(&control, Simulation::Number(), Simulation::Type(), Simulation::PositionX(), Simulation::PositionY(), Simulation::Divergence());
		break;
	case VIS::DMODE_THREE:
		Visualization::Run(&control, Simulation::Number(), Simulation::Type(), Simulation::PositionX(), Simulation::PositionY(), Simulation::Pressure());
		break;
	default:
		break;
	}
}

static void setTwVisible(TwBar* const bar, const int visible) {
	TwSetParam(GUIBar, NULL, "visible", TW_PARAM_INT32, 1, &visible);
	Render();
	glutSwapBuffers();
}

static void callBack() {
	static int outSwitchS = 0;
	static int outSwitchP = 0;
	static int count = 0;
	if (control.i_save) {
		Simulation::SaveData();
		control.i_save = 0;
	}
	if (control.i_sens || (outSwitchS && control.i_senSwitch)) {
		Simulation::SensorOut();
		outSwitchS = 0;
	}
	if (control.i_bmp || (outSwitchP && control.i_bmpSwitch)) {
		setTwVisible(GUIBar, 0);
		static Bitmap bm;
		static int i = 0;
		char name[256];
		sprintf_s(name, "./out/bm%04d.bmp", i++);
		bm.SaveAsBMP(name);
		setTwVisible(GUIBar, 1);
		outSwitchP = 0;
	}
	if (!control.i_stop) {
		Simulation::Run();
		if (count++ % 50 == 0) {
			outSwitchS = 1;
			outSwitchP = 1;
		}
		control.i_dirty = 1;
	}
}
static void fps() {

}
static void onMouse(int button, int s, int x, int y) {
	if (!TwEventMouseButtonGLUT(button, s, x, y)) {
		control.clickMouse(button, s, x, y);
		if (button == GLUT_LEFT_BUTTON && s == GLUT_DOWN) {
			const int pickID = Visualization::IntersectColorPick(&control, Simulation::Number(), x, y);
			if (pickID == 0x00FFFFFF) return;
			const Parameters::DataType* px = (Parameters::DataType*)Simulation::PositionX();
			const Parameters::DataType* py = (Parameters::DataType*)Simulation::PositionY();
			const Parameters::DataType* div = (Parameters::DataType*)Simulation::Divergence();
			const Parameters::DataType* pres = (Parameters::DataType*)Simulation::Pressure();
			std::cout << " --------------------------------------------------------------------- " << std::endl;
			std::cout << " Particle ID : " << pickID << std::endl;
			std::cout << " Coordinate (x,y) : " << px[pickID] << ", " << py[pickID] << std::endl;
			std::cout << " Divergence : " << div[pickID] << "    " << " Pressure : " << pres[pickID] << std::endl;
			std::cout << " --------------------------------------------------------------------- " << std::endl;
		}
	}
}
static void onMotion(int x, int y) {
	if (!TwEventMouseMotionGLUT(x, y)) {
		control.moveMouse(x, y);
		glutPostRedisplay();
	}
}
static void onMouseWheel(int button, int dir, int x, int y) {
	control.rollMouse(button, dir, x, y);
}
static void onReshape(int width, int height) {
	glViewport(0, 0, width, width);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluPerspective(-90.0f, float(control.u_width) / float(control.u_height), 1.0f, 100.0f);
	control.reshapeWindow();
	TwWindowSize(control.u_width, control.u_height);
}
static void onKeyboard(unsigned char key, int x, int y) {
	if (!TwEventKeyboardGLUT(key, x, y)) {
		glutPostRedisplay();
		control.pressKey(key, x, y);
		callBack();
	}
}
static void onDisplay() {
	glm::mat4 modelMatrix = glm::translate(glm::mat4(1.0f), control.m_pan)
		* glm::toMat4(control.m_rotation)
		* glm::scale(glm::mat4(1.0f), control.m_scale);

	control.m_modelMat = modelMatrix;
	control.m_viewMat = control.m_camera.GetViewMatrix();
	control.m_viewModelMat = control.m_camera.GetViewMatrix() * modelMatrix;
	control.m_projectionMat = control.m_camera.GetProjectionMatrix();
	control.m_projectionMatInv = glm::inverse(control.m_projectionMat);
	control.m_mvp = control.m_projectionMat * control.m_viewModelMat;
	control.m_mvpInv = glm::inverse(control.m_mvp);

	Render();
	TwDraw();

	glutSwapBuffers();
	glutReportErrors();

	callBack();

	if (control.i_dirty) {
		glutPostRedisplay();
		control.i_dirty = 0;
	}
	if (control.i_leave) {
		glutLeaveMainLoop();
	}
}

void TW_CALL ButtonRun_callback(void*) {
	if (control.i_stop) {
		TwDefine(" GUI/RunStop label='Stop' ");
	}
	else {
		TwDefine(" GUI/RunStop label='Run' ");
	}
	control.i_stop = !control.i_stop;
	TwDraw();
}

static void Initialize(int argc, char** argv) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(control.u_width, control.u_height);
	glutCreateWindow("RTRenderer");
	glutMouseFunc(onMouse);
	glutMotionFunc(onMotion);
	glutMouseWheelFunc(onMouseWheel);
	glutReshapeFunc(onReshape);
	glutKeyboardFunc(onKeyboard);
	glutDisplayFunc(onDisplay);

	Simulation::Initialize();
	Visualization::Initialize();

	TwInit(TW_OPENGL, NULL);
	TwWindowSize(control.u_width, control.u_height);
	GUIBar = TwNewBar("GUI");
	TwDefine(" GUI size='180 300' ");
	TwEnumVal ev[] = { { VIS::DMODE_ONE, "Temperature" }, { VIS::DMODE_TWO, "Divergence" }, { VIS::DMODE_THREE, "Pressure" }, };
	TwType quantity = TwDefineEnum("quantity", ev, 3);
	TwAddVarRW(GUIBar, "Quantity", quantity, &control.m_mode, " group='Display' ");
	TwAddVarRW(GUIBar, "Min", TW_TYPE_FLOAT, &control.f_sRangeMin, " group='Range' ");
	TwAddVarRW(GUIBar, "Max", TW_TYPE_FLOAT, &control.f_sRangeMax, " group='Range' ");
	TwDefine(" GUI/Range group='Display' ");
	TwEnumVal ev_switch[] = { { 0, "Off" }, { 1, "On" }, };
	TwType onoff = TwDefineEnum("onoff", ev_switch, 2);
	TwAddVarRW(GUIBar, "Sensors", onoff, &control.i_senSwitch, " group='Output' ");
	TwAddVarRW(GUIBar, "Bmp", onoff, &control.i_bmpSwitch, " group='Output' ");
	TwAddButton(GUIBar, "RunStop", ButtonRun_callback, NULL, " label='Run' ");
}

static void Run() {
	glutMainLoop();
}

static void Finalize() {
	TwTerminate();
}


int _tmain(int argc, _TCHAR* argv[]) {
	CreateDirectoryA(std::string(".\\out").c_str(), NULL);
	Initialize(argc, (char**)argv);
	Run();
	Finalize();
	return 0;
}

