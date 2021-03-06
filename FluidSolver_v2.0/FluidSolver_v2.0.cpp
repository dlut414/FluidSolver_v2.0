/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Released under CC BY-NC
*/
// Manipulator.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
///Manipulator.cpp main loop
#pragma once
#include <Controller.h>
#include "Header.h"
#include "Bitmap.h"
#include "../VisualizationDll/VisualizationDll.h"
#include "../SimulationDll/SimulationDll.h"
#include "PreInformation.h"

typedef VIS::VisualizationDll Visualization;
typedef SIM::SimulationDll2D Simulation;

static VIS::Controller control;
static TwBar* GUIBar;
static TwBar* StreamBar;
static void setTwVisible(TwBar* const bar, const int visible);

static void Render() {
	//Visualization::Run(&control, Parameters::Dimension, Simulation::Number(), Simulation::Type(), Simulation::Position(), Simulation::Scalar());
	switch (control.m_mode) {
	case VIS::DMODE_ONE:
		Visualization::Run(&control, Simulation::Number(), Simulation::Type(), Simulation::PositionX(), Simulation::PositionY(), Simulation::Vorticity());
		setTwVisible(StreamBar, 0);
		break;
	case VIS::DMODE_TWO:
		Visualization::Run(&control, Simulation::Number(), Simulation::Type(), Simulation::PositionX(), Simulation::PositionY(), Simulation::Divergence());
		setTwVisible(StreamBar, 0);
		break;
	case VIS::DMODE_THREE:
		Visualization::Run(&control, Simulation::Number(), Simulation::Type(), Simulation::PositionX(), Simulation::PositionY(), Simulation::Pressure());
		setTwVisible(StreamBar, 0);
		break;
	case VIS::DMODE_FOUR:
		Visualization::Run(&control, Simulation::Number(), Simulation::Type(), Simulation::PositionX(), Simulation::PositionY(), Simulation::VelocityX());
		setTwVisible(StreamBar, 0);
		break;
	case VIS::DMODE_FIVE:
		Visualization::Run(&control, Simulation::Number(), Simulation::Type(), Simulation::PositionX(), Simulation::PositionY(), Simulation::VelocityY());
		setTwVisible(StreamBar, 0);
		break;
	case VIS::DMODE_SIX:
		setTwVisible(StreamBar, 1);
		Visualization::Run_stream(&control, 0, Simulation::Interpolation);
		break;
	default:
		setTwVisible(StreamBar, 0);
		break;
	}
}

static void setTwVisible(TwBar* const bar, const int visible) {
	TwSetParam(bar, NULL, "visible", TW_PARAM_INT32, 1, &visible);
}

static void callBack() {
	static int outSwitchS = 0;
	static int outSwitchP = 0;
	static int count = 0;
	if (control.i_save) {
		Simulation::SaveData();
		control.i_save = 0;
	}
	if (count % 500 == 0 && count != 0) Simulation::SaveData();
	if (control.i_sens || (outSwitchS && control.i_senSwitch)) {
		Simulation::SensorOut();
		outSwitchS = 0;
		control.i_sens = 0;
	}
	if (control.i_bmp || (outSwitchP && control.i_bmpSwitch)) {
		setTwVisible(GUIBar, 0);
		setTwVisible(StreamBar, 0);
		Render();
		static Bitmap bm;
		static int i = 0;
		char name[256];
		sprintf_s(name, "./out/snap%04d.png", i++);
		//bm.SaveAsBMP(name);
		bm.SaveAsPNG(name);
		setTwVisible(GUIBar, 1);
		setTwVisible(StreamBar, 1);
		outSwitchP = 0;
		control.i_bmp = 0;
	}
	if (!control.i_stop) {
		Simulation::Run();
		if (count++ % 50 == 0) {
			outSwitchS = 1;
			outSwitchP = 1;
		}
		//control.i_dirty = 1;
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
			const Parameters::DataType* const px = (Parameters::DataType*)Simulation::PositionX();
			const Parameters::DataType* const py = (Parameters::DataType*)Simulation::PositionY();
			const Parameters::DataType* const div = (Parameters::DataType*)Simulation::Divergence();
			const Parameters::DataType* const pres = (Parameters::DataType*)Simulation::Pressure();
			const Parameters::DataType* const ux = (Parameters::DataType*)Simulation::VelocityX();
			const Parameters::DataType* const uy = (Parameters::DataType*)Simulation::VelocityY();
			std::cout << " --------------------------------------------------------------------- " << std::endl;
			std::cout << " Particle ID : " << pickID << std::endl;
			std::cout << " Coordinate (x,y) : " << px[pickID] << ", " << py[pickID] << std::endl;
			std::cout << " Velocity (x,y) : " << ux[pickID] << ", " << uy[pickID] << std::endl;
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
	glViewport(0, 0, width, height);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	double left, right, bottom, top;
	Simulation::BBox(left, right, bottom, top);
	control.reshapeWindow(width, height, float(left), float(right), float(bottom), float(top));
	//gluPerspective(-90.0f, float(control.u_width) / float(control.u_height), 1.0f, 100.0f);
	TwWindowSize(control.u_width, control.u_height);
}
static void onKeyboard(unsigned char key, int x, int y) {
	if (!TwEventKeyboardGLUT(key, x, y)) {
		glutPostRedisplay();
		control.pressKey(key, x, y);
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

	glutPostRedisplay();
	if (control.i_leave) {
		glutLeaveMainLoop();
	}
}

static void onIdle() {
	return;
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

void TW_CALL ButtonRenderStream_callback(void*) {
	Visualization::Run_stream(&control, 1, Simulation::Interpolation);
	TwDraw();
	glutSwapBuffers();
	glutReportErrors();
}

static void Initialize(int argc, char** argv) {
	Simulation::Initialize();
	double left, right, bottom, top;
	Simulation::BBox(left, right, bottom, top);
	double whRatio = (right - left) / (top - bottom);
	control.u_height = GLuint(control.u_width / whRatio) < control.u_height_max ? GLuint(control.u_width / whRatio) : control.u_height_max;
	control.u_width = GLuint(control.u_width / whRatio) < control.u_height_max ? control.u_width : GLuint(control.u_height_max * whRatio);
	control.setProjectionOR(float(left), float(right), float(bottom), float(top));

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
	glutIdleFunc(onIdle);

	Visualization::Initialize();

	TwInit(TW_OPENGL, NULL);
	TwWindowSize(control.u_width, control.u_height);
	GUIBar = TwNewBar("GUI");
	TwDefine(" GUI size='180 300' position='0 0' ");
	TwEnumVal ev[] = { { VIS::DMODE_ONE, "Vorticity" }, { VIS::DMODE_TWO, "Divergence" }, { VIS::DMODE_THREE, "Pressure" }, { VIS::DMODE_FOUR, "VelocityX" }, { VIS::DMODE_FIVE, "VelocityY" }, { VIS::DMODE_SIX, "Streamline" } };
	TwType quantity = TwDefineEnum("quantity", ev, 6);
	TwAddVarRW(GUIBar, "Quantity", quantity, &control.m_mode, " group='Display' ");
	TwAddVarRW(GUIBar, "Min", TW_TYPE_FLOAT, &control.f_sRangeMin, " group='Range' ");
	TwAddVarRW(GUIBar, "Max", TW_TYPE_FLOAT, &control.f_sRangeMax, " group='Range' ");
	TwDefine(" GUI/Range group='Display' ");
	TwEnumVal ev_switch[] = { { 0, "Off" }, { 1, "On" }, };
	TwType onoff = TwDefineEnum("onoff", ev_switch, 2);
	TwAddVarRW(GUIBar, "Sensors", onoff, &control.i_senSwitch, " group='Output' ");
	TwAddVarRW(GUIBar, "Snapshot", onoff, &control.i_bmpSwitch, " group='Output' ");
	TwAddButton(GUIBar, "RunStop", ButtonRun_callback, NULL, " label='Run' ");

	StreamBar = TwNewBar("Streamline");
	TwDefine(" Streamline size='250 200' position='180 0' ");
	TwAddVarRW(StreamBar, "Px", TW_TYPE_FLOAT, &control.v_p1.x, " group='P1' ");
	TwAddVarRW(StreamBar, "Py", TW_TYPE_FLOAT, &control.v_p1.y, " group='P1' ");
	TwAddVarRW(StreamBar, "Qx", TW_TYPE_FLOAT, &control.v_p2.x, " group='P2' ");
	TwAddVarRW(StreamBar, "Qy", TW_TYPE_FLOAT, &control.v_p2.y, " group='P2' ");
	TwAddVarRW(StreamBar, "# of streamlines", TW_TYPE_INT32, &control.i_nLines, "  ");
	TwAddVarRW(StreamBar, "Integration step size", TW_TYPE_FLOAT, &control.f_dStep, "  ");
	TwAddVarRW(StreamBar, "Streamline length", TW_TYPE_FLOAT, &control.f_sLength, "  ");
	TwAddButton(StreamBar, "Render", ButtonRenderStream_callback, NULL, " label='Render' ");
	setTwVisible(StreamBar, 0);
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

