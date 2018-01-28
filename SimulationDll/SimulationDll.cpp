/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Published under CC BY-NC
*/
// SimulationDll_.cpp : Defines the exported functions for the DLL application.
/// specify the schemes by defining different solver classes

#include "stdafx.h"
#include "SimulationDll.h"
#include "FractionalStep_KM_A_FSF.h"
#include <PreInformation.h>

namespace SIM {

	typedef FractionalStep_KM_A_FSF<Parameters::DataType, Parameters::Dimension, Parameters::Order> FS;
	typedef FS* FSPtr;

	static FSPtr objPtr;

	void SimulationDll2D::Initialize() {
		objPtr = new FS();
		objPtr->init();
	}

	void SimulationDll2D::Run() {
		objPtr->stepGL();
	}

	int SimulationDll2D::Number() {
		return objPtr->part->np;
	}
	NPtr SimulationDll2D::Type() {
		return NPtr(objPtr->type());
	}
	NPtr SimulationDll2D::PositionX() {
		return NPtr(objPtr->PositionX());
	}
	NPtr SimulationDll2D::PositionY() {
		return NPtr(objPtr->PositionY());
	}
	NPtr SimulationDll2D::VelocityX() {
		return NPtr(objPtr->VelocityX());
	}
	NPtr SimulationDll2D::VelocityY() {
		return NPtr(objPtr->VelocityY());
	}
	NPtr SimulationDll2D::Pressure() {
		return NPtr(objPtr->Pressure());
	}
	NPtr SimulationDll2D::Temperature() {
		return NPtr(objPtr->Temperature());
	}
	NPtr SimulationDll2D::Divergence() {
		return NPtr(objPtr->Divergence());
	}
	NPtr SimulationDll2D::Vorticity() {
		return NPtr(objPtr->Vorticity());
	}
	void SimulationDll2D::SaveData() {
		objPtr->saveData();
	}
	void SimulationDll2D::SensorOut() {
		objPtr->sensorOut();
	}
	void SimulationDll2D::BBox(double& left, double& right, double& bottom, double& top) {
		objPtr->part->getBBox(left, right, bottom, top);
	}
	void SimulationDll2D::Interpolation(const double x, const double y, double& u, double& v) {
		const auto* const part = objPtr->part;
		Eigen::Matrix<double, 2, 1> uv = part->interpolateLSA(part->vel[0].data(), part->vel[1].data(), x, y);
		u = uv[0];
		v = uv[1];
	}

}


