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
#include "FractionalStep_KM_A.h"
#include "FractionalStep_DD.h"
#include "FractionalStep_X.h"
#include <PreInformation.h>

namespace SIM {

	typedef FractionalStep_KM_A<Parameters::DataType, Parameters::Dimension, Parameters::Order> FS;
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
	void SimulationDll2D::SaveData() {
		objPtr->saveData();
	}
	void SimulationDll2D::SensorOut() {
		objPtr->sensorOut();
	}

}


