/*
 * Payload.cpp
 *
 *  Created on: Mar. 8, 2021
 *      Author: wakeats
 */

#include "Payload.h"


PayloadFlux::PayloadFlux() {

}

PayloadFlux::PayloadFlux(const PayloadFlux &flux) {

	for (auto& f: all_CV) {
		F[f] = flux.F[f];
	}
}

PayloadFlux::~PayloadFlux() {
}






GasDynFlux::GasDynFlux() {
}

GasDynFlux::GasDynFlux(const GasDynFlux &flux) : PayloadFlux(flux) {
// Copied flux vector over in base class

}



GasDynFlux::~GasDynFlux() {
}


// STATIC
PayloadFlux *PayloadFlux::Create(ProblemType pt) {

		switch (pt) {
		case GAS_DYNAMICS:
			return new GasDynFlux();
			break;
		case SHALLOW_WATER:
			;
	//		return new ShallowWaterFlux();
		}

		return NULL;
}








PayloadVar::PayloadVar() {
}

PayloadVar::PayloadVar(const PayloadVar &var) {

//	for (auto& f: all_CV) {
//		U[f] = var.U[f];
//	}

	for (int i = 0 ; i < N_CV ; i++)
		U[i] = var.U[i];

	for (int i = 0 ; i < N_PV ; i++)
		PV[i] = var.PV[i];


}

PayloadVar* PayloadVar::Create(ProblemType pt) {

	switch (pt) {
	case GAS_DYNAMICS:
		return new GasDynVar();
		break;
	case SHALLOW_WATER:
		;
//		return new ShallowWaterVar();
	}

	return NULL;


}

PayloadVar::~PayloadVar() {
}


GasDynVar::GasDynVar() {

}

GasDynVar::GasDynVar(const GasDynVar &var): PayloadVar(var) {
// Copied CV vector over in base class

}
