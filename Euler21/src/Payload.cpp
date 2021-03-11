/*
 * Payload.cpp
 *
 *  Created on: Mar. 8, 2021
 *      Author: wakeats
 */

#include "Payload.h"

GasDynFlux::GasDynFlux() {
}

GasDynFlux::GasDynFlux(const GasDynFlux &flux) : PayloadFlux(flux) {

	for (auto& f: all_CV) {
		F[f] = flux.F[f];
	}

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








PayloadFlux::PayloadFlux() {
}

PayloadFlux::PayloadFlux(const PayloadFlux &flux) {
}

PayloadFlux::~PayloadFlux() {
}

PayloadVar::PayloadVar() {
}
