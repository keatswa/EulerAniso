/*
 * Payload.cpp
 *
 *  Created on: Mar. 8, 2021
 *      Author: wakeats
 */

#include "Payload.h"


//*****************************************************************************
//*****************        CONSERVED VARIABLES            *********************
//*****************************************************************************


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



void GasDynVar::setCVFromPrimitives(cfdFloat *pv) {

	cfdFloat rho = pv[0];
	cfdFloat   u = pv[1];
	cfdFloat   v = pv[2];
	cfdFloat   p = pv[3];

	PV[PV_RHO] = pv[0];
	PV[PV_U] = pv[1];
	PV[PV_V] = pv[2];
	PV[PV_P] = pv[3];
	PV[PV_A] = pv[4];

	U[CV_DENS] = rho;
	U[CV_XMOM] = rho*u;
	U[CV_YMOM] = rho*v;
	U[CV_NRG] = ((p/(GAMMA-1.0)) + 0.5*rho*(u*u + v*v));




}


//cfdFloat GasDynVar::get_U(ConservedVariable idx) {
//	return(U[idx]);
//}
//
//
//cfdFloat GasDynVar::get_PV(PrimitiveVariable idx) {
//	return(PV[idx]);
//}

void GasDynVar::resolvePrimitives() {

	PV[PV_RHO] = U[CV_DENS];
	PV[PV_U] = U[CV_XMOM]/U[CV_DENS];
	PV[PV_V] = U[CV_YMOM]/U[CV_DENS];
	PV[PV_P] = (GAMMA-1.0)*(U[CV_NRG]-0.5*U[CV_DENS]*(PV[PV_U]*PV[PV_U] + PV[PV_V]*PV[PV_V]));
	PV[PV_A] = sqrt((GAMMA*PV[PV_P]/U[CV_DENS]));

	if (isnan(PV[PV_A]))
	{
		exit(-1);
	}

}



//*****************************************************************************
//*****************        FLUXES                         *********************
//*****************************************************************************



PayloadFlux::PayloadFlux() {

}

PayloadFlux::PayloadFlux(const PayloadFlux &flux) {

	for (auto& f: all_CV) {
		F[f] = flux.F[f];
	}
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




PayloadFlux::~PayloadFlux() {
}


void PayloadFlux::setGradient(PayloadVar *cv, ORIENTATION orient) {

	if (orient == V) {
		for (auto& idx: all_PV) {
			d_PV[idx] = cv->get_d_PV_x(idx);
		}
	} else {
		for (auto& idx: all_PV) {
			d_PV[idx] = cv->get_d_PV_y(idx);
		}
	}
}


void PayloadFlux::zeroGradient() {
		for (auto& idx: all_PV) {
			d_PV[idx] = 0.0;
		}
}


GasDynFlux::GasDynFlux() {
}

GasDynFlux::GasDynFlux(const GasDynFlux &flux) : PayloadFlux(flux) {
// Copied flux vector over in base class

}



GasDynFlux::~GasDynFlux() {
}


// Used for initial conditions
void GasDynFlux::setFluxFromPrimitives(cfdFloat *pv, ORIENTATION orient) {


	cfdFloat rho = pv[0];
	cfdFloat   u = pv[1];
	cfdFloat   v = pv[2];
	cfdFloat   p = pv[3];

	switch (orient) {
	case V:     // x-fluxes
		F[0] = rho*u;
		F[1] = rho*u*u + p;
		F[2] = rho*u*v;
		F[3] = u*(p+((p/(GAMMA-1.0)) + 0.5*rho*(u*u + v*v)));
		break;
	case H:
		F[0] = rho*v;
		F[1] = rho*v*u;
		F[2] = rho*v*v + p;
		F[3] = v*(p+((p/(GAMMA-1.0)) + 0.5*rho*(u*u + v*v)));
		break;
	}

}





