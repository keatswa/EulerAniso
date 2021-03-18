/*
 * Payload.cpp
 *
 *  Created on: Mar. 8, 2021
 *      Author: wakeats
 */

#include "Payload.h"
#include <iostream>

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
	cfdFloat vsq = pv[5];

	PV[PV_RHO] = pv[0];
	PV[PV_U] = pv[1];
	PV[PV_V] = pv[2];
	PV[PV_P] = pv[3];
	PV[PV_A] = pv[4];
	PV[PV_VSQ] = pv[5];

	U[CV_DENS] = rho;
	U[CV_XMOM] = rho*u;
	U[CV_YMOM] = rho*v;
	U[CV_NRG] = ((p/(GAMMA-1.0)) + 0.5*rho*(vsq));




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
	PV[PV_VSQ] = PV[PV_U]*PV[PV_U] + PV[PV_V]*PV[PV_V];

	if (isnan(PV[PV_A]))
	{
		cout << "PV_A is NaN!" << endl;
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


// Set boundary face gradient from cell-centred gradient
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


void PayloadFlux::calcGradient(PayloadVar *cvNeg, PayloadVar *cvPos, cfdFloat ds) {

	cfdFloat ds_inv = 1.0/ds;

	for (auto& idx: all_PV) {
		d_PV[idx] = ds_inv*(cvPos->get_PV(idx) - cvNeg->get_PV(idx));
	}


	for (auto& idx: all_CV) {
		d_U[idx] = ds_inv*(cvPos->get_U(idx) - cvNeg->get_U(idx));
	}


}



GasDynFlux::GasDynFlux() {
}

GasDynFlux::GasDynFlux(const GasDynFlux &flux) : PayloadFlux(flux) {
// Copied flux vector over in base class

}



GasDynFlux::~GasDynFlux() {
}




void GasDynFlux::calcFluxes(PayloadVar *cvNeg, PayloadVar *cvPos, ORIENTATION orient) {


	calcAUSMPlusSplitMachNumber(cvNeg, cvPos, orient);





	calcAUSMPlusSplitFluxes(cvNeg, cvPos, orient);




}


void GasDynFlux::setBoundaryFluxes(PayloadVar *cv, ORIENTATION orient, BCType bcType) {

	// Calc split mach number on boundary face
	a_IF = cv->get_PV(PV_A);
	M_plus = 0.0;
	M_minus = 0.0;

	if ((bcType == INLET) || (bcType == OUTLET)) {
		if (orient == V) {
			M_plus = calcM_plus(cv->get_PV(PV_U)/a_IF);
			M_minus = calcM_minus(cv->get_PV(PV_U)/a_IF);
		}
		else {
			M_plus = calcM_plus(cv->get_PV(PV_V)/a_IF);
			M_minus = calcM_minus(cv->get_PV(PV_V)/a_IF);
		}
	}

	M = M_plus + M_minus;
	a = a_IF;



	// Now calculate fluxes (not based on mach splitting above)


	if (bcType == WALL) {  // pressure reflection

		F[CV_DENS] = 0.0;
		F[CV_NRG]  = 0.0;

		if (orient == V) {
			F[CV_XMOM] = cv->get_PV(PV_P);
			F[CV_YMOM] = 0.0;
		}
		else {
			F[CV_XMOM] = 0.0;
			F[CV_YMOM] = cv->get_PV(PV_P);
		}

	}
	else if (bcType == INLET) {  // subsonic or supersonic

		cfdFloat rho1 = 1.0; //cv->get_PV(PV_RHO);
		cfdFloat u1 = 151.655; //cv->get_PV(PV_U);
		cfdFloat v1 = 0.0; //cv->get_PV(PV_V);
		cfdFloat p1 = 100000.0; //cv->get_PV(PV_P);

		cfdFloat inletMachSquared = rho1*(u1*u1 + v1*v1)/(GAMMA*p1);

		if (inletMachSquared < 1.0) {
			if (orient == V) {
				F[CV_DENS] = rho1*u1;
				F[CV_XMOM] = F[CV_DENS]*u1 + p1;
				F[CV_YMOM] = F[CV_DENS]*v1;
				F[CV_NRG ] = 0.5*F[CV_DENS]*(u1*u1+v1*v1) + (GAMMA/(GAMMA-1.0))*u1*p1;
			}
			else {
				F[CV_DENS] = rho1*v1;
				F[CV_XMOM] = F[CV_DENS]*u1;
				F[CV_YMOM] = F[CV_DENS]*v1 + p1;
				F[CV_NRG ] = 0.5*F[CV_DENS]*(u1*u1+v1*v1) + (GAMMA/(GAMMA-1.0))*v1*p1;
			}
		}
		else {

			// Supersonic inlet flux.  No change to existing flux spec.

		}

	}
	else if (bcType == OUTLET) {  // supersonic

		cfdFloat rho1 = cv->get_PV(PV_RHO);
		cfdFloat u1 = cv->get_PV(PV_U);
		cfdFloat v1 = cv->get_PV(PV_V);
		cfdFloat p1 = cv->get_PV(PV_P);

		if (orient == V) {
			F[CV_DENS] = rho1*u1;
			F[CV_XMOM] = F[CV_DENS]*u1 + p1;
			F[CV_YMOM] = F[CV_DENS]*v1;
			F[CV_NRG ] = 0.5*F[CV_DENS]*(u1*u1+v1*v1) + (GAMMA/(GAMMA-1.0))*u1*p1;
		}
		else {
			F[CV_DENS] = rho1*v1;
			F[CV_XMOM] = F[CV_DENS]*u1;
			F[CV_YMOM] = F[CV_DENS]*v1 + p1;
			F[CV_NRG ] = 0.5*F[CV_DENS]*(u1*u1+v1*v1) + (GAMMA/(GAMMA-1.0))*v1*p1;
		}

	}

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

void GasDynFlux::calcAUSMPlusSplitMachNumber(PayloadVar *cvNeg, PayloadVar *cvPos, ORIENTATION orient) {


	cfdFloat p_l, p_r, rho_l, rho_r, a_l, a_r;
	p_r = cvPos->get_PV(PV_P);
	p_l = cvNeg->get_PV(PV_P);
	rho_r = cvPos->get_PV(PV_RHO);
	rho_l = cvNeg->get_PV(PV_RHO);
	a_r = sqrt(GAMMA*p_r/rho_r);
	a_l = sqrt(GAMMA*p_l/rho_l);
/*
	cfdFloat a_l_hat, a_r_hat, a_star_l_sq, a_star_r_sq;

	cfdFloat vsq[2];

	vsq[NEG] = cvNeg->get_PV(PV_VSQ);

	vsq[POS] = cvPos->get_PV(PV_VSQ);

	a_star_l_sq = (2.0/GAMMA+1.0)*a_l*a_l
	          + ((GAMMA-1.0)/(GAMMA+1.0))*vsq[NEG];

	a_star_r_sq = (2.0/GAMMA+1.0)*a_r*a_r
	          + ((GAMMA-1.0)/(GAMMA+1.0))*vsq[POS];

	if (orient == V)
	{
		a_l_hat = (a_star_l_sq)/(fmax(sqrt(a_star_l_sq),
		                              (cvNeg->get_PV(PV_U))));
		a_r_hat = (a_star_r_sq)/(fmax(sqrt(a_star_r_sq),
		                              -1.0*(cvPos->get_PV(PV_U))));
	}
	else
	{
		a_l_hat = (a_star_l_sq)/(fmax(sqrt(a_star_l_sq),
		                              (cvNeg->get_PV(PV_V))));
		a_r_hat = (a_star_r_sq)/(fmax(sqrt(a_star_r_sq),
		                              -1.0*(cvPos->get_PV(PV_V))));
	}

	a_IF = fmin(a_l_hat, a_r_hat);
*/
	a_IF = sqrt(a_l*a_r);


	if ( orient == V )  // horizontal-going fluxes
	{
		M_plus = calcM_plus(cvNeg->get_PV(PV_U)/a_IF);
		M_minus = calcM_minus(cvPos->get_PV(PV_U)/a_IF);
	}
	else
	{
		M_plus = calcM_plus(cvNeg->get_PV(PV_V)/a_IF);
		M_minus = calcM_minus(cvPos->get_PV(PV_V)/a_IF);
	}

	M = M_plus + M_minus;
	a = a_IF;

}

void GasDynFlux::calcAUSMPlusSplitFluxes(PayloadVar *cvNeg, PayloadVar *cvPos, ORIENTATION orient) {

	if (orient == V) {                         // horizontal-going fluxes

		M_plus = cvNeg->get_PV(PV_U)/a;
		M_IF = fmax(0.0, M);
		if (M_IF != 0.0) {
			F[CV_DENS] = M_IF*a*cvNeg->get_U(CV_DENS);
			F[CV_XMOM] = M_IF*a*cvNeg->get_U(CV_XMOM)
					          + cvNeg->get_PV(PV_P)*calcPposCoeff(M_plus);
			F[CV_YMOM] = M_IF*a*cvNeg->get_U(CV_YMOM);
			F[CV_NRG ] = M_IF*a*cvNeg->get_U(CV_NRG) + cvNeg->get_PV(PV_P);
		}
		else {
			F[CV_DENS] = 0.0;
			F[CV_XMOM] = cvNeg->get_PV(PV_P)*calcPposCoeff(M_plus);
			F[CV_YMOM] = 0.0;
			F[CV_NRG ] = 0.0;
		}

		M_minus = cvPos->get_PV(PV_U)/a;
		M_IF = fmin(0.0, M);

		if (M_IF != 0.0) {
			F[CV_DENS] += M_IF*a*cvPos->get_U(CV_DENS);
			F[CV_XMOM] += M_IF*a*cvPos->get_U(CV_XMOM)
					           + cvPos->get_PV(PV_P)*calcPnegCoeff(M_minus);
			F[CV_YMOM] += M_IF*a*cvPos->get_U(CV_YMOM);
			F[CV_NRG ] += M_IF*a*cvPos->get_U(CV_NRG) + cvPos->get_PV(PV_P);
		}
		else {
			F[CV_DENS] += 0.0;
			F[CV_XMOM] += cvPos->get_PV(PV_P)*calcPnegCoeff(M_minus);
			F[CV_YMOM] += 0.0;
			F[CV_NRG ] += 0.0;
		}


	}
	else { // orient == H                      // vertical-going fluxes

		M_plus = cvNeg->get_PV(PV_V)/a;
		M_IF = fmax(0.0, M);
		if (M_IF != 0.0) {
			F[CV_DENS] = M_IF*a*cvNeg->get_U(CV_DENS);
			F[CV_XMOM] = M_IF*a*cvNeg->get_U(CV_XMOM);
			F[CV_YMOM] = M_IF*a*cvNeg->get_U(CV_YMOM)
							  + cvNeg->get_PV(PV_P)*calcPposCoeff(M_plus);
			F[CV_NRG ] = M_IF*a*cvNeg->get_U(CV_NRG) + cvNeg->get_PV(PV_P);
		}
		else {
			F[CV_DENS] = 0.0;
			F[CV_XMOM] = 0.0;
			F[CV_YMOM] = cvNeg->get_PV(PV_P)*calcPposCoeff(M_plus);
			F[CV_NRG ] = 0.0;
		}

		M_minus = cvPos->get_PV(PV_V)/a;
		M_IF = fmin(0.0, M);

		if (M_IF != 0.0) {
			F[CV_DENS] += M_IF*a*cvPos->get_U(CV_DENS);
			F[CV_XMOM] += M_IF*a*cvPos->get_U(CV_XMOM);
			F[CV_YMOM] += M_IF*a*cvPos->get_U(CV_YMOM)
	        				   + cvPos->get_PV(PV_P)*calcPnegCoeff(M_minus);
			F[CV_NRG ] += M_IF*a*cvPos->get_U(CV_NRG) + cvPos->get_PV(PV_P);
		}
		else {
			F[CV_DENS] += 0.0;
			F[CV_XMOM] += 0.0;
			F[CV_YMOM] += cvPos->get_PV(PV_P)*calcPnegCoeff(M_minus);
			F[CV_NRG ] += 0.0;
		}



	}








}
