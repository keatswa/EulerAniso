/*
 * Payload.h
 *
 *  Created on: Mar. 8, 2021
 *      Author: wakeats
 */

#ifndef PAYLOAD_H_
#define PAYLOAD_H_

#include "Euler.h"
#include <functional>

enum ProblemType
{
	GAS_DYNAMICS,
	SHALLOW_WATER
};





enum ConservedVariable
{
	CV_DENS = 0,
	CV_XMOM = 1,
	CV_YMOM = 2,
	CV_NRG  = 3,
	NUM_CONSERVED_VARS = 4
};

constexpr std::initializer_list<ConservedVariable> all_CV = {CV_DENS, CV_XMOM, CV_YMOM, CV_NRG};



enum PrimitiveVariable
{
	PV_XVEL,
	PV_YVEL,
	PV_PRES,
	PV_A
};




// Abstract
class PayloadVar {

private:
	inline static ProblemType problemType = GAS_DYNAMICS;  // default problemType



public:
	inline const static int N_CV = 4;  // Number of conservative vars
	inline const static int N_PV = 4;  // Number of primitive vars

	PayloadVar();
	PayloadVar(const PayloadVar& var);
	virtual ~PayloadVar();

	// "Virtual" constructor
	static PayloadVar* Create(ProblemType pt);

	// "Virtual" copy constructor
	virtual PayloadVar* Clone() = 0;


	inline static void setProblemType(ProblemType pt) { problemType = pt; }
	inline static ProblemType getProblemType() { return problemType; }




protected:
	cfdFloat U[N_CV];
	cfdFloat PV[N_PV];



};




class GasDynVar: public PayloadVar {
private:


public:
	GasDynVar();
	GasDynVar(const GasDynVar& var);

	PayloadVar* Clone() {
		return new GasDynVar(*this);
	}


};



class ShallowWaterVar: public PayloadVar {



};




//*****************************************************************************
//*****************************************************************************
//******************             FLUXES               *************************
//*****************************************************************************
//*****************************************************************************




// Abstract Base
class PayloadFlux {


protected:   // grant access to derived classes but not others
	cfdFloat F[ConservedVariable::NUM_CONSERVED_VARS];



public:
	PayloadFlux();
	PayloadFlux(const PayloadFlux& flux);
	virtual ~PayloadFlux();

	// "Virtual" constructor
	static PayloadFlux* Create(ProblemType pt);

	// "Virtual" copy constructor
	virtual PayloadFlux* Clone() = 0;

	virtual int calcFluxes(PayloadVar *cvNeg, PayloadVar *cvPos) = 0;


	virtual void setBCFluxFromPrimitives(cfdFloat *pv, ORIENTATION orient) = 0;




};






//*****************************************************************************
//*****************************************************************************
class GasDynFlux: public PayloadFlux {
private:



public:
	GasDynFlux();
	GasDynFlux(const GasDynFlux& flux);

	PayloadFlux* Clone() {
		return new GasDynFlux(*this);
	}

	~GasDynFlux();

	int calcFluxes(PayloadVar *cvNeg, PayloadVar *cvPos) {

		F[0] = 0.0;

		return 0;
	}


	void setBCFluxFromPrimitives(cfdFloat *pv, ORIENTATION orient) {

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
};







class ShallowWaterFlux: public PayloadFlux {



};








#endif /* PAYLOAD_H_ */
