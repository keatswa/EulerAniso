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






// Abstract
class PayloadVar {

private:

	inline static ProblemType problemType = GAS_DYNAMICS;  // default problemType

public:
	PayloadVar();
	virtual ~PayloadVar() = 0;

	inline static void setProblemType(ProblemType pt) {
		problemType = pt;
	}

	inline static ProblemType getProblemType() {
		return problemType;
	}

};




class GasDynVar: public PayloadVar {



};



class ShallowWaterVar: public PayloadVar {



};


//typedef GasDynVar  ConsVar;


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

};







class ShallowWaterFlux: public PayloadFlux {



};








#endif /* PAYLOAD_H_ */
