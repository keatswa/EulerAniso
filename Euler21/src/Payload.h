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
	PV_RHO = 0,
	PV_U   = 1,
	PV_V   = 2,
	PV_P   = 3,
	PV_A   = 4,
	PV_VSQ = 5,
	NUM_PRIMITIVE_VARS = 6
};
constexpr std::initializer_list<PrimitiveVariable> all_PV = {PV_RHO, PV_U, PV_V, PV_P, PV_A};


enum AnalyticVariable
{
	AV_SCHLIEREN = 0,
	AV_PSI0 = 1,
	AV_PSI1 = 2,
	NUM_ANALYTIC_VARS = 3
};




class PayloadVar {

private:
	inline static ProblemType problemType = GAS_DYNAMICS;  // default problemType



public:
	inline const static int N_CV = NUM_CONSERVED_VARS;  // Number of conservative vars
	inline const static int N_PV = NUM_PRIMITIVE_VARS;  // Number of primitive vars

	PayloadVar();
	PayloadVar(const PayloadVar& var);
	virtual ~PayloadVar();

	// "Virtual" constructor
	static PayloadVar* Create(ProblemType pt);

	// "Virtual" copy constructor
	virtual PayloadVar* Clone() = 0;


	inline static void setProblemType(ProblemType pt) { problemType = pt; }
	inline static ProblemType getProblemType() { return problemType; }

	virtual cfdFloat get_WaveSpeed_x() = 0;
	virtual cfdFloat get_WaveSpeed_y() = 0;


	virtual void setCVFromPrimitives(cfdFloat *pv) = 0;

	virtual cfdFloat get_U(ConservedVariable idx) { return U[idx]; }
	virtual cfdFloat get_PV(PrimitiveVariable idx) { return PV[idx]; }
	virtual cfdFloat get_d_PV_x(PrimitiveVariable idx) {return d_PV_x[idx]; }
	virtual cfdFloat get_d_PV_y(PrimitiveVariable idx) {return d_PV_y[idx]; }

	virtual void update_U(cfdFloat *newU) { for (auto& idx: all_CV) { U[idx] = newU[idx]; } }
	virtual void update_PV(cfdFloat *newPV) { for (auto& idx: all_PV) { PV[idx] = newPV[idx]; } }

	virtual void update_dU(AXIS axis, ConservedVariable idx, cfdFloat value) {
		switch (axis) {
		case X:
			d_U_x[idx] = value;
			break;
		case Y:
			d_U_y[idx] = value;
			break;
		}
	}


	virtual void update_dPV(AXIS axis, PrimitiveVariable idx, cfdFloat value) {
		switch (axis) {
		case X:
			d_PV_x[idx] = value;
			break;
		case Y:
			d_PV_y[idx] = value;
			break;
		}
	}


	virtual cfdFloat get_dU(AXIS axis, ConservedVariable idx) {
		switch (axis) {
			case X:
				return d_U_x[idx];
				break;
			case Y:
				return d_U_y[idx];
				break;
			}
		return 0;
	}

	virtual cfdFloat get_dPV(AXIS axis, PrimitiveVariable idx) {
		switch (axis) {
			case X:
				return d_PV_x[idx];
				break;
			case Y:
				return d_PV_y[idx];
				break;
			}
		return 0;
	}


	virtual void resolvePrimitives() = 0;

	virtual void set_PHI(AnalyticVariable idx, cfdFloat value) { PHI[idx] = value; }
	virtual cfdFloat get_PHI(AnalyticVariable idx) { return PHI[idx]; }

	virtual void set_U_RK(int timestep, ConservedVariable idx, cfdFloat value) {
		U_RK[timestep][idx] = value;
	}

	virtual cfdFloat get_U_RK(int timestep, ConservedVariable idx) {
		return U_RK[timestep][idx];
	}

	virtual void calc_average_from_src(PayloadVar *srcU);


protected:
	cfdFloat U[NUM_CONSERVED_VARS];
	cfdFloat PV[NUM_PRIMITIVE_VARS];
	cfdFloat d_U_x[NUM_CONSERVED_VARS];
	cfdFloat d_U_y[NUM_CONSERVED_VARS];
	cfdFloat d_PV_x[NUM_PRIMITIVE_VARS];
	cfdFloat d_PV_y[NUM_PRIMITIVE_VARS];

	cfdFloat PHI[NUM_ANALYTIC_VARS];

	// For R-K time stepping
	cfdFloat U_RK[2][NUM_CONSERVED_VARS];

};




class GasDynVar: public PayloadVar {
private:


public:
	GasDynVar();
	GasDynVar(const GasDynVar& var);

	PayloadVar* Clone() {
		return new GasDynVar(*this);
	}

	void setCVFromPrimitives(cfdFloat *pv);

//	cfdFloat get_U(ConservedVariable idx);
//	cfdFloat get_PV(PrimitiveVariable idx);

	void resolvePrimitives();
	cfdFloat get_WaveSpeed_x() { return (abs(PV[PV_U]) + PV[PV_A]); }
	cfdFloat get_WaveSpeed_y() { return (abs(PV[PV_V]) + PV[PV_A]); }

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
	cfdFloat F[NUM_CONSERVED_VARS];
	cfdFloat d_PV[NUM_PRIMITIVE_VARS];
	cfdFloat d_U[NUM_CONSERVED_VARS];


public:
	PayloadFlux();
	PayloadFlux(const PayloadFlux& flux);
	virtual ~PayloadFlux();

	// "Virtual" constructor
	static PayloadFlux* Create(ProblemType pt);

	// "Virtual" copy constructor
	virtual PayloadFlux* Clone() = 0;

	virtual void setFluxFromPrimitives(cfdFloat *pv, ORIENTATION orient) = 0;

	virtual void calcFluxes(PayloadVar *cvNeg, PayloadVar *cvPos, ORIENTATION orient, cfdFloat dsNeg, cfdFloat dsPos) = 0;

	virtual void setBoundaryFluxes(PayloadVar *cv, ORIENTATION orient, BCType bcType) = 0;

	// Set boundary face gradient from cell-centred gradient
	void setGradient(PayloadVar *cv, ORIENTATION orient);

	void zeroGradient();

	void set_d_PV(PrimitiveVariable idx, cfdFloat value) {d_PV[idx] = value; }

	cfdFloat get_d_PV(PrimitiveVariable idx) {return d_PV[idx]; }
	cfdFloat get_d_U(ConservedVariable idx) {return d_U[idx]; }

	void calcGradient(PayloadVar *cvNeg, PayloadVar *cvPos, cfdFloat ds);

	cfdFloat *get_F() { return F; }

};






//*****************************************************************************
//*****************************************************************************
class GasDynFlux: public PayloadFlux {
private:

	cfdFloat a;          // face-averaged speed of sound
	cfdFloat M;          // fusion of M_plus and M_minus
	cfdFloat M_plus;
	cfdFloat M_minus;
	cfdFloat M_IF;
	cfdFloat a_IF;


public:
	GasDynFlux();
	GasDynFlux(const GasDynFlux& flux);

	PayloadFlux* Clone() {
		return new GasDynFlux(*this);
	}

	~GasDynFlux();

	void calcFluxes(PayloadVar *cvNeg, PayloadVar *cvPos, ORIENTATION orient, cfdFloat dsNeg, cfdFloat dsPos) ;
	void setBoundaryFluxes(PayloadVar *cv, ORIENTATION orient, BCType bcType) ;

	// Used for initial conditions
	void setFluxFromPrimitives(cfdFloat *pv, ORIENTATION orient);

	void calcAUSMPlusSplitMachNumber(PayloadVar *cvNeg, PayloadVar *cvPos, ORIENTATION orient);

	void calcAUSMPlusSplitFluxes(PayloadVar *cvNeg, PayloadVar *cvPos, ORIENTATION orient);



	inline cfdFloat calcM_plus(cfdFloat M)
	{
		return 	(fabs(M) >= 1.0 ?
				0.5*(M+fabs(M)) :
		        0.25*(M+1.0)*(M+1.0)+0.125*(M*M-1.0)*(M*M-1.0));
	}

	inline cfdFloat calcM_minus(cfdFloat M)
	{
		return 	(fabs(M) >= 1.0 ?
				0.5*(M-fabs(M)) :
				-0.25*(M-1.0)*(M-1.0)-0.125*(M*M-1.0)*(M*M-1.0));
	}


	inline cfdFloat calcPposCoeff(cfdFloat M)
	{
		return fabs(M) >= 1.0 ?
		       0.5*(1.0+(M/fabs(M))) :
		       0.25*(M+1.0)*(M+1.0)*(2.0-M) + 0.1875*M*(M*M-1.0)*(M*M-1.0);

	}

	inline cfdFloat calcPnegCoeff(cfdFloat M)
	{
		return fabs(M) >= 1.0 ?
		       0.5*(1.0-(M/fabs(M))) :
		       0.25*(M-1.0)*(M-1.0)*(2.0+M) - 0.1875*M*(M*M-1.0)*(M*M-1.0);
	}






};







class ShallowWaterFlux: public PayloadFlux {



};








#endif /* PAYLOAD_H_ */
