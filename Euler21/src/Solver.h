/*
 * Solver.h
 *
 *  Created on: Mar. 12, 2021
 *      Author: wakeats
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include "Mesh.h"
#include <iostream>

enum ProblemCase {
	BACKWARD_STEP,
	FORWARD_STEP,
	SHOCK_TUBE,
	SHALLOW_TANK
};



typedef struct BackwardStepStates {
	cfdFloat rho0;
	cfdFloat u0;
	cfdFloat v0;
	cfdFloat p0;
	cfdFloat rho1;
	cfdFloat u1;
	cfdFloat v1;
	cfdFloat p1;
} BackwardStepStates;




class Solver {


private:

	Mesh *mesh;

	cfdInt reflevelMax = -1;
	cfdInt reflevelMin = -1;
	cfdFloat errorTol = -1;
	cfdFloat dt = -1;
	cfdFloat tEnd = -1;
	cfdFloat maxCFL = -1;
	cfdInt maxIter = -1;

	BackwardStepStates bss;

	EulerDisplay *display;  // Could be a pointer to EulerDisplay's parent class if we have multiple possible 'display' types

	CB_EulerDisplay_drawFn drawFn;


	void calcGradientsAtCells();

	// solve() calls:
	void resolveCellPrimitives();   // Conservative Variables to Primitive Variables
	void calcPVGradientsAtFaces();       // e.g. finite difference across the face
	void limitCellGradients();      // e.g. minmod of face gradients
	void calcFaceFluxes();			// TBD: parameterize by scheme { AUSM+, VANLEER, HLLC, etc.. }
	void doCellTimestep(cfdFloat dt, cfdInt rkStage);

	cfdFloat calcMinWavePeriod();   // Loop over all cells to limit the timestep.

public:
	Solver();
	Solver(Mesh *m) { mesh = m; }

	// Specify the object that will be used to display results
	void set_display(EulerDisplay *_display) { display = _display; }

	// Set callback to draw function in 'display' object
	void set_cb_drawFn(CB_EulerDisplay_drawFn& df) { drawFn = df; }

	void setSolverParams(int _reflevelMax, int _reflevelMin, cfdFloat _errorTol, cfdFloat _dt, cfdFloat _tEnd, cfdFloat _maxCFL, int _maxIter) {
		reflevelMax = _reflevelMax;
		reflevelMin = _reflevelMin;
		errorTol = _errorTol;
		dt = _dt;
		tEnd = _tEnd;
		maxCFL = _maxCFL;
		maxIter = _maxIter;
	}

	void setBackwardStepParams(cfdFloat _r0, cfdFloat _u0, cfdFloat _v0, cfdFloat _p0, cfdFloat _r1, cfdFloat _u1, cfdFloat _v1, cfdFloat _p1) {
		bss.rho0 = _r0;
		bss.u0   = _u0;
		bss.v0   = _v0;
		bss.p0   = _p0;
		bss.rho1 = _r1;
		bss.u1   = _u1;
		bss.v1   = _v1;
		bss.p1   = _p1;
	}


//	BackwardStepStates getBackwardStepParams() { return bss; }


	// Apply the two-state initial condition to the grid of cells that makes up the Backward Step.
	// If desired, pre-refine the mesh before applying initial conditions based on geometry.
	// TBD!
	void assignBackwardStepICs();


	void setInletCondition(unsigned long fid, cfdFloat rho, cfdFloat u, cfdFloat v, cfdFloat p) {

		cout << "Inlet fid: " << fid << ", u = " << u << ", v = " << v << endl;
		cfdFloat inletPVarr[PayloadVar::N_PV] = {rho, u, v, p};
		mesh->faceMap[fid]->setBCFlux(inletPVarr);

	}


	void solve();



	virtual ~Solver();
};

#endif /* SOLVER_H_ */
