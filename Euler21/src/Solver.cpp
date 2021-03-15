/*
 * Solver.cpp
 *
 *  Created on: Mar. 12, 2021
 *      Author: wakeats
 */

#include "Solver.h"
#include <iostream>

Solver::Solver() {
	// TODO Auto-generated constructor stub

}

void Solver::resolveCellPrimitives() {

	for (auto& cp: mesh->cellMap) {
		Cell *c = cp.second;
		c->get_U()->resolvePrimitives();
	}

}

void Solver::calcPVGradientsAtFaces() {

	for (auto& fp: mesh->faceMap) {
		Face *f = fp.second;

		if (f->get_is_bc()) {

			if (f->get_bcType() == OUTLET) {

				f->get_F()->zeroGradient();


			} else {  // WALL, INLET, ?

				f->get_F()->setGradient(f->getBoundaryCell()->get_U(), f->get_orient());


			}




		} else {





		}





	}





}

void Solver::limitCellGradients() {

	for (auto& cp: mesh->cellMap) {
		Cell *c = cp.second;
//		c->do_stuff();
	}


}




void Solver::calcFaceFluxes() {
	for (auto& fp: mesh->faceMap) {
		Face *f = fp.second;

		// do stuff

	}

}

void Solver::doCellTimestep(cfdFloat dt, cfdInt rkStage) {
}

cfdFloat Solver::calcMinWavePeriod() {

	cfdFloat minPeriod  = 1.0e6;	// the minimum outcome thus far.
	cfdFloat tmpPeriod;
	cfdFloat tmpPeriodX;
	cfdFloat tmpPeriodY;


	for (auto& cp: mesh->cellMap) {
		Cell *c = cp.second;

		// Avoid knowledge of Gasdyn impl in the Solver class.  Get wavespeed from Payload
		tmpPeriodX = c->get_dx()/c->get_U()->get_WaveSpeed_x();
		tmpPeriodY = c->get_dx()/c->get_U()->get_WaveSpeed_y();

		tmpPeriod = min(tmpPeriodX, tmpPeriodY);

		if (tmpPeriod < minPeriod)
			minPeriod = tmpPeriod;

	}


	return minPeriod;


}


//for (auto& cp: mesh->cellMap) {
//	Cell *c = cp.second;
//	// do stuff
//}




Solver::~Solver() {
	// TODO Auto-generated destructor stub
}

void Solver::assignBackwardStepICs() {

	for (auto& c: mesh->cellMap) {

		Cell *tmpCell = c.second;

		if (tmpCell->get_x() < 0) {
			cfdFloat a0 = sqrt(GAMMA*bss.p0/bss.rho0);
			cfdFloat icPVarr[PayloadVar::N_PV] = {bss.rho0, bss.u0, bss.v0, bss.p0, a0};
			tmpCell->setCVfromPV(icPVarr);
		} else {
			cfdFloat a1 = sqrt(GAMMA*bss.p1/bss.rho1);
			cfdFloat icPVarr[PayloadVar::N_PV] = {bss.rho1, bss.u1, bss.v1, bss.p1, a1};
			tmpCell->setCVfromPV(icPVarr);
		}
	}
}


void Solver::solve() {

	cfdFloat tElapsed = 0.0;
	cfdInt   iteration = 0;
	cfdInt numOutputFile = 0;


	mesh->doUniformRefine(reflevelMin);


	do {

		// Calculate timestep
		dt = maxCFL*calcMinWavePeriod();
		if ((dt + tElapsed) > tEnd)
			dt = tEnd - tElapsed;


		resolveCellPrimitives();
		calcPVGradientsAtFaces();
		limitCellGradients();
		calcFaceFluxes();
		doCellTimestep(dt, 0);		// rkStage = 0: Euler scheme.  rkStage in {1,2,3}: Runge-Kutta


		tElapsed += dt;
		iteration++;


		if (iteration % 10 == 0)
			cout << "Iter: " << iteration << endl;


	}
	while ((tElapsed < tEnd) && (iteration < maxIter));



}
