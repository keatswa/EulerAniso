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

	Cell *c = NULL;
	PayloadVar *cvNeg, *cvPos;
	cfdFloat ds, ds_neg, ds_pos;
	ORIENTATION orient;


	for (auto& fp: mesh->faceMap) {
		Face *f = fp.second;

		orient = f->get_orient();

		if (f->get_is_bc()) {
			if (f->get_bcType() == OUTLET) {
				f->get_F()->zeroGradient();
			} else {  // WALL, INLET, ?
				f->get_F()->setGradient(f->getBoundaryCell()->get_U(), orient);
			}

		} else {   // FACE IS NOT A BOUNDARY - calc gradient from nb cells

			if (f->get_nbCell(c, NEG)) {
				cvNeg = c->get_U();
				ds_neg = c->get_ds_normal(orient);
			}
			else
				exit(-1);

			if (f->get_nbCell(c, POS)) {
				cvPos = c->get_U();
				ds_pos = c->get_ds_normal(orient);
			}
			else
				exit(-1);

			ds = ds_pos + ds_neg;

			f->get_F()->calcGradient(cvNeg, cvPos, ds);
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

	ORIENTATION orient;
	Cell *c = NULL;
	PayloadVar *cvNeg, *cvPos;


	for (auto& fp: mesh->faceMap) {
		Face *f = fp.second;


		orient = f->get_orient();

		if (f->get_is_bc()) {

			f->get_F()->setBoundaryFluxes(f->getBoundaryCell()->get_U(), orient, f->get_bcType());



		} else {   // FACE IS NOT A BOUNDARY

			if (f->get_nbCell(c, NEG)) {
				cvNeg = c->get_U();
			}
			else
				exit(-1);

			if (f->get_nbCell(c, POS)) {
				cvPos = c->get_U();
			}
			else
				exit(-1);


			f->get_F()->calcFluxes(cvNeg, cvPos);


		}

	}

}

void Solver::doCellTimestep(cfdFloat dt, cfdInt rkStage) {

	cfdFloat dx, dy;
	int fcnt;
	cfdFloat flux[4][PayloadVar::N_CV];
	cfdFloat newU[PayloadVar::N_CV];

	for (auto& cp: mesh->cellMap) {
		Cell *c = cp.second;

		dx = c->get_dx();
		dy = c->get_dy();

		for (auto& d: all_DIR) {

			for (auto& idx: all_CV) { flux[d][idx] = 0.0; }

			fcnt = 1;
			for_each(c->get_nbFaces().at(d)->begin(), c->get_nbFaces().at(d)->end(),
						[&c, &flux, &d, &fcnt] (Face *f)
			{
				cfdFloat *tmpFlux = f->get_F()->get_F();
				for (auto& idx: all_CV) {
					flux[d][idx] += tmpFlux[idx];
					flux[d][idx] /= fcnt;
				}
				fcnt++;
			});

		} // for DIR


		for (auto& idx: all_CV) {
			newU[idx] = c->get_U(idx) - (dt/dx)*(flux[E][idx] - flux[W][idx])
					                  - (dt/dy)*(flux[N][idx] - flux[S][idx]);
		}

		c->get_U()->update_U(newU);



	}  // for Cells







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
		tmpPeriodY = c->get_dy()/c->get_U()->get_WaveSpeed_y();

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

		Cell::set_dt(dt);

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
		{
			cout << "Iter: " << iteration << endl;
//			((display)->*(drawFn))(mesh->cellMap, mesh->faceMap);
		}

	}
	while ((tElapsed < tEnd) && (iteration < maxIter));



}
