/*
 * Solver.cpp
 *
 *  Created on: Mar. 12, 2021
 *      Author: wakeats
 */

#include "Solver.h"
#include <iostream>
#include <omp.h>

Solver::Solver() {
	// TODO Auto-generated constructor stub

}

void Solver::resolveCellPrimitives() {

	mesh->init_vecCellIDs();
	unsigned long Ncells = mesh->vecCellIDs.size();
	unsigned long i;
	Cell *c;
	PayloadVar *tmpU;

	BROKEN
#pragma omp parallel for shared(mesh) private(c, tmpU) num_threads(8)
	for (i = 0 ; i < Ncells ; i++) {

//	for (auto& cp: mesh->cellMap) {
//		Cell *c = cp.second;

		c = mesh->cellMap.at(mesh->vecCellIDs[i]);
		tmpU = c->get_U();
		tmpU->resolvePrimitives();
//		c->get_U()->resolvePrimitives();
	}

}



void Solver::calcGradientsAtCells() {

	Cell *tmpC = NULL;
	cfdFloat dx, dy;
	cfdFloat rx, ry;


	cfdFloat sUWXX; //[NUM_CONSERVED_VARS];
	cfdFloat sUWXY; //[NUM_CONSERVED_VARS];
	cfdFloat sUWYY; //[NUM_CONSERVED_VARS];

	cfdFloat sUWX[NUM_CONSERVED_VARS];
	cfdFloat sUWY[NUM_CONSERVED_VARS];
	cfdFloat wtU ; //[NUM_CONSERVED_VARS]; // weight

	cfdFloat dU[NUM_CONSERVED_VARS];

	int weightX, weightY;






	cfdFloat sPWXX[NUM_PRIMITIVE_VARS];
	cfdFloat sPWXY[NUM_PRIMITIVE_VARS];
	cfdFloat sPWYY[NUM_PRIMITIVE_VARS];

	cfdFloat sPWX[NUM_PRIMITIVE_VARS];


	for (auto& cp: mesh->cellMap) {
		Cell *c = cp.second;

		dx = c->get_dx();
		dy = c->get_dy();

		wtU = 0.0;
		sUWXX = 0.0;
		sUWXY = 0.0;
		sUWYY = 0.0;

		int i;

		for (i = 0 ; i < NUM_CONSERVED_VARS ; i++) {
			sUWX[i] = 0.0;
			sUWY[i] = 0.0;
		}





		for (auto& d: all_DIR) {

			for_each(c->get_nbFaces().at(d)->begin(), c->get_nbFaces().at(d)->end(),
						[&c, &d, &rx, &ry, &wtU, &sUWXX, &sUWXY, &sUWYY, &sUWX, &sUWY, &dU, &weightX, &weightY ] (Face *f)
			{
				if (!f->get_is_bc()) {

					Cell *nb = NULL;
					if (!f->get_nbCell(nb, d)) {exit(-1);}
					rx = nb->get_x() - c->get_x();
					ry = nb->get_y() - c->get_y();
					cfdFloat tmpWeight = 1.0/(rx*rx + ry*ry);

					wtU = tmpWeight;
					sUWXX += wtU*rx*rx;
					sUWXY += wtU*rx*ry;
					sUWYY += wtU*ry*ry;

					weightX = dirIsXAxis(d);
					weightY = dirIsYAxis(d);

					for (auto& idx: all_CV) {
						dU[idx] = nb->get_U(idx) - c->get_U(idx);
						sUWX[idx] += weightX*wtU*rx*dU[idx];
						sUWY[idx] += weightY*wtU*ry*dU[idx];

					} // all_CV
				}  // if face is not a BC
			});  // faces at DIR d

		} // for DIR


		cfdFloat denom = sUWXY*sUWXY - sUWXX*sUWYY;

		// Denominator nonzero:
		if (fabs(denom) > 1e-6) {
			for (auto& idx: all_CV) {
				c->get_U()->update_dU(X, idx, (sUWXY*sUWY[idx] - sUWYY*sUWX[idx])/denom);
				c->get_U()->update_dU(Y, idx, (sUWXY*sUWX[idx] - sUWXX*sUWY[idx])/denom);
			} // all_CV

		}
		else {
			for (auto& idx: all_CV) {
				c->get_U()->update_dU(X, idx, 0.0);
				c->get_U()->update_dU(Y, idx, 0.0);
			}
		}




	}  // for each Cell





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
	Face *f = NULL;
	PayloadVar *cvNeg, *cvPos;
	unsigned long i;
	unsigned long Nfaces = mesh->vecFaceIDs.size();


	mesh->init_vecFaceIDs();


//	for (auto& fp: mesh->faceMap) {
//		f = fp.second;

#pragma omp parallel for shared(mesh) private(c, f, cvNeg, cvPos) num_threads(8)

	for (i = 0 ; i < Nfaces ; i++) {

		f = mesh->faceMap.at(mesh->vecFaceIDs[i]);

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


			f->get_F()->calcFluxes(cvNeg, cvPos, orient);


		}

	}

}


void Solver::doCellTimestep(cfdFloat dt, cfdInt rkStage) {

	cfdFloat dx, dy;
	int fcnt;
	cfdFloat flux[4][PayloadVar::N_CV];
	cfdFloat newU[PayloadVar::N_CV];

	mesh->init_vecCellIDs();

	unsigned long Ncells = mesh->vecCellIDs.size();

//	for (auto& cp: mesh->cellMap) {
//		Cell *c = cp.second;

	unsigned long i;
	Cell *c;
	DIR d;
	ConservedVariable idx;
//	CellMap *cm = &(mesh->cellMap);


#pragma omp parallel for shared(mesh/*, cm*/) private(c, flux, newU, idx, d, fcnt) num_threads(8)
	for (i = 0 ; i < Ncells ; i++) {

//		c = cm->at(mesh->vecCellIDs[i]);

		c = mesh->cellMap.at(mesh->vecCellIDs[i]);


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
			cfdFloat vsq0 = bss.u0*bss.u0 +bss.v0*bss.v0;
			cfdFloat icPVarr[PayloadVar::N_PV] = {bss.rho0, bss.u0, bss.v0, bss.p0, a0, vsq0};
			tmpCell->setCVfromPV(icPVarr);
		} else {
			cfdFloat a1 = sqrt(GAMMA*bss.p1/bss.rho1);
			cfdFloat vsq1 = bss.u1*bss.u1 +bss.v1*bss.v1;
			cfdFloat icPVarr[PayloadVar::N_PV] = {bss.rho1, bss.u1, bss.v1, bss.p1, a1, vsq1};
			tmpCell->setCVfromPV(icPVarr);
		}
	}
}


void Solver::solve() {

	cfdFloat tElapsed = 0.0;
	cfdInt   iteration = 0;
	cfdInt numOutputFile = 0;


	mesh->doUniformRefine(reflevelMin);
//	mesh->doUniformRefine(2);

	resolveCellPrimitives();

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

		resolveCellPrimitives();

		calcGradientsAtCells();

		tElapsed += dt;
		iteration++;

		if (((display)->*(drawFn))(mesh->cellMap, mesh->faceMap) == -1)
			exit(1);

		if (iteration % 10 == 0)
		{
			cout << "Iter: " << iteration << endl;

			ulong testID0 = Cell::generate_id(0, 16, reflevelMin, reflevelMin);

			for (auto& c: mesh->cellMap) {
				Cell *tmpCell = c.second;
				if (tmpCell->get_i_idx() == 0 && tmpCell->get_j_idx() == 16) {
					cout << "testID: " << testID0  << ", j_idx: " << tmpCell->get_j_idx() << endl;
					cout << " genID: " << Cell::generate_id(tmpCell->get_i_idx(), tmpCell->get_j_idx(),
											  tmpCell->get_li(),    tmpCell->get_lj()) << endl;
				}

				if (tmpCell->get_i_idx() == 10 && tmpCell->get_j_idx() == 17) {
					cout << "testID: " << testID0  << ", j_idx: " << tmpCell->get_j_idx() << endl;
					cout << " genID: " << Cell::generate_id(tmpCell->get_i_idx(), tmpCell->get_j_idx(),
											  tmpCell->get_li(),    tmpCell->get_lj()) << endl;
				}


			}
		}

	}
	while ((tElapsed < tEnd) && (iteration < maxIter));



}
