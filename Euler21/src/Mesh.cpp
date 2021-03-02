/*
 * Mesh.cpp
 *
 *  Created on: Feb. 3, 2021
 *      Author: wakeats
 */

#include <functional>
#include "Mesh.h"
#include <vector>
#include <iostream>

#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))



Mesh::Mesh() {
	// TODO Auto-generated constructor stub

}

Mesh::~Mesh() {
	// TODO Auto-generated destructor stub
}

//void testfn() { cout << "test.." << endl; }


// Refine from level 0 to target 'reflvl'
void Mesh::doUniformRefine(unsigned int reflvl) {

	int numRefinedAtLvl = 0;

	for (unsigned int lvl = 0 ; lvl < reflvl ; lvl++) {

		CellMap::iterator citor = cellMap.begin();
		while (citor != cellMap.end()) {
			Cell *tgtCell = citor->second;
			tgtCell->refFlags.set(doRefineX);
			tgtCell->refFlags.set(doRefineY);
			tgtCell->refFlags.reset(wasXRefined);
			tgtCell->refFlags.reset(wasYRefined);
			tgtCell->refFlags.reset(doCoarsenX);
			tgtCell->refFlags.reset(doCoarsenY);
			*citor++;
		}

		numRefinedAtLvl = refineCells(lvl);

		if (numRefinedAtLvl > 0)
			maxRefLvl = lvl;        // Better would be to track this independently by polling the Faces/Cells in the Mesh
	}
}



// Attempt to coarsen cells to reach 'reflvl' refinement level
void Mesh::doUniformCoarsen(unsigned int reflvl) {

	int numCoarsenedAtLvl = 0;

	// Highest current refinement level = Mesh::maxRefLvl

	FaceMap::iterator fitor = faceMap.begin();
	while (fitor != faceMap.end()) {
		Face *tgtFace = fitor->second;
		tgtFace->faceRefFlags.set(doCoarsen);
		*fitor++;
	}

	numCoarsenedAtLvl = coarsenCells(reflvl);

	cout << "Coarsened " << numCoarsenedAtLvl << " Faces." << endl;

}



// Refine cells whose x/y refinement level is reflvl.  Refine in X and then in Y.
unsigned int Mesh::refineCells(unsigned int reflvl) {

	unsigned int numRefined = 0;
	unsigned int xRefined = 0;
	unsigned int yRefined = 0;
	unsigned int i;
	vector<unsigned int> *cellIDVec = new vector<unsigned int>();

	// Define and populate vector with cell identifiers present in cellMap.
	// We will operate only on these IDs since new ones will be added to the map.



	// refine in Y
	cellIDVec->resize(cellMap.size());
	i = 0;
	for (auto& c: cellMap) {
		(*cellIDVec)[i] = c.first;
		i++;
	}

	for (auto c: (*cellIDVec)) {
		Cell *tgtCell = cellMap[(*cellIDVec)[c]];

		if (tgtCell->refFlags.test(doRefineY) && (tgtCell->get_lj() == reflvl)){
			if (refineY(tgtCell)) {
				yRefined++;
			}
		}

		((display)->*(drawFn))(cellMap, faceMap);
// Alternatively, but requires g++-10:
//		invoke(drawFn, display, cellMap, faceMap); // @suppress("Invalid arguments")

	}
	cout << "Y-refined " << yRefined << " cells." << endl;

	// refine in X
	cellIDVec->resize(cellMap.size());
	i = 0;
	for (auto& c: cellMap) {
		(*cellIDVec)[i] = c.first;
		i++;
	}

	for (auto c: (*cellIDVec)) {
		Cell *tgtCell = cellMap[(*cellIDVec)[c]];

		if (tgtCell->refFlags.test(doRefineX) && (tgtCell->get_li() == reflvl)){
			if (refineX(tgtCell)) {
				xRefined++;
			}
		}
		((display)->*(drawFn))(cellMap, faceMap);
//		std::invoke(drawFn, display, cellMap, faceMap); // @suppress("Invalid arguments")

	}

	cout << "X-refined " << xRefined << " cells." << endl;

	numRefined = xRefined + yRefined;

	// Clean up
	cellIDVec->clear();
	delete(cellIDVec);

	cout << "Refined " << numRefined << " cells." << endl;
	return numRefined;
}

// Attempt to coarsen cells currently at 'reflvl' refinement level
unsigned int Mesh::coarsenCells(unsigned int reflvl) {

	unsigned int numCoarsened = 0;
	unsigned int i;

	vector<unsigned int> *faceIDVec = new vector<unsigned int>();

	// Define and populate vector with cell identifiers present in cellMap.
	// We will operate only on these IDs since new ones will be added to the map.



	// Coarsen in X
	faceIDVec->resize(faceMap.size());
	i = 0;
	for (auto& f: faceMap) {
		(*faceIDVec)[i] = f.first;
		i++;
	}

	for (auto f: (*faceIDVec)) {
		Face *tgtFace = faceMap[(*faceIDVec)[f]];

		if (tgtFace->faceRefFlags.test(doCoarsen)
				&& ( tgtFace->get_orient() == V)
				&& (!tgtFace->faceRefFlags.test(doRecycleFace)) ) {
			if (coarsenX(tgtFace, reflvl)) {
				tgtFace->faceRefFlags.reset(doCoarsen);
				tgtFace->faceRefFlags.set(doRecycleFace);
				numCoarsened++;
			}
		}

		((display)->*(drawFn))(cellMap, faceMap);
// Alternatively, but requires g++-10:
//		invoke(drawFn, display, cellMap, faceMap); // @suppress("Invalid arguments")

	}
	cout << "Coarsened " << numCoarsened << " Faces." << endl;





	return numCoarsened;
}


bool Mesh::refineX(Cell *c0) {

	cout << "X-Refining cell id " << c0->get_id() << endl;

	// Test for ability to refine
	// 1. Check refinement flags
	if (c0->refFlags.test(doCoarsenX) || c0->refFlags.test(wasXRefined))
		return false;

	// 2. Check that neighbouring cells N,S,E,W are not preventing refinement
	Cell *testNbCell = NULL;

	for (auto& d: all_DIR) {
		if (c0->get_nbFaces().at(d)->front()->get_nbCell(testNbCell, d)) // if a neighbour to direction 'd' exists
			if (testNbCell->get_li() == c0->get_li()-1)                  // if the neighbour is twice as wide
				if (refineX(testNbCell) == false)                        // if unable to refine wider neighbour
					return false;                                        // give up attempt to refine this cell
	}


	// 3. If above checks are passed, refine cell.
	//    Halve the size of this cell and add another cell to the South or West.
	//  a) First, obtain a new cell identifier

	//  b) Create new cell to the West
	//     Static factory method will be used to initialize, set flags, and remap/create faces
	unsigned int newCellID = provideNewCellID();
	RefinedCellFaceGroup newCellAndFaces = Cell::createRefinedCell(H, c0, newCellID);
	cellMap[newCellID] = newCellAndFaces.second;

	// Assign IDs to the new faces and add to the faceMap
	// Also TBD:  Initialize cell x, y, dx, dy values.

	for (auto& f: (*newCellAndFaces.first)) {
		unsigned int newFaceID = provideNewFaceID();
		f->set_id(newFaceID);
		faceMap[newFaceID] = f;
		cout << "Added fid " << newFaceID << endl;
	}

	newCellAndFaces.first->clear();
	delete (newCellAndFaces.first);

	//  c) Extrapolate cell-centered values
	// TBD

	return true;
}

bool Mesh::refineY(Cell *c0) {
	cout << "Y-Refining cell id " << c0->get_id() << endl;

	// Test for ability to refine
	// 1. Check refinement flags
	if (c0->refFlags.test(doCoarsenY) || c0->refFlags.test(wasYRefined))
		return false;

	// 2. Check that neighbouring cells N,S,E,W are not preventing refinement
	Cell *testNbCell = NULL;

	for (auto& d: all_DIR) {
		if (c0->get_nbFaces().at(d)->front()->get_nbCell(testNbCell, d)) // if a neighbour to the North exists
			if (testNbCell->get_lj() == c0->get_lj()-1)                  // if the northern neighbour is twice as wide
				if (refineY(testNbCell) == false)                        // if unable to refine wider northern neighbour
					return false;                                        // give up attempt to refine this cell
	}


	// 3. If above checks are passed, refine cell.
	//    Halve the size of this cell and add another cell to the South or West.
	//  a) First, obtain a new cell identifier

	//  b) Create new cell to the West
	//     Static factory method will be used to initialize, set flags, and remap/create faces
	unsigned int newCellID = provideNewCellID();
	RefinedCellFaceGroup newCellAndFaces = Cell::createRefinedCell(V, c0, newCellID);
	cellMap[newCellID] = newCellAndFaces.second;

	// Assign IDs to the new faces and add to the faceMap
	// Also TBD:  Initialize cell x, y, dx, dy values.

	for (auto& f: (*newCellAndFaces.first)) {
		unsigned int newFaceID = provideNewFaceID();
		f->set_id(newFaceID);
		faceMap[newFaceID] = f;
	}

	newCellAndFaces.first->clear();
	delete (newCellAndFaces.first);

	//  c) Extrapolate cell-centered values
	// TBD
	return true;

}

bool Mesh::coarsenX(Face *f, unsigned int lvl) {





	return true;
}

bool Mesh::coarsenY(Face *f, unsigned int lvl) {



	return true;
}
