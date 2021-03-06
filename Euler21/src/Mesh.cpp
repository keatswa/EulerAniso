/*
 * Mesh.cpp
 *
 *  Created on: Feb. 3, 2021
 *      Author: wakeats
 */

#include <functional>
#include <algorithm>
#include "Mesh.h"
#include <vector>
#include <iostream>
#include <climits>

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
		fitor++;
	}

	numCoarsenedAtLvl = coarsenCells(reflvl);

	cout << "Mesh::doUniformCoarsen: Coarsened " << numCoarsenedAtLvl << " Faces." << endl;

}



// Refine cells whose x/y refinement level is reflvl.  Refine in X and then in Y.
unsigned int Mesh::refineCells(unsigned int reflvl) {

	unsigned int numRefined = 0;
	unsigned int xRefined = 0;
	unsigned int yRefined = 0;
	unsigned int i;
	vector<unsigned long> *cellIDVec = new vector<unsigned long>();

	// Define and populate vector with cell identifiers present in cellMap.
	// We will operate only on these IDs since new ones will be added to the map.



	// refine in Y
//	cout << "Cells to Y-refine: " << endl;
	cellIDVec->resize(cellMap.size());
	i = 0;
	for (auto& c: cellMap) {
		(*cellIDVec)[i] = c.first;
//		cout << (*cellIDVec)[i] << " , " ;
		i++;
	}
//	cout << endl;

	for (auto c: (*cellIDVec)) {
		Cell *tgtCell = cellMap[c];
//		Cell *tgtCell = cellMap[(*cellIDVec)[c]];

		if (tgtCell->refFlags.test(doRefineY) && (tgtCell->get_lj() == reflvl)){
			if (refineY(tgtCell)) {
				yRefined++;
			}
		}

//		((display)->*(drawFn))(cellMap, faceMap);
// Alternatively, but requires g++-10:
//		invoke(drawFn, display, cellMap, faceMap); // @suppress("Invalid arguments")

	}
	cout << "Y-refined " << yRefined << " cells." << endl;

	cout << "cellMap.size() == " << cellMap.size() << endl;

	// refine in X
//	cout << "Cells to X-refine: " << endl;
	cellIDVec->resize(cellMap.size());
	i = 0;
	for (auto& c: cellMap) {
		(*cellIDVec)[i] = c.first;
//		cout << (*cellIDVec)[i] << " , " ;
		i++;
	}
//	cout << endl;

	for (auto c: (*cellIDVec)) {

//		cout << "c: " << c << endl;

//		Cell *tgtCell = cellMap[(*cellIDVec)[c]];
		Cell *tgtCell = cellMap[c];

		if (tgtCell->refFlags.test(doRefineX) && (tgtCell->get_li() == reflvl)){
			if (refineX(tgtCell)) {
				xRefined++;
			}
		}
//		((display)->*(drawFn))(cellMap, faceMap);
//		std::invoke(drawFn, display, cellMap, faceMap); // @suppress("Invalid arguments")

	}

	cout << "X-refined " << xRefined << " cells." << endl;
	cout << "cellMap.size() == " << cellMap.size() << endl;

	numRefined = xRefined + yRefined;

	// Clean up
	cellIDVec->clear();
	delete(cellIDVec);
	cout << "Drawing " << cellMap.size() << " cells, " << faceMap.size() << " faces." << endl;
	((display)->*(drawFn))(cellMap, faceMap);

	cout << "Refined " << numRefined << " cells." << endl;
	return numRefined;
}

// Attempt to coarsen cells currently at 'reflvl' refinement level
unsigned int Mesh::coarsenCells(unsigned int reflvl) {

	unsigned int numCoarsened = 0;
	unsigned int numCoarsenedX = 0;
	unsigned int numCoarsenedY = 0;
	forward_list<unsigned long>::iterator recycledFaceIter; // = nullptr;

	// Coarsen in Y

	for (auto& f: faceMap) {
		Face *tgtFace = f.second;

		if (tgtFace->faceRefFlags.test(doCoarsen)
				&& ( tgtFace->get_orient() == H)
				&& (!tgtFace->faceRefFlags.test(doRecycleFace)) ) {
//			if (coarsenX(tgtFace, reflvl)) {
			if (coarsenFace(tgtFace, reflvl, H)) {
				tgtFace->faceRefFlags.reset(doCoarsen);
				tgtFace->faceRefFlags.set(doRecycleFace);
				numCoarsened++;
			}
		}
//		((display)->*(drawFn))(cellMap, faceMap);
	}

	numCoarsenedY = numCoarsened;
	cout << "Y sweep:  Coarsened " << numCoarsenedY << " Faces." << endl;
	cout << "cellMap.size() == " << cellMap.size() << endl;
	cout << "faceMap.size() == " << faceMap.size() << endl;

	recycledFaceIter = recycledFaceIDs.begin();
	while (recycledFaceIter != recycledFaceIDs.end()) {
		if (faceMap.find(*recycledFaceIter) != faceMap.end()) {
			delete faceMap[*recycledFaceIter];
			faceMap.erase(*recycledFaceIter);
		}
		recycledFaceIter++;
	}
	recycledFaceIDs.clear();

	// Coarsen in X

	for (auto& f: faceMap) {
		Face *tgtFace = f.second;

		if (tgtFace->faceRefFlags.test(doCoarsen)
				&& ( tgtFace->get_orient() == V)
				&& (!tgtFace->faceRefFlags.test(doRecycleFace)) ) {
//			if (coarsenX(tgtFace, reflvl)) {
			if (coarsenFace(tgtFace, reflvl, V)) {
				tgtFace->faceRefFlags.reset(doCoarsen);
				tgtFace->faceRefFlags.set(doRecycleFace);
				numCoarsened++;
			}
		}
//		((display)->*(drawFn))(cellMap, faceMap);
	}

	numCoarsenedX = numCoarsened - numCoarsenedY;
	cout << "X sweep:  Coarsened " << numCoarsenedX << " Faces." << endl;
	cout << "cellMap.size() == " << cellMap.size() << endl;
	cout << "faceMap.size() == " << faceMap.size() << endl;


	recycledFaceIter = recycledFaceIDs.begin();
	while (recycledFaceIter != recycledFaceIDs.end()) {
		if (faceMap.find(*recycledFaceIter) != faceMap.end()) {
			delete faceMap[*recycledFaceIter];
			faceMap.erase(*recycledFaceIter);
		}
		recycledFaceIter++;
	}
	recycledFaceIDs.clear();




//	forward_list<unsigned int>::iterator recycledFaceIter = recycledFaceIDs.begin();
//	while (recycledFaceIter != recycledFaceIDs.end()) {
//		if (faceMap.find(*recycledFaceIter) != faceMap.end()) {
//			delete faceMap[*recycledFaceIter];
//			faceMap.erase(*recycledFaceIter);
//		}
//		recycledFaceIter++;
//	}


//	unsigned int i;
//	vector<unsigned int> *faceIDVec = new vector<unsigned int>();
//
//	// Define and populate vector with cell identifiers present in cellMap.
//	// We will operate only on these IDs since new ones will be added to the map.
//	faceIDVec->resize(faceMap.size());
//	i = 0;
//	for (auto& f: faceMap) {
//		(*faceIDVec)[i] = f.first;
//		i++;
//	}
//
//	for (auto f: (*faceIDVec)) {
//		Face *tgtFace = faceMap[(*faceIDVec)[f]];
//
//		if (tgtFace->faceRefFlags.test(doCoarsen)
//				&& ( tgtFace->get_orient() == V)
//				&& (!tgtFace->faceRefFlags.test(doRecycleFace)) ) {
//			if (coarsenX(tgtFace, reflvl)) {
//				tgtFace->faceRefFlags.reset(doCoarsen);
//				tgtFace->faceRefFlags.set(doRecycleFace);
//				numCoarsened++;
//			}
//		}
//	}
	cout << "Mesh::coarsenCells: Coarsened " << numCoarsenedX + numCoarsenedY << " Faces." << endl;

	cout << "Drawing " << cellMap.size() << " cells, " << faceMap.size() << " faces." << endl;
	((display)->*(drawFn))(cellMap, faceMap);


	return numCoarsened;
}


bool Mesh::refineX(Cell *c0) {

//	cout << "X-Refining cell id " << c0->get_id() << endl;

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
	RefinedCellFaceGroup newCellAndFaces = Cell::createRefinedCell(H, c0); //, newCellID);

	unsigned long newCellID = newCellAndFaces.second->get_id();

	if (cellMap.find(newCellID) != cellMap.end())
		cout << "X: Squashing cid " << newCellID << endl;

	cellMap[newCellID] = newCellAndFaces.second;

	// Assign IDs to the new faces and add to the faceMap
	// Also TBD:  Initialize cell x, y, dx, dy values.

	for (auto& f: (*newCellAndFaces.first)) {
		unsigned long newFaceID = provideNewFaceID(f->get_x(), f->get_y());
		f->set_id(newFaceID);

		if (faceMap.find(newFaceID) != faceMap.end())
			cout << "X: Squashing fid " << newFaceID << endl;

		faceMap[newFaceID] = f;
//		cout << "Added fid " << newFaceID << endl;
	}

	newCellAndFaces.first->clear();
	delete (newCellAndFaces.first);

	//  c) Extrapolate cell-centered values
	// TBD

	return true;
}

bool Mesh::refineY(Cell *c0) {
//	cout << "Y-Refining cell id " << c0->get_id() << endl;

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
	RefinedCellFaceGroup newCellAndFaces = Cell::createRefinedCell(V, c0); //, newCellID);
	unsigned long newCellID = newCellAndFaces.second->get_id();
	if (cellMap.find(newCellID) != cellMap.end())
		cout << "Y: Squashing cid " << newCellID << endl;
	cellMap[newCellID] = newCellAndFaces.second;

	// Assign IDs to the new faces and add to the faceMap
	// Also TBD:  Initialize cell x, y, dx, dy values.

	for (auto& f: (*newCellAndFaces.first)) {
		unsigned long newFaceID = provideNewFaceID(f->get_x(), f->get_y());
		f->set_id(newFaceID);
		if (faceMap.find(newFaceID) != faceMap.end())
			cout << "Y: Squashing fid " << newFaceID << endl;
		faceMap[newFaceID] = f;
	}

	newCellAndFaces.first->clear();
	delete (newCellAndFaces.first);

	//  c) Extrapolate cell-centered values
	// TBD
	return true;

}


bool Mesh::coarsenFace(Face *f, unsigned int lvl, ORIENTATION orient) {

	Cell *c[2] = {NULL, NULL};   // [0] - c0 - NEG;  [1] - c1 - POS
	unsigned int i_idx[2] = {UINT_MAX, UINT_MAX};
	unsigned int j_idx[2] = {UINT_MAX, UINT_MAX};
	unsigned int li[2]    = {UINT_MAX, UINT_MAX};
	unsigned int lj[2]    = {UINT_MAX, UINT_MAX};

	unsigned int *idx = NULL;
	unsigned int *lev = NULL;

	std::initializer_list<DIR> xy_DIR;
	DIR axisDir[2];

	switch (orient) {
		case V:      // Vertical face -> Horizontal coarsening
		{
			idx = i_idx;   // Cell index points to i-indices
			lev = li;
			xy_DIR = std::initializer_list<DIR>(y_DIR);
			axisDir[0] = W;
			axisDir[1] = E;
//			axisDir = {W,E};
			break;
		}
		case H:		// Horizontal face -> Vertical coarsening
		{
			idx = j_idx;	// Cell index points to j-indices
			lev = lj;
			xy_DIR = std::initializer_list<DIR>(x_DIR);
//			xy_DIR = x_DIR;
			axisDir[0] = S;
			axisDir[1] = N;
//			axisDir = {S,N};
			break;
		}
	}


	for (auto& s: all_SIGN) {

		if (!(f->get_nbCell(c[s], s))) {
			f->faceRefFlags.reset(doCoarsen);
			return false;   // No cell exists; must be a boundary face.  Cannot coarsen.
		} else {

			i_idx[s] = c[s]->get_i_idx();
			j_idx[s] = c[s]->get_j_idx();
			li[s] = c[s]->get_li();
			lj[s] = c[s]->get_lj();

			// Abort if neg cell i-index or j-index not even
			if (idx[0] % 2 > 0)    // points to i or j index
				return false;

			// Check neighbouring cell refinement levels
			// a) North and South if coarsening horizontally
			for (auto& d: xy_DIR) {
				if (c[s]->nbFaces[d]->size() > 1)
					return false;
			}
			// b) West cell(s)
			Cell *tmpCell = NULL;

			for_each(c[s]->nbFaces.at(axisDir[s])->begin(), c[s]->nbFaces.at(axisDir[s])->end(),
					[&tmpCell,lvl,orient,s](Face *f) -> bool { // Lambda fn
				if (f->get_nbCell(tmpCell, s)) {
					if (orient == V) {
						if (tmpCell->get_li() > lvl) {
							return false;
						}
					} else {
						if (tmpCell->get_lj() > lvl) {
							return false;
						}
					}
				}
				return true;
			} );   // end for_each
		} // end if nbCell
	}  // end for SIGN


	// Abort if current refinement level not at target level
	if (lev[0] != lvl) {
		return false;
	}

	// Abort if dx not equal
	if (li[0] != li[1]) {
		return false;
	}

	// Abort if dy not equal
	if (lj[0] != lj[1]) {
		return false;
	}


	// Finally, Coarsen:


	Cell *c0nb = NULL;
	Cell *c1nb = NULL;

	// 1. Supply c0 with c1's faces
	for (auto& d: xy_DIR) {
		// Case 1:  c0 and c1 share a single cell neighbour in direction 'd'.
		//          -> disconnect cell neighbour from c1's N face
		//			-> recycle c1's N face, and extend c0's N face

		bool c0_has_nb = c[0]->nbFaces[d]->back()->get_nbCell(c0nb, d);
		bool c1_has_nb = c[1]->nbFaces[d]->back()->get_nbCell(c1nb, d);


		if (c0_has_nb && c1_has_nb) {                   // c0 and c1 both have a neighbour in direction 'd'
			if (c0nb->get_id() == c1nb->get_id()) {        // Bingo - c0 and c1 share a single neighbour

				c1nb->nbFaces.at(oppositeDir(d))->pop_back();   // disconnect cell neighbour from c1's face

				Face *f_toRecycle = c[1]->nbFaces[d]->back();
				f_toRecycle->faceRefFlags.set(doRecycleFace);
				recycledFaceIDs.push_front(f_toRecycle->get_id());
//				faceMap.erase(f_toRecycle->get_id()); // will invalidate faceMap iterator
//				delete f_toRecycle;
				c[0]->nbFaces[d]->back()->extend();



			}
			else {
	// Case 2:  c0 and c1 have separate cell neighbours in direction 'd'.
	//          -> reassign c1's N face to c0.  Connect c0 to face and connect face to c0.

				c[1]->nbFaces[d]->back()->connectCell(c[0], d);
				c[0]->nbFaces[d]->push_back(c[1]->nbFaces[d]->back());


			}
		}
		else if (!c0_has_nb && !c1_has_nb) {             // Neither cell has a neighbour.  Must lie on a boundary.

			Face *f_toRecycle = c[1]->nbFaces[d]->back();
			f_toRecycle->faceRefFlags.set(doRecycleFace);
			recycledFaceIDs.push_front(f_toRecycle->get_id());
//			faceMap.erase(f_toRecycle->get_id()); // will invalidate faceMap iterator
//			delete f_toRecycle;
			c[0]->nbFaces[d]->back()->extend();

		}
		else {
			cout << "Coarsening ERROR: one of cells " << c[0]->get_id() << " , " << c[1]->get_id() << " has no neighbour!" << endl;
			return false;
		}

	}

	// 2. Map c1's E faces to c0
	swap(*(c[0]->nbFaces.at(axisDir[1])), *(c[1]->nbFaces.at(axisDir[1])));
	// c0 now has c1's E faces.  Go through these faces and connect c0 to them.
	for_each(c[0]->nbFaces.at(axisDir[1])->begin(), c[0]->nbFaces.at(axisDir[1])->end(),
			[&c,axisDir](Face *testNbFace) { testNbFace->connectCell(c[0], axisDir[1]); });


	// 3. Geometric parameters
	if (orient == V) {
		c[0]->set_x(c[0]->get_x() + 0.5*c[0]->get_dx());
		c[0]->set_dx(c[0]->get_dx()*2);
		c[0]->set_li(lvl-1);
		c[0]->set_i_idx(i_idx[0]/2);
	} else {
		c[0]->set_y(c[0]->get_y() + 0.5*c[0]->get_dy());
		c[0]->set_dy(c[0]->get_dy()*2);
		c[0]->set_lj(lvl-1);
		c[0]->set_j_idx(j_idx[0]/2);
	}





	// 4.  Recycle c1 and f
	c[1]->refFlags.set(doRecycleCell);
	f->faceRefFlags.set(doRecycleFace);

//	cout << "Coarsen, recycle cid " << c[1]->get_id()  << endl;
	recycledCellIDs.push_front(c[1]->get_id());
	cellMap.erase(c[1]->get_id());
	delete c[1];


	recycledFaceIDs.push_front(f->get_id());
//	faceMap.erase(f->get_id());					// Erasing will invalidate iterator
//	delete f;



	return true;
}




bool Mesh::coarsenX(Face *f, unsigned int lvl) {

	// Face is vertical.  Check cells on either side:
	Cell *c0 = NULL;  // Cell on the NEG side
	Cell *c1 = NULL;  // Cell on the POS side
	// TBD: Simplify using array of [NEG,POS]
	unsigned int i_neg = UINT_MAX;
	unsigned int li_neg = UINT_MAX;
	unsigned int lj_neg = UINT_MAX;
//	unsigned int i_pos = UINT_MAX;
	unsigned int li_pos = UINT_MAX;
	unsigned int lj_pos = UINT_MAX;


	// Abort if face is a BC
	if (f->get_nbCell(c0, NEG)) {   // If a cell exists to the neg side (West if X-refining, South if Y-refining)
		i_neg = c0->get_i_idx();
		li_neg = c0->get_li();
		lj_neg = c0->get_lj();

		// Abort if neg cell i-index not even
		if (i_neg % 2 > 0)
			return false;

		// Check neighbouring cell refinement levels
		// a) North and South
		for (auto& d: y_DIR) {
			if (c0->nbFaces[d]->size() > 1)
				return false;
		}
		// b) West cell(s)
		Cell *tmpCell = NULL;

		for_each(c0->nbFaces.at(W)->begin(), c0->nbFaces.at(W)->end(),
				[&tmpCell,lvl](Face *f) -> bool { // Lambda fn
			if (f->get_nbCell(tmpCell, NEG)) {
				if (tmpCell->get_li() > lvl) {
					return false;
				}
			}
			return true;
		} );   // end for_each

//		for (auto& f: tmpCellNeg->nbFaces.at(W)) {
//			if (f->get_nbCell(tmpCell, NEG)) {
//				if (tmpCell->get_li() > lvl)
//					return false;
//			}
//		}



	} else {
		f->faceRefFlags.reset(doCoarsen);
		return false;   // No cell exists; must be a boundary face.  Cannot x-coarsen.
	}

	// Abort if face is a BC
	if (f->get_nbCell(c1, POS)) {
//		i_pos = c1->get_i_idx();
		li_pos = c1->get_li();
		lj_pos = c1->get_lj();

		// Similar to code from above..
		// Check neighbouring cell refinement levels
		// a) North and South
		for (auto& d: y_DIR) {
			if (c1->nbFaces[d]->size() > 1)
				return false;
		}
		// b) East cell(s)
		Cell *tmpCell = NULL;
		for_each(c1->nbFaces.at(E)->begin(), c1->nbFaces.at(E)->end(),
				[&tmpCell,lvl](Face *f) -> bool { // Lambda fn
			if (f->get_nbCell(tmpCell, POS)) {
				if (tmpCell->get_li() > lvl) {
					return false;
				}
			}
			return true;
		} );   // end for_each

//		for (auto& f: tmpCellPos->nbFaces.at(E)) {
//			if (f->get_nbCell(tmpCell, POS)) {
//				if (tmpCell->get_li() > lvl)
//					return false;
//			}
//		}
	} else {
		return false;   // No cell exists; must be a boundary face.  Cannot x-coarsen.
	}

	// Abort if dx not equal
	if (li_neg != li_pos) {
		return false;
	}

	// Abort if current refinement level not at target level
	if (li_neg != lvl) {
		return false;
	}

	// Abort if dy not equal
	if (lj_neg != lj_pos) {
		return false;
	}


	// Finally, Coarsen:


	Cell *c0nb = NULL;
	Cell *c1nb = NULL;

	// 1. Supply c0 with c1's faces
	for (auto& d: y_DIR) {
		// Case 1:  c0 and c1 share a single cell neighbour in direction 'd'.
		//          -> disconnect cell neighbour from c1's N face
		//			-> recycle c1's N face, and extend c0's N face

		bool c0_has_nb = c0->nbFaces[d]->back()->get_nbCell(c0nb, d);
		bool c1_has_nb = c1->nbFaces[d]->back()->get_nbCell(c1nb, d);


		if (c0_has_nb && c1_has_nb) {                   // c0 and c1 both have a neighbour in direction 'd'
			if (c0nb->get_id() == c1nb->get_id()) {        // Bingo - c0 and c1 share a single neighbour

				c1nb->nbFaces.at(oppositeDir(d))->pop_back();   // disconnect cell neighbour from c1's face

				Face *f_toRecycle = c1->nbFaces[d]->back();
				f_toRecycle->faceRefFlags.set(doRecycleFace);
				recycledFaceIDs.push_front(f_toRecycle->get_id());
//				faceMap.erase(f_toRecycle->get_id()); // will invalidate faceMap iterator
//				delete f_toRecycle;
				c0->nbFaces[d]->back()->extend();



			}
			else {
	// Case 2:  c0 and c1 have separate cell neighbours in direction 'd'.
	//          -> reassign c1's N face to c0.  Connect c0 to face and connect face to c0.

				c1->nbFaces[d]->back()->connectCell(c0, d);
				c0->nbFaces[d]->push_back(c1->nbFaces[d]->back());


			}
		}
		else if (!c0_has_nb && !c1_has_nb) {             // Neither cell has a neighbour.  Must lie on a boundary.

			Face *f_toRecycle = c1->nbFaces[d]->back();
			f_toRecycle->faceRefFlags.set(doRecycleFace);
			recycledFaceIDs.push_front(f_toRecycle->get_id());
//			faceMap.erase(f_toRecycle->get_id()); // will invalidate faceMap iterator
//			delete f_toRecycle;
			c0->nbFaces[d]->back()->extend();

		}
		else {
			cout << "Coarsening ERROR: one of cells " << c0->get_id() << " , " << c1->get_id() << " has no neighbour!" << endl;
			return false;
		}

	}

	// 2. Map c1's E faces to c0
	swap(*(c0->nbFaces.at(E)), *(c1->nbFaces.at(E)));
	// c0 now has c1's E faces.  Go through these faces and connect c0 to them.
	for_each(c0->nbFaces.at(E)->begin(), c0->nbFaces.at(E)->end(),
			[&c0](Face *testNbFace) { testNbFace->connectCell(c0, E); });


	// 3. Geometric parameters
	c0->set_x(c0->get_x() + 0.5*c0->get_dx());
	c0->set_dx(c0->get_dx()*2);
	c0->set_li(lvl-1);
	c0->set_i_idx(i_neg/2);


	// 4.  Recycle c1 and f
	c1->refFlags.set(doRecycleCell);
	f->faceRefFlags.set(doRecycleFace);

	recycledCellIDs.push_front(c1->get_id());
	cellMap.erase(c1->get_id());
	delete c1;


	recycledFaceIDs.push_front(f->get_id());
//	faceMap.erase(f->get_id());					// Erasing will invalidate iterator
//	delete f;



	return true;
}

bool Mesh::coarsenY(Face *f, unsigned int lvl) {



	return true;
}
