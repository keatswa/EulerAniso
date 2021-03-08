/*
 * Mesh.h
 *
 *  Created on: Feb. 3, 2021
 *      Author: wakeats
 */
using namespace std;

#ifndef MESH_H_
#define MESH_H_

#include <iostream>
#include <cmath>
#include <valarray>
#include <vector>
#include <deque>
#include <list>
#include <unordered_map>
#include <forward_list>
#include <functional>
#include <climits>
#include "Euler.h"
#include "Cell.h"
#include "Face.h"
//#include "EulerDisplay.h"

class EulerDisplay;


class Mesh {
private:

	unsigned int maxRefLvl = 0;
	unsigned int iteration;
	cfdFloat maxCFL;

	forward_list<unsigned long> recycledCellIDs;
	forward_list<unsigned long> recycledFaceIDs;

	cfdFloat x0, y0, x1, y1;   // geometric bounds



	EulerDisplay *display;  // Could be a pointer to EulerDisplay's parent class if we have multiple possible 'display' types

	CB_EulerDisplay_drawFn drawFn;


public:
	Mesh();
	virtual ~Mesh();

	// Specify the object that will be used to display results
	void set_display(EulerDisplay *_display) { display = _display; }

	// Set callback to draw function in 'display' object
	void set_cb_drawFn(CB_EulerDisplay_drawFn& df) { drawFn = df; }

	CellMap cellMap;
	FaceMap faceMap;

	void set_iteration(int n) { iteration = n; }

	void set_bounds(cfdFloat _x0, cfdFloat _y0, cfdFloat _x1, cfdFloat _y1) {
		x0 = _x0;
		y0 = _y0;
		x1 = _x1;
		y1 = _y1;
	}

	bool refineX(Cell *c0);
	bool refineY(Cell *c0);
	bool coarsenX(Face *f, unsigned int lvl);
	bool coarsenY(Face *f, unsigned int lvl);

	bool coarsenFace(Face *f, unsigned int lvl, ORIENTATION orient);

	// Refine cells whose x/y refinement level is currently 'reflvl'.  Refine in X and then in Y.
	unsigned int refineCells(unsigned int reflvl);

	// Attempt to coarsen cells currently at 'reflvl' refinement level
	unsigned int coarsenCells(unsigned int reflvl);

	/*
	 * ID handling for cell and face refinement
	 *
	 * - Maintain 2 <forward_list>s, one for recycled Cell IDs and the other for recycled Face IDs.
	 * - Add Cell IDs and Face IDs to these sets when a coarsening operation takes place.
	 * - When refining, if forward_list is not empty, then pop a Cell ID/Face ID out the front.
	 * - If forward_list is empty, then the next available Cell ID should be the size of cellMap.
	 *
	 */


//	unsigned int provideNewCellID() {
//		unsigned int return_id;
//		if (!recycledCellIDs.empty()) {
//			return_id = recycledCellIDs.front();
//			recycledCellIDs.pop_front();
//		}
//		else {
//			return_id = (unsigned int)(cellMap.size());
//		}
//		return(return_id);
//	}
//
//	unsigned int provideNewFaceID() {
//		unsigned int return_id;
//		if (!recycledFaceIDs.empty()) {
//			return_id = recycledFaceIDs.front();
//			recycledFaceIDs.pop_front();
//		}
//		else {
//			return_id = (unsigned int)(faceMap.size());
//		}
//		return(return_id);
//	}



//	ulong provideNewCellID(ulong i, ulong li, ulong j, ulong lj) {
//		ulong return_id = 0;
//		return_id =  ( (i*(2<<li)) << 0x20 ) + (j*(2<<lj));
//		return return_id;
//	}

	ulong provideNewFaceID(cfdFloat x, cfdFloat y) {
		ulong return_id = 0;
		return_id =  ( (ulong)(UINT_MAX*((x-x0)/(x1-x0))) << 0x1F ) + (ulong)(UINT_MAX*((y-y0)/(y1-y0)));
		return return_id;
	}






	void addAndRemapFaceToCell(Cell *c0, Cell *c1, DIR d);

	// Refine from level 0 to target 'reflvl'
	void doUniformRefine(unsigned int reflvl);

	// Attempt to coarsen cells to reach 'reflvl' refinement level
	void doUniformCoarsen(unsigned int reflvl);

	// Mesh initialization - change name to indicate for initialization only
	int init_addCell(const Cell& c) {
		cellMap.emplace(c.get_id(), new Cell(c));
		return(0);
	}

	// Mesh initialization - change name to indicate for initialization only
	int init_addFace(const Face& f) {
		faceMap.emplace(f.get_id(), new Face(f));
		return(0);
	}

	// Mesh initialization - Assume Cell's nbFaces has been constructed/initialized
	void init_mapFaceToCell(unsigned long cellID, unsigned long faceID, DIR d) {

//		Cell *tmpCell = cellMap.at(cellID);

//		cout << "map face " << faceID << " to tgt cell id: " << cellMap.at(cellID)->get_id() << endl;
//		cout << "unordered_map NeighbouringFaces size: " << cellMap.at(cellID)->get_nbFaces().size() << endl;
//		cellMap.at(cellID)->get_nbFaces().emplace(d, new FaceDeque());                         // TBD CAN WE INITIALIZE IN THE DEFAULT CONSTRUCTOR?

		if (cellMap.at(cellID)->get_nbFaces().at(d)->size() > 0)
			if (cellMap.at(cellID)->get_nbFaces().at(d)->back()->get_id() == faceID)
				cellMap.at(cellID)->get_nbFaces().at(d)->pop_back();

		cellMap.at(cellID)->get_nbFaces().at(d)->push_back(faceMap.at(faceID));

		// Now provide Face objects with a reference to each Cell
		cellMap.at(cellID)->get_nbFaces().at(d)->back()->connectCell(cellMap.at(cellID), d);

//		cout << "gave " << cellMap.at(cellID)->get_nbFaces().at(d)->front()->get_id() << " a reference to cell " << cellID << endl;;

		if (d == N) {
			cellMap.at(cellID)->init_x(faceMap.at(faceID)->get_x());
			cellMap.at(cellID)->init_dx(faceMap.at(faceID)->get_length());
		}

		if (d == E) {
			cellMap.at(cellID)->init_y(faceMap.at(faceID)->get_y());
			cellMap.at(cellID)->init_dy(faceMap.at(faceID)->get_length());
		}

		if (d == N || d == S)
			faceMap.at(faceID)->set_orient(H);
		else
			faceMap.at(faceID)->set_orient(V);

	}

// DEPRECATE
//	void init_cell_and_face_IDs() {
//
//		for (auto& f: faceMap){
//			Face *tf = f.second;
//			auto nh = faceMap.extract(f.first);
//			tf->set_id(provideNewFaceID(tf->get_x(), tf->get_y()));
//			nh.key() = tf->get_id();
//			faceMap.insert(move(nh));
//		}
//
//		for (auto& c: cellMap) {
//			Cell *tc = c.second;
//			auto nh = cellMap.extract(c.first);
//			tc->init_id();
//			nh.key() = tc->get_id();
//			cellMap.insert(move(nh));
//		}
//	}



	void printMesh() {

//		for (auto& x: cellMap){
//			cout << x.first << ": " << x.second->get_i_idx() << " , " << x.second->get_j_idx() << endl;
//		}
//
//		for (auto& x: cellMap){
//
//			for (auto& d : all_DIR) {
//				cout << "Cell: " << x.first << ": " << d << " face: " << x.second->get_nbFaces().at(d)->back()->get_id() << endl;
//			}
//		}

		for (auto& x: faceMap){
			Cell *tmpCell = NULL;
			cout << "fid " << x.first << ": " << x.second->get_id() <<
					", length: " << x.second->get_length() <<
					", orient: " << (x.second->get_orient() == H ? "H" : "V") <<  endl;
			if (x.second->get_nbCell(tmpCell, NEG))
				cout << "  neg cid: " << tmpCell->get_id() << endl;
			if (x.second->get_nbCell(tmpCell, POS))
				cout << "  pos cid: " << tmpCell->get_id() << endl;

			cout << "BC Type: " << x.second->get_bcType() << endl;

		}


		// HOWTO ERASE CELLS
//		while (cellMap.size() > 0) {
//			cout << "erasing cell " << cellMap.begin()->second->get_id() << endl;
//			delete(cellMap.begin()->second);
//			cellMap.erase(cellMap.begin());
//		}

	}


	void calcFaceFluxes();

	void doCellTimestep (const cfdFloat& dt, int rkStage );

	cfdFloat getMaxCFL() { return maxCFL; }

};

#endif /* MESH_H_ */
