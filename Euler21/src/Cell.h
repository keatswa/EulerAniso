/*
 * Cell.h
 *
 *  Created on: Feb. 3, 2021
 *      Author: wakeats
 */

#ifndef CELL_H_
#define CELL_H_


#include "Euler.h"
#include "Face.h"
#include "Payload.h"
#include <bitset>
#include <iostream>


class Cell {

private:
	friend class Mesh;
	unsigned long id;				// Cell ID
	unsigned int li, lj;			// Refinement level
	unsigned int i_idx, j_idx;		// Integer position index

	cfdFloat x, y;					// Geometric location
	cfdFloat dx, dy;				// Cell dimensions
	cfdFloat errX, errY;			// Error indicators  TBD: vectorize?

	void generate_id() {
//		cout << "id: " << id << " -> " ;
		id = (( ((ulong)i_idx+(li<<0x18)) << 0x20 ) + ((ulong)j_idx+(lj<<0x18)));
//		cout << id << endl;
	}

	NeighbouringFaces nbFaces;		// unordered_map of neighbour faces, indexed by DIR { N, S, W, E }

	PayloadVar *U;

	// Keep a static copy of dt to avoid passing it
	inline static cfdFloat dt = -1;

public:
	Cell();
	Cell(const Cell& c);
	virtual ~Cell();

	inline static void set_dt(cfdFloat _dt) { dt = _dt; }

	PayloadVar* get_U() { return U; }

	cfdFloat get_U(ConservedVariable idx) { return U->get_U(idx); }
	cfdFloat get_PV(PrimitiveVariable idx) { return U->get_PV(idx); }

	cfdFloat get_dU_EuclideanNorm(ConservedVariable idx) {
		return sqrt(pow(U->get_dU(X, idx), 2) + pow(U->get_dU(Y, idx), 2));
	}

	void setCVfromPV(cfdFloat *icPVarr){ U->setCVFromPrimitives(icPVarr); }


	static RefinedCellFaceGroup createRefinedCell(ORIENTATION orient, Cell* c0); //, unsigned int newCellID);


	inline static ulong generate_id(uint _i_idx, uint _j_idx, uint _li, uint _lj) {
		return (( ((ulong)_i_idx+(_li<<0x18)) << 0x20 ) + ((ulong)_j_idx+(_lj<<0x18)));
	}



	bitset<RefinementFlags::NUM_CELLREFINEMENT_FLAGS> refFlags;

	NeighbouringFaces& get_nbFaces() { return(nbFaces); }

	const unsigned long get_id()    const { return(id); }
	const unsigned  int get_i_idx() const { return(i_idx); }
	const unsigned  int get_j_idx() const { return(j_idx); }
	const unsigned  int get_li()    const { return(li); }
	const unsigned  int get_lj()    const { return(lj); }

	cfdFloat get_x() { return(x); }
	cfdFloat get_y() { return(y); }
	cfdFloat get_dx() { return(dx); }
	cfdFloat get_dy() { return(dy); }

	// Get cell dimension {dx,dy} normal to face at given orientation
	cfdFloat get_ds_normal(ORIENTATION orient) {
		switch (orient) {
		case V:
			return dx;
			break;
		case H:
			return dy;
			break;
		}
		return 0;
	}

	void set_id(unsigned long _id) { id = _id; }
	void set_i_idx(unsigned int _i_idx) { i_idx = _i_idx; }
	void set_j_idx(unsigned int _j_idx) { j_idx = _j_idx; }
	void set_li(unsigned int _li) { li = _li;  }
	void set_lj(unsigned int _lj) { lj = _lj; }
	void set_x(cfdFloat _x) { x = _x;  }
	void set_y(cfdFloat _y) { y = _y; }
	void set_dx(cfdFloat _dx) { dx = _dx;  }
	void set_dy(cfdFloat _dy) { dy = _dy; }

	void init_x(cfdFloat _x) { x = _x; }
	void init_y(cfdFloat _y) { y = _y; }
	void init_dx(cfdFloat _dx) { dx = _dx; }
	void init_dy(cfdFloat _dy) { dy = _dy; }
	void init_id() { generate_id(); }


};

#endif /* CELL_H_ */
