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
#include <bitset>


class Cell {

private:
	unsigned int id;				// Cell ID
	unsigned int li, lj;			// Refinement level
	unsigned int i_idx, j_idx;		// Integer position index

	cfdFloat x, y;					// Geometric location
	cfdFloat dx, dy;				// Cell dimensions
	cfdFloat errX, errY;			// Error indicators  TBD: vectorize?




public:
	Cell();
	Cell(const Cell& c);
	virtual ~Cell();

	static RefinedCellFaceGroup createRefinedCell(ORIENTATION orient, Cell* c0, unsigned int newCellID);
	NeighbouringFaces nbFaces;		// unordered_map of neighbour faces, indexed by DIR { N, S, W, E }

	bitset<RefinementFlags::NUM_CELLREFINEMENT_FLAGS> refFlags;

	NeighbouringFaces& get_nbFaces() { return(nbFaces); }

	const unsigned int get_id()    const { return(id); }
	const unsigned int get_i_idx() const { return(i_idx); }
	const unsigned int get_j_idx() const { return(j_idx); }
	const unsigned int get_li()    const { return(li); }
	const unsigned int get_lj()    const { return(lj); }

	cfdFloat get_x() { return(x); }
	cfdFloat get_y() { return(y); }
	cfdFloat get_dx() { return(dx); }
	cfdFloat get_dy() { return(dy); }

	void set_id(unsigned int _id) { id = _id; }
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

};

#endif /* CELL_H_ */
