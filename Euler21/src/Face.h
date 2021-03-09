/*
 * Face.h
 *
 *  Created on: Feb. 3, 2021
 *      Author: wakeats
 */

#ifndef FACE_H_
#define FACE_H_

#include "Euler.h"
#include <bitset>


class Face {

private:
	friend class Cell;
	// network topology
	unsigned long id;				// Face ID
	NeighbouringCells nbCells;  	// mapped by [NEG, POS]

	// geometry
	cfdFloat length;
	cfdFloat x, y;
	ORIENTATION orient;

	// payload
	BCType bc_type;
	bool is_bc;




public:
	Face();
	Face(const Face& f);
	virtual ~Face();

	void set_bcType(BCType _bc) { bc_type = _bc; }
	BCType get_bcType() { return bc_type; }
	void set_is_bc(bool _bc) { is_bc = _bc; }
	bool get_is_bc() { return is_bc; }

	bitset<FaceRefinementFlags::NUM_FACEREFINEMENT_FLAGS> faceRefFlags;


	void set_id(unsigned long _id) { id = _id; }
	const unsigned long get_id() const { return(id); }

	void set_length(cfdFloat _length) { length = _length; }
	const cfdFloat get_length() const { return(length); }

	ORIENTATION get_orient() { return(orient); }
	void set_orient(ORIENTATION _orient) { orient = _orient; }

	cfdFloat get_x() { return(x); }
	cfdFloat get_y() { return(y); }

	void set_x(cfdFloat _x) { x = _x; }
	void set_y(cfdFloat _y) { y = _y; }

	void extend() {
		if (orient == H)
			x += 0.5*length;
		else
			y += 0.5*length;

		length *= 2.0;
	}

	void connectCell(Cell* c, DIR d);

	// Pass a reference to the pointer in order to guarantee
	// being able to redirect its pointer address

	bool get_nbCell(Cell*& c, DIR d) {
		c = NULL;
		if (d == N || d == E) {
			if (nbCells.find(POS) == nbCells.end())
				return(false);
			else {
				c = nbCells.find(POS)->second;
				return(true);
			}
		}
		else {
			if (nbCells.find(NEG) == nbCells.end())
				return(false);
			else {
				c = nbCells.find(NEG)->second;
				return(true);
			}
		}
		return(false);
	}


	bool get_nbCell(Cell*& c, SIGN sgn) {
		c = NULL;
		if (nbCells.find(sgn) == nbCells.end())
			return(false);
		else {
			c = nbCells.find(sgn)->second;
			return(true);
		}
		return(false);
	}



	static Face* createRefinedFace(ORIENTATION orient, Face *f0);

	static Face* createIntermediateFace(ORIENTATION orient, Cell *c0, Cell *c1);

};

#endif /* FACE_H_ */
