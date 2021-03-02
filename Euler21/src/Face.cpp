/*
 * Face.cpp
 *
 *  Created on: Feb. 3, 2021
 *      Author: wakeats
 */

#include "Face.h"
#include "Euler.h"
#include <iostream>

Face::Face() {
	// TODO Auto-generated constructor stub

}

Face::~Face() {
	// TODO Auto-generated destructor stub
	cout << "in destructor, fid: " << id << endl;
}

Face::Face(const Face& f):
		id(f.id),
    length(f.length),
		 x(f.x),
		 y(f.y),
	orient(f.orient)
{
	cout << "in copy constructor, fid: " << id << endl;
}

// Cell *c is requesting that this face know about it
// DIR d:  direction of Face from Cell *c's POV
void Face::connectCell(Cell *c, DIR d) {

	switch (d) {
		case N: nbCells[NEG]=c; break;   // Cell specifying North face; Cell is to the negative side of Face
		case E: nbCells[NEG]=c; break;
		case S: nbCells[POS]=c; break;	// Cell specifying South face; Cell is to the positive side of Face
		case W: nbCells[POS]=c; break;
	}

}

// STATIC
Face* Face::createRefinedFace(ORIENTATION _orient, Face *f0) {

	Face *f1 = new Face(*f0);

	f0->length *= 0.5;
	f1->length *= 0.5;
	f1->orient = _orient;

	cfdFloat halfLength = 0.5*(f0->length);

	switch (_orient) {
	case H:
		f0->x -= halfLength;
		f1->x += halfLength;
		break;
	case V:
		f0->y -= halfLength;
		f1->y += halfLength;
		break;

	}

	return (f1);

}

// 'orient' describes the orientation of the face, not the neighbouring cells
Face* Face::createIntermediateFace(ORIENTATION _orient, Cell *c0, Cell *c1) {

	Face *f0 = new Face();

	f0->orient = _orient;

	switch (_orient) {
	case V:
	{
		f0->connectCell(c0, E);
		f0->connectCell(c1, W);
		break;
	}
	case H:
	{
		f0->connectCell(c0, N);
		f0->connectCell(c1, S);
		break;

	}
	}

	return (f0);
}
