/*
 * Cell.cpp
 *
 *  Created on: Feb. 3, 2021
 *      Author: wakeats
 */

#include "Cell.h"
#include <iostream>
#include <deque>
#include <functional>
#include <algorithm>

DIR oppositeDir(DIR d) {
	DIR oppDir;

	switch (d) {
	case N:
		oppDir = S;
		break;
	case E:
		oppDir = W;
		break;
	case S:
		oppDir = N;
		break;
	case W:
		oppDir = E;
		break;
	}
	return(oppDir);
}


Cell::Cell() {
//	cout << "in constructor, cid: " << id << endl;
//	cout << "nbFaces.size(): " << nbFaces.size() << endl;
}


//Initialize Cell from a copy but do not connect faces
Cell::Cell(const Cell& c):
   id(c.id),
   li(c.li),
   lj(c.lj),
i_idx(c.i_idx),
j_idx(c.j_idx),
   x(c.x),
   y(c.y),
  dx(c.dx),
  dy(c.dy),
errX(c.errX),
errY(c.errY),
refFlags(c.refFlags)
{

//	cout << "in copy constructor, cid: " << id << endl;
	for (auto d : all_DIR) {
		// Do job with e
		nbFaces.emplace(d, new FaceDeque());
//		nbFaces[d] = new FaceDeque();
	}

}



Cell::~Cell() {

//	cout << "in destructor, cid: " << id << endl;

	if (nbFaces.size() > 0)
	for (auto d : all_DIR) {
//		cout << "Deleting neighbour dir " << d << ": fid: " << nbFaces.at(d)->front()->get_id() << endl;
		nbFaces.at(d)->clear();
		delete(nbFaces.at(d));
	}


}

//STATIC
RefinedCellFaceGroup Cell::createRefinedCell(ORIENTATION orient, Cell *c0, unsigned int newCellID) {

	FaceDeque *newFaces = new FaceDeque();

	Cell *c1 = new Cell(*c0);
	c1->id = newCellID;

	c0->refFlags.reset(doRecycleCell);
	c1->refFlags.reset(doRecycleCell);

	switch (orient) {
	case H:
	{
		// Reset refinement flags
		c0->refFlags.reset(doRefineX);
		c1->refFlags.reset(doRefineX);

		c0->refFlags.reset(doCoarsenX);
		c1->refFlags.reset(doCoarsenX);

		c0->refFlags.set(wasXRefined);
		c1->refFlags.set(wasXRefined);

		// Alter cell geometry
		c0->dx *= 0.5;
		c1->dx *= 0.5;
		c0->x -= 0.5*c0->dx;
		c1->x += 0.5*c1->dx;
		c0->li++;
		c1->li++;
		c0->i_idx *= 2;
		c1->i_idx = c0->i_idx+1;

		// Look at North and South faces:
		//  If cell has only one face to the north/south, split that face in two.
		// Create the new face to the right (back of the deque) of the initial face.
		for (auto& d: y_DIR) {   // d in {N, S}

//			cout << "dir: " << d << ", # faces: " << c0->nbFaces.at(d)->size() << endl;

			if (c0->nbFaces.at(d)->size() == 1) {

				Face *f = Face::createRefinedFace(orient, c0->nbFaces.at(d)->front());
				c1->nbFaces.at(d)->push_back(f);
				newFaces->push_back(f);
				f->connectCell(c1, d);

				// If there exists a cell to the north/south, connect that cell to the newly created faces.
				Cell *testNbCell = NULL;
				if (c0->nbFaces.at(d)->front()->get_nbCell(testNbCell, d)) {
					// Connect Face to Cell
					f->connectCell(testNbCell, oppositeDir(d));
					// Connect Cell to Face
					testNbCell->nbFaces.at(oppositeDir(d))->push_back(f);
				}
			}
			//  If cell already has two faces to the north, reassign rightmost face to c1.
			else {
				// Connect Face to Cell c1
				c0->nbFaces.at(d)->back()->connectCell(c1, d);
				// Connect Cell c1 to Face
				c1->nbFaces.at(d)->push_back(c0->nbFaces.at(d)->back());
				// Disconnect Cell c0 from Face
				c0->nbFaces.at(d)->pop_back();
				// no need to adjust the 2nd face's N neighbour, it should already be mapped.
			}
		}
		// Connect c0's East face(s) to c1
//		c0->nbFaces.at(E)->front()->connectCell(c1, E);
//		if (c0->nbFaces.at(E)->size() > 1)
//			c0->nbFaces.at(E)->back()->connectCell(c1, E);
// same as above, uses lambda:
		for_each(c0->nbFaces.at(E)->begin(), c0->nbFaces.at(E)->end(),
				[&c1](Face *testNbFace) { testNbFace->connectCell(c1, E); });

		// Reassign c0's East face(s)
		// c0's East face(s) become(s) c1's East face(s)
		swap(*(c0->nbFaces.at(E)), *(c1->nbFaces.at(E)));

//		c1->nbFaces.at(E)->push_back(c0->nbFaces.at(E)->back());
//		c0->nbFaces.at(E)->pop_back();

//		c1->nbFaces.at(E)->assign(c0->nbFaces.at(E)->begin(), c0->nbFaces.at(E)->end());


		// New face to lie between West and East newly-refined cells c0 and c1:
		// Factory method connects the Face to its neighbouring cells
		Face *fc = Face::createIntermediateFace(V, c0, c1);
		fc->set_length(c0->dy);
		fc->set_x(0.5*(c0->x+c1->x));
		fc->set_y(c0->y);

		c0->nbFaces.at(E)->clear();
		c1->nbFaces.at(W)->clear();

		// Now connect the Cells to the intermediate Face
		c0->nbFaces.at(E)->push_back(fc);
		c1->nbFaces.at(W)->push_back(fc);

		newFaces->push_back(fc);

		break;  // end case H:
	}

	case V:
	{

		// Reset refinement flags
		c0->refFlags.reset(doRefineY);
		c1->refFlags.reset(doRefineY);

		c0->refFlags.reset(doCoarsenY);
		c1->refFlags.reset(doCoarsenY);

		c0->refFlags.set(wasYRefined);
		c1->refFlags.set(wasYRefined);

		// Alter cell geometry
		c0->dy *= 0.5;
		c1->dy *= 0.5;
		c0->y -= 0.5*c0->dy;
		c1->y += 0.5*c1->dy;
		c0->lj++;
		c1->lj++;
		c0->j_idx *= 2;
		c1->j_idx = c0->j_idx+1;

		// Look at North and South faces:
		//  If cell has only one face to the north/south, split that face in two.
		// Create the new face to the right (back of the deque) of the initial face.
		for (auto& d: x_DIR) {   // d in {W,E}
			if (c0->nbFaces.at(d)->size() == 1) {

				Face *f = Face::createRefinedFace(orient, c0->nbFaces.at(d)->front());
				c1->nbFaces.at(d)->push_back(f);
				newFaces->push_back(f);
				f->connectCell(c1, d);

				// If there exists a cell to the west/east, connect that cell to the newly created faces.
				Cell *testNbCell = NULL;
				if (c0->nbFaces.at(d)->front()->get_nbCell(testNbCell, d)) {
					// Connect Face to Cell
					f->connectCell(testNbCell, oppositeDir(d));
					// Connect Cell to Face
					testNbCell->nbFaces.at(oppositeDir(d))->push_back(f);
				}
			}
			//  If cell already has two faces to the west/east, reassign [pos]most face to c1.
			else {
				// Connect Face to Cell c1
				c0->nbFaces.at(d)->back()->connectCell(c1, d);
				// Connect Cell c1 to Face
				c1->nbFaces.at(d)->push_back(c0->nbFaces.at(d)->back());
				// Disconnect Cell c0 from Face
				c0->nbFaces.at(d)->pop_back();
				// no need to adjust the 2nd face's N neighbour, it should already be mapped.
			}
		}

		// Connect c0's North face(s) to c1
//		c0->nbFaces.at(N)->front()->connectCell(c1, N);
//		if (c0->nbFaces.at(N)->size() > 1)
//			c0->nbFaces.at(N)->back()->connectCell(c1, N);

// Alt approach, iterate across faceDeque and apply connectCell using lambda fn:
		for_each(c0->nbFaces.at(N)->begin(), c0->nbFaces.at(N)->end(),
				[&c1](Face *testNbFace) { testNbFace->connectCell(c1, N); });
		// lambda: ^  capture 'c1' from enclosing scope [& : by reference, = : by value, etc.]

		// Reassign c0's North face(s) reference(s)
		// c0's North face(s) become(s) c1's North face(s)
		swap(*(c0->nbFaces.at(N)), *(c1->nbFaces.at(N)));
// Alt1
//		c1->nbFaces.at(N)->push_back(c0->nbFaces.at(N)->back());
//		c0->nbFaces.at(N)->pop_back();
// Alt2
//		c1->nbFaces.at(N)->assign(c0->nbFaces.at(N)->begin(), c0->nbFaces.at(N)->end());

		// New face to lie between West and East newly-refined cells c0 and c1:
		// Factory method connects the Face to its neighbouring cells
		Face *fc = Face::createIntermediateFace(H, c0, c1);
		fc->set_length(c0->dx);
		fc->set_x(c0->x);
		fc->set_y(0.5*(c0->y+c1->y));

		c0->nbFaces.at(N)->clear();
		c1->nbFaces.at(S)->clear();

		// Now connect the Cells to the intermediate Face
		c0->nbFaces.at(N)->push_back(fc);
		c1->nbFaces.at(S)->push_back(fc);

		newFaces->push_back(fc);


		break;  // end case V:
	}

	} // end switch(orient)

//	cout << "Cell::createRefinedCell: " << c0->id << ", " << c1->id << ": " << newFaces->size() << " faces" << endl;

	return(make_pair(newFaces, c1));

}
