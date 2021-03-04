/*
 * Euler.h
 *
 *  Created on: Feb. 14, 2021
 *      Author: wakeats
 */
using namespace std;

#ifndef EULER_H_
#define EULER_H_

#include <deque>
#include <unordered_map>
#include <tuple>
#include <functional>


class Face;
class Cell;
class EulerDisplay;

enum ORIENTATION {
	H=0, V=1
};

enum DIR {
    N=0,
    S=1,
	W=2,
	E=3
};

constexpr std::initializer_list<DIR> all_DIR = {N,S,W,E};
constexpr std::initializer_list<DIR>   x_DIR = {    W,E};
constexpr std::initializer_list<DIR>   y_DIR = {S,N    };

// Defined in Cell.cpp:
extern DIR oppositeDir(DIR d);

//DIR oppositeDir(DIR d) {
//	DIR oppDir;
//
//	switch (d) {
//	case N:
//		oppDir = S;
//		break;
//	case E:
//		oppDir = W;
//		break;
//	case S:
//		oppDir = N;
//		break;
//	case W:
//		oppDir = E;
//		break;
//	}
//
//	return(oppDir);
//
//}

//for (auto& d : all_DIR) {
//    // Do job with d
//}

enum SIGN {
	NEG=0, POS=1
};

constexpr std::initializer_list<SIGN>   all_SIGN = {NEG,POS};


enum RefinementFlags
{
	doRefineX,
	doRefineY,
	doCoarsenX,
	doCoarsenY,
	wasXRefined,
	wasYRefined,
	doRecycleCell,
//	prevCyclXRefd,	// was x-refined in the previous Cycle
//	prevCyclYRefd,
	flaggedRefineX,
	flaggedRefineY,
	NUM_CELLREFINEMENT_FLAGS
};


enum FaceRefinementFlags
{
	doCoarsen,
	doRecycleFace,
	NUM_FACEREFINEMENT_FLAGS
};




typedef float cfdFloat;
typedef deque<Face*> FaceDeque;
typedef unordered_map<DIR, FaceDeque*> NeighbouringFaces;
typedef unordered_map<SIGN, Cell*> NeighbouringCells;

typedef unordered_map<unsigned int, Cell*> CellMap;
typedef unordered_map<unsigned int, Face*> FaceMap;

typedef pair<FaceDeque*, Cell*> RefinedCellFaceGroup;

typedef int (EulerDisplay::*CB_EulerDisplay_drawFn) (CellMap& cm, FaceMap& fm);

#endif /* EULER_H_ */
