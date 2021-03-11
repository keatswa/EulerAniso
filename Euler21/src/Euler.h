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
#include <map>
#include <tsl/robin_map.h>
#include <tuple>
#include <functional>


class Face;
class Cell;
class EulerDisplay;
class GasDynFlux;
class GasDynVar;

const int nEqn = 4;


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


enum SIGN {
	NEG=0, POS=1
};

constexpr std::initializer_list<SIGN> all_SIGN = {NEG,POS};


enum RefinementFlags
{
	doRefineX,
	doRefineY,
	doCoarsenX,
	doCoarsenY,
	wasXRefined,
	wasYRefined,
	doRecycleCell,
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



enum BCType
{
	WALL,
	INLET,
	OUTLET,
	SYMMETRY,
	PERIODIC,
	FREESTREAM,
	NONE,
	UNKNOWN
};



typedef float cfdFloat;
typedef deque<Face*> FaceDeque;
typedef map<DIR, FaceDeque*> NeighbouringFaces;
typedef map<SIGN, Cell*> NeighbouringCells;

//typedef tsl::robin_map<unsigned long, Cell*> CellMap;
//typedef tsl::robin_map<unsigned long, Face*> FaceMap;
typedef unordered_map<unsigned long, Cell*> CellMap;
typedef unordered_map<unsigned long, Face*> FaceMap;

// The process of refining a Cell will return a new Cell and one or more new Face objects.
// RefinedCellFaceGroup is used as a return type.
typedef pair<FaceDeque*, Cell*> RefinedCellFaceGroup;


typedef int (EulerDisplay::*CB_EulerDisplay_drawFn) (CellMap& cm, FaceMap& fm);


// Cell payload typedefs
//  ConservedVariables - deque[# timesteps] of vector [# variables = 4 in 2D]

//typedef GasDynFlux Flux;
//typedef GasDynVar  ConsVar;


#endif /* EULER_H_ */
