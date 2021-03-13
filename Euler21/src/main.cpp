//============================================================================
// Name        : main.cpp
// Author      : wak
// Version     : 2021-02-03                - initial take
// Copyright   : Your copyright notice
// Description : Execute Euler21 Solver
//============================================================================

using namespace std;

// Check memory usage:
//
// pmap -x `ps x | grep Euler21 | grep meshtest | awk '{print $1}'` | grep total
//


#include "Solver.h"
#include "Mesh.h"
#include "EulerIO.h"
#include "EulerDisplay.h"



int main(int argc, char* argv[]) {

	PayloadVar::setProblemType(GAS_DYNAMICS);

	cout << "argv[1]: " << argv[1] << endl;

	EulerDisplay *display = new EulerDisplay();
	Mesh *mesh = new Mesh();
	mesh->set_display(display);
	CB_EulerDisplay_drawFn drawFn = &EulerDisplay::drawMesh;
	mesh->set_cb_drawFn(drawFn);
	EulerIO *io = new EulerIO();

	io->readJsonIntoMesh(argv[1], mesh);


	Solver *solver = new Solver(mesh);
	io->readPhysParamsIntoSolver(argv[2], solver);





//	mesh->printMesh();
//	mesh->init_cell_and_face_IDs();
//	mesh->printMesh();

	for (int i = 0 ; i < 10 ; i++) {
		mesh->doUniformRefine(4);  //

		mesh->doUniformCoarsen(4);     // coarsen cells whose reflvl = 4
		mesh->doUniformCoarsen(3);     // coarsen cells whose reflvl = 3
		mesh->doUniformCoarsen(2);     // coarsen cells whose reflvl = 2
		mesh->doUniformCoarsen(1);     // coarsen cells whose reflvl = 1
	}

	mesh->doUniformRefine(4);

	display->drawMesh(mesh->cellMap, mesh->faceMap);

	cout << "sizeof(Cell): " << sizeof(Cell) << endl;
	cout << "sizeof(Face): " << sizeof(Face) << endl;
	cout << "sizeof(Mesh): " << sizeof(Mesh) << endl;



	while(display->checkQuitEvent() == 0) {
		SDL_Delay(5);
	}


	cout << "SUCCESS" << endl;



	return 0;
}


