//============================================================================
// Name        : main.cpp
// Author      : wak
// Version     : 2021-02-03                - initial take
// Copyright   : Your copyright notice
// Description : Execute Euler21 Solver
//============================================================================

using namespace std;





#include "Mesh.h"
#include "EulerIO.h"
#include "EulerDisplay.h"



int main(int argc, char* argv[]) {


	cout << "argv[1]: " << argv[1] << endl;

	EulerDisplay *display = new EulerDisplay();

	Mesh *mesh = new Mesh();

	mesh->set_display(display);

	CB_EulerDisplay_drawFn drawFn = &EulerDisplay::drawMesh;

	mesh->set_cb_drawFn(drawFn);
	EulerIO *io = new EulerIO();

	io->readJsonIntoMesh(argv[1], mesh);

	mesh->printMesh();
//	mesh->init_cell_and_face_IDs();
//	mesh->printMesh();

for (int i = 0 ; i < 6 ; i++) {
	mesh->doUniformRefine(4);  //

	mesh->doUniformCoarsen(4);     // coarsen cells whose reflvl = 4
	mesh->doUniformCoarsen(3);     // coarsen cells whose reflvl = 3
	mesh->doUniformCoarsen(2);     // coarsen cells whose reflvl = 2
	mesh->doUniformCoarsen(1);     // coarsen cells whose reflvl = 1
}


	//	mesh->printMesh();

	mesh->doUniformRefine(4);  //

//	mesh->doUniformCoarsen(4);     // coarsen cells whose reflvl = 4
//	mesh->doUniformCoarsen(3);     // coarsen cells whose reflvl = 3
//	mesh->doUniformCoarsen(2);     // coarsen cells whose reflvl = 2
//	mesh->doUniformCoarsen(1);     // coarsen cells whose reflvl = 1

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

//
//
//	Json::Value root;
//	fstream ifs;
//	ifs.open(argv[1]);
//
//
//	Json::CharReaderBuilder builder;
//	builder["collectComments"] = true;
//	JSONCPP_STRING errs;
//	if (!parseFromStream(builder, ifs, &root, &errs)) {
//	    std::cout << errs << std::endl;
//	    return EXIT_FAILURE;
//	}
////	std::cout << root << std::endl;
//
//	Json::Value jsFaces = root["faces"];
//	Json::Value jsCells = root["cells"];
//
//	for (unsigned int i = 0 ; i < jsFaces.size() ; i++)
//	{
//		cout << jsFaces[i] << endl;
//	}
//
//	for (unsigned int i = 0 ; i < jsCells.size() ; i++)
//	{
//		cout << jsCells[i]["cid"] << endl;
//		cout << jsCells[i]["i"] << endl;
//		cout << jsCells[i]["j"] << endl;
//
//		Json::Value jsFaceIDs = jsCells[i]["faceIDs"];
//		Json::Value jsFaceN = jsFaceIDs["fN"];
//		for (unsigned int j = 0 ; j < jsFaceN.size() ; j++)
//		{
////			cout << jsFaceN[j] << endl;
//			cout << "fid: " << jsFaceN[j].get("fid", -1).asInt() << endl;
//		}
//	}
//
//
//
//	cout << "Done." << endl;
//
//
////	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!
//	return 0;
//}
