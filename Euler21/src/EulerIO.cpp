/*
 * EulerIO.cpp
 *
 *  Created on: Feb. 5, 2021
 *      Author: wakeats
 */

#include "EulerIO.h"

EulerIO::EulerIO() {
	// TODO Auto-generated constructor stub

}

EulerIO::~EulerIO() {
	// TODO Auto-generated destructor stub
}




int EulerIO::readJsonIntoMesh(char *jsonFileName, Mesh *mesh) {

	Json::Value root;
	fstream ifs;
	ifs.open(jsonFileName);


	Json::CharReaderBuilder builder;
	builder["collectComments"] = true;
	JSONCPP_STRING errs;
	if (!parseFromStream(builder, ifs, &root, &errs)) {
		std::cout << errs << std::endl;
		return EXIT_FAILURE;
	}


	mesh->set_bounds(root["x0"].asDouble(), root["y0"].asDouble(),
			         root["x1"].asDouble(), root["y1"].asDouble());


	Json::Value jsFaces = root["faces"];
	Json::Value jsCells = root["cells"];

	unsigned int nCells = jsCells.size();
	unsigned int nFaces = jsFaces.size();

	// Add Cells to mesh
	Cell tmpCell;

	for (unsigned int i = 0 ; i < nCells ; i++)
	{
		tmpCell.set_id( jsCells[i]["cid"].asLargestUInt());   // OVERRIDDEN BELOW
		tmpCell.set_i_idx(jsCells[i]["i"].asUInt());
		tmpCell.set_j_idx(jsCells[i]["j"].asUInt());
		tmpCell.set_li(0);
		tmpCell.set_lj(0);
		cout << "cid: " << tmpCell.get_id() ;
		tmpCell.init_id();                                    // OVERRIDES ABOVE
		cout << " -> " << tmpCell.get_id() << endl;
		mesh->init_addCell(tmpCell);
	}
	// Add Faces to mesh
	Face tmpFace;
	for (unsigned int i = 0 ; i < nFaces ; i++)
	{
		tmpFace.set_id( jsFaces[i]["fid"].asLargestUInt());   // OVERRIDDEN BELOW

		tmpFace.set_length(jsFaces[i]["length"].asDouble());
		tmpFace.set_x( jsFaces[i]["x"].asDouble());
		tmpFace.set_y( jsFaces[i]["y"].asDouble());
		tmpFace.set_is_bc(true);
		string strBCType = jsFaces[i]["bcType"].asString();
		if (strBCType.compare("NA") == 0) {
			tmpFace.set_bcType(NONE);
			tmpFace.set_is_bc(false);
		}
		else if (strBCType.compare("WALL") == 0)
			tmpFace.set_bcType(WALL);
		else if (strBCType.compare("INLET") == 0)
			tmpFace.set_bcType(INLET);
		else if (strBCType.compare("OUTLET") == 0)
			tmpFace.set_bcType(OUTLET);
		else if (strBCType.compare("SYMMETRY") == 0)
			tmpFace.set_bcType(SYMMETRY);
		else if (strBCType.compare("PERIODIC") == 0)
			tmpFace.set_bcType(PERIODIC);
		else if (strBCType.compare("FREESTREAM") == 0)
			tmpFace.set_bcType(FREESTREAM);
		else
			tmpFace.set_bcType(UNKNOWN);



		tmpFace.set_id(mesh->provideNewFaceID(tmpFace.get_x(), tmpFace.get_y()));   // OVERRIDES ABOVE
		mesh->init_addFace(tmpFace);  // calls copy constructor

		cout << "fid " << tmpFace.get_id() << " -> " << mesh->provideNewFaceID(tmpFace.get_x(), tmpFace.get_y()) << endl;
	}


	// Connect Cells to Faces
	for (unsigned int i = 0 ; i < nCells ; i++)
	{
		Json::Value jsFaceIDs = jsCells[i]["faceIDs"];
		Json::Value jsFaceDir;
		jsFaceDir = jsFaceIDs["fN"];
		for (unsigned int j = 0 ; j < jsFaceDir.size() ; j++) {
			mesh->init_mapFaceToCell(jsCells[i]["cid"].asLargestUInt(), jsFaceDir[j].get("fid", -1).asLargestUInt(), N);
		}
		jsFaceDir = jsFaceIDs["fS"];
		for (unsigned int j = 0 ; j < jsFaceDir.size() ; j++) {
			mesh->init_mapFaceToCell(jsCells[i]["cid"].asLargestUInt(), jsFaceDir[j].get("fid", -1).asLargestUInt(), S);
		}
		jsFaceDir = jsFaceIDs["fW"];
		for (unsigned int j = 0 ; j < jsFaceDir.size() ; j++) {
			mesh->init_mapFaceToCell(jsCells[i]["cid"].asLargestUInt(), jsFaceDir[j].get("fid", -1).asLargestUInt(), W);
		}
		jsFaceDir = jsFaceIDs["fE"];
		for (unsigned int j = 0 ; j < jsFaceDir.size() ; j++) {
			mesh->init_mapFaceToCell(jsCells[i]["cid"].asLargestUInt(), jsFaceDir[j].get("fid", -1).asLargestUInt(), E);
		}


	}




	cout << "Done." << endl;





	return(0);

}



int EulerIO::writeMeshToJson(char *jsonFileName, const Mesh &mesh) {

	return(0);
}







int EulerIO::readPhysParamsIntoSolver(char *jsonFileName, Solver *solver) {


	Json::Value root;
	fstream ifs;
	ifs.open(jsonFileName);


	Json::CharReaderBuilder builder;
	builder["collectComments"] = true;
	JSONCPP_STRING errs;
	if (!parseFromStream(builder, ifs, &root, &errs)) {
		std::cout << errs << std::endl;
		return EXIT_FAILURE;
	}

	solver->setSolverParams(root["reflevelMax"].asInt(),
							root["reflevelMin"].asInt(),
							root["errorTol"].asDouble(),
							root["dt"].asDouble(),
							root["tEnd"].asDouble(),
							root["maxCFL"].asDouble(),
							root["maxIter"].asInt()        );




	string caseName = root["case"].asString();

	Json::Value jsInlets = root["inlets"];
	unsigned int nInlets = jsInlets.size();

	if (caseName.compare("BackwardStep") == 0) {

		// Assign inlets to faces first, before assigning initial conditions.
		// Assignment of initial conditions may require mesh to be refined, so
		// additional faces would not receive inlet conditions.

		for (unsigned int i = 0 ; i < nInlets ; i++)
		{

			solver->setInletCondition(jsInlets[i]["fid"].asLargestUInt(),
									  jsInlets[i]["rho"].asDouble(),
									  jsInlets[i]["u"].asDouble(),
									  jsInlets[i]["v"].asDouble(),
									  jsInlets[i]["p"].asDouble());

		}





		// Set up backward step params in Solver
		solver->setBackwardStepParams(root["rho0"].asDouble(),
									    root["u0"].asDouble(),
									    root["v0"].asDouble(),
									    root["p0"].asDouble(),
									  root["rho1"].asDouble(),
									    root["u1"].asDouble(),
									    root["v1"].asDouble(),
									    root["p1"].asDouble() );



		solver->assignBackwardStepICs();


	} else if (caseName.compare("ForwardStep") == 0) {




	} else {
		cout << "Unknown case. " << endl;
		return EXIT_FAILURE;
	}








	return 0;

}


















