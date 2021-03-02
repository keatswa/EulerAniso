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

	Json::Value jsFaces = root["faces"];
	Json::Value jsCells = root["cells"];


	unsigned int nCells = jsCells.size();
	unsigned int nFaces = jsFaces.size();

	// Add Cells to mesh
	Cell tmpCell;

	for (unsigned int i = 0 ; i < nCells ; i++)
	{
		tmpCell.set_id( jsCells[i]["cid"].asUInt());
		tmpCell.set_i_idx(jsCells[i]["i"].asUInt());
		tmpCell.set_j_idx(jsCells[i]["j"].asUInt());
		tmpCell.set_li(0);
		tmpCell.set_lj(0);
		mesh->init_addCell(tmpCell);
	}

	// Add Faces to mesh
	Face tmpFace;
	for (unsigned int i = 0 ; i < nFaces ; i++)
	{
		tmpFace.set_id( jsFaces[i]["fid"].asUInt());
		tmpFace.set_length(jsFaces[i]["length"].asDouble());
		tmpFace.set_x( jsFaces[i]["x"].asDouble());
		tmpFace.set_y( jsFaces[i]["y"].asDouble());
		mesh->init_addFace(tmpFace);
	}


	// Connect Cells to Faces
	for (unsigned int i = 0 ; i < nCells ; i++)
	{
		Json::Value jsFaceIDs = jsCells[i]["faceIDs"];
		Json::Value jsFaceDir;
		jsFaceDir = jsFaceIDs["fN"];
		for (unsigned int j = 0 ; j < jsFaceDir.size() ; j++) {
			mesh->init_mapFaceToCell(jsCells[i]["cid"].asUInt(), jsFaceDir[j].get("fid", -1).asUInt(), N);
		}
		jsFaceDir = jsFaceIDs["fS"];
		for (unsigned int j = 0 ; j < jsFaceDir.size() ; j++) {
			mesh->init_mapFaceToCell(jsCells[i]["cid"].asUInt(), jsFaceDir[j].get("fid", -1).asUInt(), S);
		}
		jsFaceDir = jsFaceIDs["fW"];
		for (unsigned int j = 0 ; j < jsFaceDir.size() ; j++) {
			mesh->init_mapFaceToCell(jsCells[i]["cid"].asUInt(), jsFaceDir[j].get("fid", -1).asUInt(), W);
		}
		jsFaceDir = jsFaceIDs["fE"];
		for (unsigned int j = 0 ; j < jsFaceDir.size() ; j++) {
			mesh->init_mapFaceToCell(jsCells[i]["cid"].asUInt(), jsFaceDir[j].get("fid", -1).asUInt(), E);
//TEST			mesh->mapFaceToCell(jsCells[i]["cid"].asUInt(), jsFaceDir[j].get("fid", -1).asUInt(), E);
//TEST			mesh->mapFaceToCell(jsCells[i]["cid"].asUInt(), jsFaceDir[j].get("fid", -1).asUInt(), E);
		}


	}






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
//	//			cout << jsFaceN[j] << endl;
//			cout << "fid: " << jsFaceN[j].get("fid", -1).asInt() << endl;
//		}
//	}



	cout << "Done." << endl;





	return(0);

}





int EulerIO::writeMeshToJson(char *jsonFileName, const Mesh &mesh) {

	return(0);
}
