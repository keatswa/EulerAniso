/*
 * EulerIO.h
 *
 *  Created on: Feb. 5, 2021
 *      Author: wakeats
 */

#ifndef EULERIO_H_
#define EULERIO_H_


#include "Mesh.h"
#include "json/json.h"
#include <fstream>
#include <iostream>


class EulerIO {
public:
	EulerIO();
	virtual ~EulerIO();


	int readJsonIntoMesh(char* jsonFileName, Mesh* mesh);

	int writeMeshToJson(char* jsonFileName, const Mesh& mesh);






};

#endif /* EULERIO_H_ */
