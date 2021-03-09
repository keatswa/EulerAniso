/*
 * EulerDisplay.h
 *
 *  Created on: Feb. 23, 2021
 *      Author: wakeats
 */

#ifndef EULERDISPLAY_H_
#define EULERDISPLAY_H_

#include "Euler.h"
#include "Mesh.h"
#include <SDL2/SDL.h>
#include <SDL2/SDL2_gfxPrimitives.h>


class EulerDisplay {
public:
	EulerDisplay();
	virtual ~EulerDisplay();

	//Screen dimension constants
	const static int SCREEN_WIDTH = 800;
	const static int SCREEN_HEIGHT = 800;

	const static int X_SHIFT = 50;
	const static int Y_SHIFT = 50;

	const static int X_SCALE = 8;
	const static int Y_SCALE = 8;

//	int drawMesh(Mesh &m);

	int drawMesh(CellMap &cm, FaceMap &fm);
	int drawMeshNull(CellMap &cm, FaceMap &fm);

	int testFn() { cout << "hello" << endl; return(0); }

	int checkQuitEvent() {

		if (SDL_PollEvent(&e)) {
			if (e.type == SDL_QUIT || e.type == SDL_MOUSEBUTTONUP || e.type == SDL_MOUSEBUTTONDOWN) {
				return(1);
			}
		}
		return(0);
	}

private:


	SDL_Window   *window = NULL;
	SDL_Renderer *renderer = NULL;
	SDL_Event    e;


};

#endif /* EULERDISPLAY_H_ */
