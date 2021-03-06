/*
 * EulerDisplay.cpp
 *
 *  Created on: Feb. 23, 2021
 *      Author: wakeats
 */

#include "EulerDisplay.h"
#include <iostream>
#include <SDL2/SDL.h>
#include <SDL2/SDL2_gfxPrimitives.h>

EulerDisplay::EulerDisplay() {

	if (SDL_Init(SDL_INIT_VIDEO) < 0 ) {
		printf( "SDL could not initialize! SDL_Error: %s\n", SDL_GetError() );
	}
	else
	{
		//Create window
		window = SDL_CreateWindow( "Euler21", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
				                    SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN );

		if( window == NULL ) {
			printf( "Window could not be created! SDL_Error: %s\n", SDL_GetError() );
		}
		else {
			//Get window surface
//			screenSurface = SDL_GetWindowSurface( window );

			renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
//			renderer = SDL_CreateSoftwareRenderer(window);

			if (renderer == NULL) {
				printf( "Renderer could not be created! SDL Error: %s\n", SDL_GetError() );
			}
		}
	}

}

EulerDisplay::~EulerDisplay() {

	SDL_DestroyRenderer(renderer);
	//Destroy window
	SDL_DestroyWindow( window );
	//Quit SDL subsystems
	SDL_Quit();


}


int EulerDisplay::drawMeshNull(CellMap &cm, FaceMap &fm) {
	return 0;
}


int EulerDisplay::drawMesh(CellMap &cm, FaceMap &fm) {

	SDL_SetRenderDrawColor( renderer, 0xFF, 0xFF, 0xFF, 0xFF );
	SDL_RenderClear(renderer);
	SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);

//	cout << "rendering cells" << endl;
	for (auto& c: cm){

		int x = X_SCALE*(c.second->get_x() - 0.5*c.second->get_dx() + X_SHIFT);
		int y = Y_SCALE*(-c.second->get_y() - 0.5*c.second->get_dy() + Y_SHIFT);
		int dx = c.second->get_dx()*X_SCALE;
		int dy = c.second->get_dy()*Y_SCALE;

		SDL_SetRenderDrawColor( renderer, 0xFF, 0x00, 0x00, 0x88 );

//		unsigned char rc = 1000*(fabs(2*c.second->get_U()->get_dU(X, CV_DENS)));
//		unsigned char gc = 1000*(fabs(2*c.second->get_U()->get_dU(Y, CV_DENS)));
//		unsigned char bc = 200*c.second->get_U(CV_DENS); //->get_dU(X, CV_DENS);

//		unsigned char rc = 0.002*c.second->get_PV(PV_P);
//		unsigned char gc = 1.002*c.second->get_PV(PV_U);
//		unsigned char bc = 1.002*c.second->get_PV(PV_V);


		unsigned char rc = 200*c.second->get_U()->get_PHI(AV_SCHLIEREN);
		unsigned char gc = 200*c.second->get_U()->get_PHI(AV_SCHLIEREN);
		unsigned char bc = 200*c.second->get_U()->get_PHI(AV_SCHLIEREN);




		unsigned char al = 0.90*(1-c.second->get_PV(PV_RHO));

		SDL_SetRenderDrawColor( renderer, rc, gc, bc, 0x88 );

		SDL_Rect fillRect = {x, y, (int)dx, (int)dy};
		SDL_RenderDrawRect( renderer, &fillRect );
		SDL_RenderFillRect( renderer, &fillRect );

//		cout << c.first << ": (" << x << " , " << y << " , " << dx << " , " << dy << ")" << endl;

	}

//	cout << "rendering face points" << endl;

	for (auto& f: fm){


		int x = X_SCALE*(f.second->get_x()+X_SHIFT);
		int y = Y_SCALE*(-f.second->get_y()+Y_SHIFT);

		SDL_SetRenderDrawColor( renderer, 0x11, 0x11, 0xFF, 0x10 );
		SDL_Rect fillRect = {x-1, y-1, 3, 3};

//		if (f.second->faceRefFlags.test(doRecycleFace)) {
//			SDL_SetRenderDrawColor( renderer, 0xFF, 0x11, 0xFF, 0xFF );
//			fillRect.h = 5;
//			fillRect.w = 5;
//			fillRect.x = x-3;
//			fillRect.y = y-3;
//
//		}

		SDL_RenderDrawRect( renderer, &fillRect );
//		SDL_RenderDrawPoint( renderer, x, y );

//		cout << f.first << ": (" << x << " , " << y << ")" << endl;

	}

//	while (getSDLEvent() != SDL_KEYDOWN) {

	{

//		if(checkQuitEvent() == 1) {
		if (getSDLEvent() == SDL_MOUSEBUTTONDOWN) {
			return(-1);
		}

		SDL_Delay(10 );
	}

	SDL_RenderPresent(renderer);
	return(0);

}
