/*
 * EulerDisplay.cpp
 *
 *  Created on: Feb. 23, 2021
 *      Author: wakeats
 */

#include "EulerDisplay.h"
#include <cerrno>
#include <csignal>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <string>
#include <sys/types.h>
#include <sys/wait.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL2_gfxPrimitives.h>
#include <unistd.h>
#include <vector>

namespace {

void reportExecError(int fd, int error) {
	const ssize_t bytesWritten = write(fd, &error, sizeof(error));
	(void)bytesWritten;
}

}

class EulerDisplay::VideoRecorder {
public:
	~VideoRecorder() {
		stop();
	}

	bool start(int width, int height, const std::string& outputPath) {
		int inputPipe[2] = {-1, -1};
		int execPipe[2] = {-1, -1};
		if (pipe(inputPipe) != 0 || pipe(execPipe) != 0) {
			std::cerr << "Could not create FFmpeg pipes: " << std::strerror(errno) << std::endl;
			closePipes(inputPipe, execPipe);
			return false;
		}

		if (fcntl(execPipe[1], F_SETFD, FD_CLOEXEC) == -1) {
			std::cerr << "Could not configure FFmpeg startup pipe: "
			          << std::strerror(errno) << std::endl;
			closePipes(inputPipe, execPipe);
			return false;
		}

		struct sigaction ignoreSigpipe {};
		ignoreSigpipe.sa_handler = SIG_IGN;
		sigemptyset(&ignoreSigpipe.sa_mask);
		if (sigaction(SIGPIPE, &ignoreSigpipe, &previousSigpipeAction) != 0) {
			std::cerr << "Could not configure FFmpeg pipe handling: "
			          << std::strerror(errno) << std::endl;
			closePipes(inputPipe, execPipe);
			return false;
		}
		sigpipeIgnored = true;

		pid_t child = fork();
		if (child == -1) {
			std::cerr << "Could not start FFmpeg: " << std::strerror(errno) << std::endl;
			closePipes(inputPipe, execPipe);
			restoreSigpipe();
			return false;
		}

		if (child == 0) {
			close(inputPipe[1]);
			close(execPipe[0]);
			if (dup2(inputPipe[0], STDIN_FILENO) == -1) {
				const int childError = errno;
				reportExecError(execPipe[1], childError);
				_exit(127);
			}
			close(inputPipe[0]);

			struct sigaction defaultSigpipe {};
			defaultSigpipe.sa_handler = SIG_DFL;
			sigemptyset(&defaultSigpipe.sa_mask);
			sigaction(SIGPIPE, &defaultSigpipe, nullptr);

			const std::string videoSize = std::to_string(width) + "x" + std::to_string(height);
			char* args[] = {
				const_cast<char*>("ffmpeg"),
				const_cast<char*>("-hide_banner"),
				const_cast<char*>("-loglevel"),
				const_cast<char*>("error"),
				const_cast<char*>("-y"),
				const_cast<char*>("-f"),
				const_cast<char*>("rawvideo"),
				const_cast<char*>("-pixel_format"),
				const_cast<char*>("rgb24"),
				const_cast<char*>("-video_size"),
				const_cast<char*>(videoSize.c_str()),
				const_cast<char*>("-framerate"),
				const_cast<char*>("20"),
				const_cast<char*>("-i"),
				const_cast<char*>("pipe:0"),
				const_cast<char*>("-an"),
				const_cast<char*>("-c:v"),
				const_cast<char*>("libx264"),
				const_cast<char*>("-preset"),
				const_cast<char*>("veryfast"),
				const_cast<char*>("-crf"),
				const_cast<char*>("20"),
				const_cast<char*>("-pix_fmt"),
				const_cast<char*>("yuv420p"),
				const_cast<char*>("-movflags"),
				const_cast<char*>("+faststart"),
				const_cast<char*>(outputPath.c_str()),
				nullptr
			};
			execvp(args[0], args);

			const int childError = errno;
			reportExecError(execPipe[1], childError);
			_exit(127);
		}

		close(inputPipe[0]);
		close(execPipe[1]);

		int childError = 0;
		ssize_t bytesRead;
		do {
			bytesRead = read(execPipe[0], &childError, sizeof(childError));
		} while (bytesRead == -1 && errno == EINTR);
		close(execPipe[0]);

		if (bytesRead != 0) {
			if (bytesRead == -1) {
				childError = errno;
			}
			close(inputPipe[1]);
			int status = 0;
			waitpid(child, &status, 0);
			std::cerr << "Could not execute FFmpeg: " << std::strerror(childError) << std::endl;
			restoreSigpipe();
			return false;
		}

		inputFd = inputPipe[1];
		pid = child;
		frameSize = static_cast<std::size_t>(width) * height * 3;
		pixels.resize(frameSize);
		std::cout << "Recording animation to " << outputPath << " at 10 fps." << std::endl;
		return true;
	}

	void capture(SDL_Renderer* renderer, int width, int height) {
		if (inputFd == -1) {
			return;
		}

		if (SDL_RenderReadPixels(renderer, nullptr, SDL_PIXELFORMAT_RGB24,
		                         pixels.data(), width * 3) != 0) {
			std::cerr << "Could not capture animation frame: " << SDL_GetError() << std::endl;
			stop();
			return;
		}

		std::size_t offset = 0;
		while (offset < pixels.size()) {
			const ssize_t bytesWritten = write(inputFd, pixels.data() + offset, pixels.size() - offset);
			if (bytesWritten == -1 && errno == EINTR) {
				continue;
			}
			if (bytesWritten <= 0) {
				std::cerr << "Could not send animation frame to FFmpeg: "
				          << std::strerror(errno) << std::endl;
				stop();
				return;
			}
			offset += static_cast<std::size_t>(bytesWritten);
		}
	}

private:
	static void closePipes(const int inputPipe[2], const int execPipe[2]) {
		for (int fd : {inputPipe[0], inputPipe[1], execPipe[0], execPipe[1]}) {
			if (fd != -1) {
				close(fd);
			}
		}
	}

	void restoreSigpipe() {
		if (sigpipeIgnored) {
			sigaction(SIGPIPE, &previousSigpipeAction, nullptr);
			sigpipeIgnored = false;
		}
	}

	void stop() {
		if (inputFd != -1) {
			close(inputFd);
			inputFd = -1;
		}

		if (pid != -1) {
			int status = 0;
			pid_t result;
			do {
				result = waitpid(pid, &status, 0);
			} while (result == -1 && errno == EINTR);

			if (result == -1) {
				std::cerr << "Could not wait for FFmpeg: " << std::strerror(errno) << std::endl;
			} else if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
				std::cerr << "FFmpeg failed while writing the animation." << std::endl;
			}
			pid = -1;
		}

		restoreSigpipe();
	}

	int inputFd = -1;
	pid_t pid = -1;
	std::size_t frameSize = 0;
	std::vector<unsigned char> pixels;
	struct sigaction previousSigpipeAction {};
	bool sigpipeIgnored = false;
};

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
			else {
				const char* configuredPath = std::getenv("EULER21_VIDEO");
				const std::string outputPath =
					configuredPath != nullptr && configuredPath[0] != '\0'
						? configuredPath
						: "Euler21.mp4";
				videoRecorder = std::make_unique<VideoRecorder>();
				videoRecorder->start(SCREEN_WIDTH, SCREEN_HEIGHT, outputPath);
			}
		}
	}

}

EulerDisplay::~EulerDisplay() {

	videoRecorder.reset();
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

	if (videoRecorder) {
		videoRecorder->capture(renderer, SCREEN_WIDTH, SCREEN_HEIGHT);
	}
	SDL_RenderPresent(renderer);
	return(0);

}
