# EulerAniso
2D anisotropic mesh refinement for the compressible Euler equations

Modernized implementation of the method in my 2004 thesis: http://hdl.handle.net/10012/948

Graphical output using SDL2

JSON input/output using jsoncpp

## Build and run

Install the SDL2, SDL2_gfx, and JsonCpp development libraries, then build from
the repository root:

```sh
make
```

The executable is written to `build/euler21`. Run it with the mesh and
physical-parameters JSON files:

```sh
./build/euler21 Euler21/meshtest.json Euler21/physparams.json
```
