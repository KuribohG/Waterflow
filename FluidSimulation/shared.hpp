#pragma once

#include <iostream>
#include <cstdio>
#include <tuple>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <cassert>
#include <vector>
#include <list>
#include <GL/glut.h>
using namespace std;
typedef char BYT;
typedef float Float;

#ifdef NDEBUG
#define LOGM(...) fprintf(stderr,  __VA_ARGS__)
#else
#define LOGM(...)  fprintf(stderr,  __VA_ARGS__)
#endif

enum GridMaterial { WATER, SOLID, AIR };

constexpr Float  EPS = 1e-3f;

const Float g = 9.8;

const int FPS = 10;

const Float TIME_DELTA = 1.0 / (Float)FPS;

const int TOTAL_SCREEN = 3;

const int LINSOLVER_ITER = 100;

const int SHOW_SIZE_X = 512;
const int SHOW_SIZE_Y = 512;

const int GRID_SIZE_X = 10;
const int GRID_SIZE_Y = 64;
const int GRID_SIZE_Z = 64;
//#define assert(x) if (!(x)) { asm("int $3"); }
inline int ID(int x, int y, int z) {
    assert(0 <= x&&x < GRID_SIZE_X);
    assert(0 <= y&&y < GRID_SIZE_Y);
    assert(0 <= z&&z < GRID_SIZE_Z);
    return x*GRID_SIZE_Y*GRID_SIZE_Z + y*GRID_SIZE_Z + z;
}

inline int random(int x, int y) {
    return rand() % (y - x + 1) + x;
}

inline Float randomF() {
    return random(0, 9999) / (10000.0);
}

inline bool validX(int x) {
    return 0 <= x&&x < GRID_SIZE_X;
}

inline bool validY(int y) {
    return 0 <= y&&y < GRID_SIZE_Y;
}

inline bool validZ(int z) {
    return 0 <= z&&z < GRID_SIZE_Z;
}

inline bool Valid_Water(int x, int y, int z, int *mask) {
	if (!validX(x)) return false;
	if (!validY(y)) return false;
	if (!validZ(y)) return false;
	return mask[ID(x, y, z)] == WATER;
}

inline Float Interpolation_In_Water_3D(Float *f, Float x, Float y, Float z, int *mask) {
	/*int x0 = floor(x-0.5), x1 = x0 + 1;
	int y0 = floor(y-0.5), y1 = y0 + 1;
	int z0 = floor(z-0.5), z1 = z0 + 1;
	Float s0 = x - x0 - 0.5, s1 = x1 - x + 0.5;
	Float u0 = y - y0 - 0.5, u1 = y1 - y + 0.5;
	Float v0 = z - z0 - 0.5, v1 = z1 - z + 0.5;
	Float v = 0, sw = 0;*/
	x -= 0.5, y -= 0.5, z -= 0.5;
	int x0 = floor(x), x1 = x0 + 1;
	int y0 = floor(y), y1 = y0 + 1;
	int z0 = floor(z), z1 = z0 + 1;
	Float s0 = x - x0, s1 = x1 - x;
	Float u0 = y - y0, u1 = y1 - y;
	Float v0 = z - z0, v1 = z1 - z;
	Float v = 0, sw = 0;
	if (Valid_Water(x0, y0, z0, mask)) v += s1*u1*v1*f[ID(x0, y0, z0)], sw += s1*u1*v1;
	if (Valid_Water(x0, y0, z0, mask)) v += s1*u1*v0*f[ID(x0, y0, z1)], sw += s1*u1*v0;
	if (Valid_Water(x0, y0, z0, mask)) v += s1*u0*v1*f[ID(x0, y1, z0)], sw += s1*u0*v1;
	if (Valid_Water(x0, y0, z0, mask)) v += s1*u0*v0*f[ID(x0, y1, z1)], sw += s1*u0*v0;
	if (Valid_Water(x0, y0, z0, mask)) v += s0*u1*v1*f[ID(x1, y0, z0)], sw += s0*u1*v1;
	if (Valid_Water(x0, y0, z0, mask)) v += s0*u1*v0*f[ID(x1, y0, z1)], sw += s0*u1*v0;
	if (Valid_Water(x0, y0, z0, mask)) v += s0*u0*v1*f[ID(x1, y1, z0)], sw += s0*u0*v1;
	if (Valid_Water(x0, y0, z0, mask)) v += s0*u0*v0*f[ID(x1, y1, z1)], sw += s0*u0*v0;
	if (sw == 0) return 0;
	return v / sw;
}

/*inline Float Interpolation_3D(Float *f, Float x, Float y, Float z) {
	int x0 = floor(x), x1 = x0 + 1;
	int y0 = floor(y), y1 = y0 + 1;
	int z0 = floor(z), z1 = z0 + 1;
	Float s0 = x - x0, s1 = x1 - x;
	Float u0 = y - y0, u1 = y1 - y;
	Float v0 = z - z0, v1 = z1 - z;
	Float v = 0, sw = 0;
	if (validX(x0) && validY(y0) && validZ(z0)) v += s1*u1*v1*f[ID(x0, y0, z0)], sw += s1*u1*v1;
	if (validX(x0) && validY(y0) && validZ(z1)) v += s1*u1*v0*f[ID(x0, y0, z1)], sw += s1*u1*v0;
	if (validX(x0) && validY(y1) && validZ(z0)) v += s1*u0*v1*f[ID(x0, y1, z0)], sw += s1*u0*v1;
	if (validX(x0) && validY(y1) && validZ(z1)) v += s1*u0*v0*f[ID(x0, y1, z1)], sw += s1*u0*v0;
	if (validX(x1) && validY(y0) && validZ(z0)) v += s0*u1*v1*f[ID(x1, y0, z0)], sw += s0*u1*v1;
	if (validX(x1) && validY(y0) && validZ(z1)) v += s0*u1*v0*f[ID(x1, y0, z1)], sw += s0*u1*v0;
	if (validX(x1) && validY(y1) && validZ(z0)) v += s0*u0*v1*f[ID(x1, y1, z0)], sw += s0*u0*v1;
	if (validX(x1) && validY(y1) && validZ(z1)) v += s0*u0*v0*f[ID(x1, y1, z1)], sw += s0*u0*v0;
	if (sw == 0) return 0;
	return v / sw;
}*/