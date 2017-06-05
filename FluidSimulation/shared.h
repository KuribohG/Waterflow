#pragma once


#define OPENGL
#define OPENMP

#include <iostream>
#include <cstdio>
#include <tuple>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <cassert>
#include <vector>
#include <list>
#include <float.h>
#include <ctime>
#include <omp.h>


#ifdef OPENGL

#include <GL/glut.h>

#endif

using namespace std;
typedef char BYT;
typedef float Float;

#ifdef NDEBUG
#define LOGM(...) fprintf(stderr,  __VA_ARGS__)
//#define LOGM(...) void(0);
#else
#define LOGM(...)  fprintf(stderr,  __VA_ARGS__)
#endif

#ifdef _WIN32
#define finite _finite
#endif

//#define ARYTPL template<int n,int m,int w> // without typename T
//#define ARYTPLT template<typename T, int n, int m, int w> // with typename T
//#define ARYDEF Array3D<T,n,m,w>
enum AXES { _X, _Y, _Z };
enum GridMaterial { WATER, SOLID, AIR };

const double INF = 1e15;

const Float TSDF_EPS = 1e-6f;
constexpr Float  EPS = 1e-3f;

const Float DENSITY = 20;

const Float g = 9.8;

const int FPS = 20;

const Float TIME_DELTA = 1.0 / (Float)FPS;

const int TOTAL_SCREEN = 3;

const double TOLERANCE = 1e-4;
const int LINSOLVER_ITER = 100;

const int SHOW_SIZE_X = 512;
const int SHOW_SIZE_Y = 512;

const int MAXGRID = 256;

extern int GRIDX;
extern int GRIDY;
extern int GRIDZ;

template<typename T>
inline T Scale_Along(T x, AXES axis) {
	if (axis == _X) return (T)(x + 0.0)*GRIDX / MAXGRID;
	else if (axis == _Y) return (T)(x + 0.0)*GRIDY / MAXGRID;
	else if (axis == _Z) return (T)(x + 0.0)*GRIDZ / MAXGRID;
}

template<typename T>
inline T Clip(T x, T low, T high) {
	if (x < low) return low;
	if (x > high) return high;
	return x;
}

//#define assert(x) if (!(x)) { asm("int $3"); }
/*inline int ID(int x, int y, int z) {
    assert(0 <= x&&x < GRIDX);
    assert(0 <= y&&y < GRIDY);
    assert(0 <= z&&z < GRIDZ);
    return x*GRIDY*GRIDZ + y*GRIDZ + z;
}*/

inline int random(int x, int y) {
    return rand() % (y - x + 1) + x;
}

inline Float randomF() {
    return random(0, 9999) / (10000.0);
}

inline bool validX(int x) {
    return 0 <= x&&x < GRIDX;
}

inline bool validY(int y) {
    return 0 <= y&&y < GRIDY;
}

inline bool validZ(int z) {
    return 0 <= z&&z < GRIDZ;
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
