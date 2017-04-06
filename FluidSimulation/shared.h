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
#define LOGM(...) void(0)
#else
#define LOGM(...)  fprintf(stderr,  __VA_ARGS__)
#endif

constexpr Float  EPS = 1e-3f;


const int FPS = 100;

const Float TIME_DELTA = 1.0 / (Float)FPS;

const int TOTAL_SCREEN = 2;

const int LINSOLVER_ITER = 4;

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