#pragma once

#include "shared.hpp"
#include "gridmath.h"

class MarkerParticle {
public:
	Float x, y, z;
	MarkerParticle(Float _x,Float _y,Float _z):x(_x),y(_y),z(_z){}
};

void Mark_Water_By(vector<MarkerParticle> &particles, aryi &mask);
void Place_Particles(vector<MarkerParticle> &particles, aryi &mask);
void Advect_Particles(vector<MarkerParticle> &particles, aryf &vx, aryf &vy, aryf &vz, aryi &mask);