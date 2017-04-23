#pragma once

#include "shared.hpp"

class MarkerParticle {
public:
	Float x, y, z;
	MarkerParticle(Float _x,Float _y,Float _z):x(_x),y(_y),z(_z){}
};

void Mark_Water_By(vector<MarkerParticle> &particles, int *mask);
void Place_Particles(vector<MarkerParticle> &particles, int *mask);
void Advect_Particles(vector<MarkerParticle> &particles,Float *vx, Float *vy, Float *vz,int *mask);