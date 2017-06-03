#pragma once

#include "shared.h"
#include "gridmath.h"

class WaterSource {
public:
	int x0, x1, y0, y1, z0, z1;
	Float gen_rate;//every second, generate how many water per cell
	Float init_vx, init_vy, init_vz;
	int pourend;
	WaterSource(int _x0, int _x1, int _y0, int _y1, int _z0, int _z1, Float _gen_rate, Float _init_vx, Float _init_vy, Float _init_vz, int _pourend);
};

class MarkerParticle {
public:
	Float x, y, z, vx, vy, vz;
	MarkerParticle(Float _x,Float _y,Float _z):x(_x),y(_y),z(_z){}
	MarkerParticle(Float _x, Float _y, Float _z, Float _vx, Float _vy, Float _vz) :x(_x), y(_y), z(_z), vx(_vx), vy(_vy), vz(_vz){}
};

void Mark_Water_By(vector<MarkerParticle> &particles, aryi &mask);
void Init_Particles(vector<MarkerParticle> &particles, aryi &mask);
void Add_Single_Particle(vector<MarkerParticle> &particles, aryi &mask, int i, int j, int k);
void Add_Single_Particle(vector<MarkerParticle> &particles, aryi &mask, int i, int j, int k, Float vx, Float vy, Float vz);
void Get_Particles_Velocity(vector<MarkerParticle> &particles, aryf &vx, aryf &vy, aryf &vz, aryi &mask);
void Advect_Particles(vector<MarkerParticle> &particles, aryf &vx, aryf &vy, aryf &vz, aryi &mask);