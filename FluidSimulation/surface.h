#pragma once

#include "shared.h"
#include "gridmath.h"

class WaterSource {
public:
	int x0, x1, y0, y1, z0, z1;
	Float init_vx, init_vy, init_vz;
	int pourend;
	WaterSource(int _x0, int _x1, int _y0, int _y1, int _z0, int _z1, Float _init_vx, Float _init_vy, Float _init_vz, int _pourend);
};

class PeriodBox {
public:
	int x0, x1, y0, y1, z0, z1;
	Float semi_period;
	Float vx0, vy0, vz0, vx1, vy1, vz1;
	int end;
	PeriodBox(int _x0, int _x1, int _y0, int _y1, int _z0, int _z1, int _semi_period, Float _vx0, Float _vy0, Float _vz0, Float _vx1, Float _vy1, Float _vz1, int _end);
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
void Get_Particles_Velocity_FLIP(vector<MarkerParticle> &particles, aryf &vx, aryf &vy, aryf &vz, aryf &vx0, aryf &vy0, aryf &vz0, aryi &mask, Float flip);
void Advect_Particles(vector<MarkerParticle> &particles, aryf &vx, aryf &vy, aryf &vz, aryi &mask);