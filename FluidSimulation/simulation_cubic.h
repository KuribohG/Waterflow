#pragma once

#include "shared.hpp"
#include "surface.h"
#include "gridmath.h"

class SimulationCubic {
public:
	Float viscosity=0.0;
	Float diff=0.0;
	aryf vx, vx0, vy, vy0, vz, vz0;
	aryf p, p0;
	aryi mask;
	/*Float *vx;
	Float *vy;
	Float *vz;
	Float *vx0;
	Float *vy0;
	Float *vz0;
	Float *s;
	//Float *density;

	int *mask;*/

	vector<MarkerParticle> particles;

	SimulationCubic(void);

	~SimulationCubic();

	void Apply_External_Forces(void);

	void Linear_Solve(int axis, aryf & x, aryf & x0, Float a, Float wtsum, int iter);

	void Calc_Divergence(aryf & vx, aryf & vy, aryf & vz, aryf & div);

	void Project(aryf & vx, aryf & vy, aryf & vz, aryf & p, aryf & div);

	void Runge_Kutta(int i, int j, int k, Float delta, int iter, aryf & vx, aryf & vy, aryf & vz, Float & x, Float & y, Float & z);

	void Advect(int axis, aryf & f, aryf & f0, aryf & vx, aryf & vy, aryf & vz);

	void Step_Time();

};