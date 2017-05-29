#pragma once

#include "shared.hpp"
#include "surface.h"
#include "gridmath.h"
#include "pressuresolver.h"

class SimulationCubic {
private:
	aryf temp_dis;
public:
	Float viscosity=0.0;
	Float diff=0.0;
	aryf vx, vx0, vy, vy0, vz, vz0;
	aryf p, p0;
	aryi mask;
	PressureSolver solver;
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

	void Mark_Single_Water(int i, int j, int k);

	void Mark_Water(int x0, int x1, int y0, int y1, int z0, int z1);

	void Bound_Solid(void);

	void Apply_External_Forces(void);

	void Cancel_Air_Velocity(void);

	//void Linear_Solve(int axis, aryf & x, aryf & x0, Float a, Float wtsum, int iter);

	void Calc_Divergence(aryf & vx, aryf & vy, aryf & vz, aryf & div);

	void Project(aryf & vx, aryf & vy, aryf & vz, aryf & p, aryf & div);

	void Runge_Kutta(int axis, int i, int j, int k, Float delta, int iter, const aryf & vx, const aryf & vy, const aryf & vz, Float & x, Float & y, Float & z);

	void Advect_Velocity(int axis, aryf &f, const aryf &f0, const aryf &vx, const aryf &vy, const aryf &vz);

	void Step_Time();

};