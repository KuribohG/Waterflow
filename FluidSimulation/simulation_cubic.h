#pragma once

#include "shared.h"

enum GridMaterial { WATER, SOLID };

class SimulationCubic {
public:
	Float viscosity=0.0;
	Float diff=0.0;
	Float *vx;
	Float *vy;
	Float *vz;
	Float *vx0;
	Float *vy0;
	Float *vz0;
	Float *s;
	Float *density;

	int *mask;

	SimulationCubic(void);

	~SimulationCubic();

	void Linear_Solve(int axis, Float * x, Float * x0, Float a, Float wtsum, int iter);

	void Diffuse(int axis, Float * x, Float * x0, Float diff, Float dt, int iter);

	void Project(Float * vx, Float * vy, Float * vz, Float * pressure, Float * dv);

	void Advect(int axis, Float * density, Float * density0, Float * vx, Float * vy, Float * vz);

	void Step_Time();

};