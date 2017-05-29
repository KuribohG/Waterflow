#pragma once

#include "shared.hpp"
#include "simulation_cubic.h"
#include "drawer.h"
#include "marchingcubes.h"
#include "surface.h"

class WaterSource {
public:
	int x0, x1, y0, y1, z0, z1;
	Float gen_rate;//every second, generate how many water per cell
	Float init_vx, init_vy, init_vz;
	WaterSource(int _x0, int _x1, int _y0, int _y1, int _z0, int _z1, Float _gen_rate, Float _init_vx, Float _init_vy, Float _init_vz);
};

class FluidSimulation {
private:
	int framenum = 0;
	SimulationCubic cubic;
	vector<WaterSource> sources;
	//vector<pair<int, Float> > lightPath[SHOW_SIZE_X][SHOW_SIZE_Y];
	aryf signed_dis;
	MarchingCubes meshcubes;
	Float dis[GRIDX][GRIDY][GRIDZ];
	MarkerParticle *nearest[GRIDX][GRIDY][GRIDZ];
public:
	FluidSimulation();
	void init_Density_3d();
	//void Draw_Vel(int i);
	//void Draw_Pixels(void);

	void Draw_On_Screen(void);
	void Read_Scene_File(const char * filename);
	void Pour_Source(void);
	void Step_Time(void);
	void Calculate_Signed_Distance(void);
	void Calculate_Nearest_Particle(void);
	void Get_Full_Velocity(void);
};