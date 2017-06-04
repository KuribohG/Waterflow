#pragma once

#include "shared.h"
#include "simulation_cubic.h"
#include "marchingcubes.h"
#include "surface.h"

#ifdef OPENGL
#include "drawer.h"
#endif

class FluidSimulation {
private:
	int framenum = 0;
	int endframe = 0x7fffffff;
	char dumppref[100] = "objs/";
	SimulationCubic cubic;
	vector<WaterSource> sources;
	//vector<pair<int, Float> > lightPath[SHOW_SIZE_X][SHOW_SIZE_Y];
	aryf signed_dis;
	MarchingCubes meshcubes;
	Float dis[MAXGRID][MAXGRID][MAXGRID];
	MarkerParticle *nearest[MAXGRID][MAXGRID][MAXGRID];
public:
	FluidSimulation();
	void init_Density_3d();
	//void Draw_Vel(int i);
	//void Draw_Pixels(void);

	void Draw_On_Screen(void);
	void Read_Scene_File(const char * filename);
	void Step_Time(void);
	void Calculate_Signed_Distance(void);
	void Calculate_Nearest_Particle(void);
	void Get_Full_Velocity(void);
};