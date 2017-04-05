#pragma once

#include "shared.h"
#include "simulation_cubic.h"
#include "drawer.h"

class FluidSimulation {
private:
	SimulationCubic cubic;
	vector<pair<int, Float> > lightPath[SHOW_SIZE_X][SHOW_SIZE_Y];

public:
	FluidSimulation();
	void init_Density_3d();
	void Draw_Vel(int i);
	void Draw_Pixels(void);
	void Draw_On_Screen(void);
	void Step_Time(void);
};