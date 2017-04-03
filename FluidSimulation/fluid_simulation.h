#pragma once

#include "shared.h"
#include "simulation_cubic.h"
#include "drawer.h"

class FluidSimulation {
private:
	SimulationCubic cubic;
	
public:
	void Draw_Vel(int i);
	void Draw_Pixels(void);
	void Draw_On_Screen(void);
	void Step_Time(void);
};