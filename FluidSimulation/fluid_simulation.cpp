#include "fluid_simulation.h"



void FluidSimulation::Draw_On_Screen(void){
	//printf("draw_on_screen\n");
	Draw_Density_2d(cubic.density, RIGHT_SCREEN);
	Draw_Velocity_2d(2, cubic.vx, cubic.vy, cubic.vz, LEFT_SCREEN);
}

void FluidSimulation::Step_Time(void){
	cubic.Step_Time();
}

