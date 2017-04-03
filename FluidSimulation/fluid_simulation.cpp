#include "fluid_simulation.h"

void FluidSimulation::Draw_Vel(int i) {
	Draw_Velocity(2, cubic.vx, cubic.vy, cubic.vz, 1);
}

void FluidSimulation::Draw_Pixels(void) {
	Draw_Density(cubic.density);
}

void FluidSimulation::Draw_On_Screen(void){
	//printf("draw_on_screen\n");
	Draw_Pixels();
	Draw_Vel(2);
}

void FluidSimulation::Step_Time(void){
	cubic.Step_Time();
}

