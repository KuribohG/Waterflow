#pragma once
#include "shared.h"

enum SCREENID_T { LEFT_SCREEN, RIGHT_SCREEN };
void Draw_Density_2d(Float *density, SCREENID_T screenid);
void Draw_Velocity_2d(int i, Float *vxs, Float *vys, Float *vzs, SCREENID_T screenid);