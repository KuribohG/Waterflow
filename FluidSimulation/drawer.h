#pragma once
#include "shared.hpp"
#include "surface.h"

enum SCREENID_T { LEFT_SCREEN, RIGHT_SCREEN, THIRD_SCREEN };
void Draw_Density_2d(Float *density, SCREENID_T screenid);
void Draw_Particle_2d(vector<MarkerParticle> &particles, SCREENID_T screenid);
void Draw_Velocity_2d(int i, Float *vxs, Float *vys, Float *vzs, SCREENID_T screenid);
void Draw_Density_3d(Float *density, SCREENID_T screen,
                     const vector<pair<int, Float> > (*lightPath)[SHOW_SIZE_Y]);
void Draw_Mask_2d(int *mask, SCREENID_T screenid);
