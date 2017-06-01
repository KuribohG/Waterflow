
#include "drawer.h"
#ifdef OPENGL
#include "fluid_simulation.h"
#include <iostream>

static void ScreenCoor_to_ClipCoor(Float &x, Float &y, SCREENID_T screenid) {
    x = (x + screenid*SHOW_SIZE_X) * 2 / (SHOW_SIZE_X*TOTAL_SCREEN) - 1.0f;
    y = y * 2 / (SHOW_SIZE_Y) - 1.0f;
}
static void GridCoor_to_ClipCoor(Float &x, Float &y, SCREENID_T screenid) {
    x *= (SHOW_SIZE_X / GRIDY);
    y *= (SHOW_SIZE_Y / GRIDZ);
    x = (x + screenid*SHOW_SIZE_X) * 2 / (SHOW_SIZE_X*TOTAL_SCREEN) - 1.0f;
    y = y * 2 / (SHOW_SIZE_Y)- 1.0f;
}

void Draw_Particle_2d(vector<MarkerParticle> &particles, SCREENID_T screenid) {
	LOGM("draw particles\n");
	glBegin(GL_POINTS);
	glColor3f(1.0, 1.0, 1.0);
	for (MarkerParticle &p : particles) {
		int ix = floor(p.x);
		Float y = p.y, z = p.z;
		if (ix == GRIDX / 2) {
			GridCoor_to_ClipCoor(y, z, screenid);
			glVertex2f(y, z);
		}
	}
	glEnd();
	glFlush();
}

void Draw_Mask_2d(aryi &mask, SCREENID_T screenid) {
	LOGM("draw density\n");
	assert(screenid < 5);
	glBegin(GL_POINTS);
	for (int is = 0; is < SHOW_SIZE_X; is++) {
		for (int js = 0; js < SHOW_SIZE_Y; js++) {
			int i = is*GRIDY / SHOW_SIZE_X;
			int j = js*GRIDZ / SHOW_SIZE_Y;
			Float x = is, y = js;
			ScreenCoor_to_ClipCoor(x, y, screenid);
			int d = mask(GRIDX / 2, i, j);
			if (d == WATER) glColor3f(1, 1, 1);
			else if (d == AIR) glColor3f(0, 0, 0);
			else if (d == SOLID) glColor3f(0, 0, 1);
			glVertex2f(x, y);
		}
	}
	glEnd();
	glFlush();
}

void Draw_Density_2d(aryf &density, SCREENID_T screenid) {
    LOGM("draw density\n");
    assert(screenid < 5);
    glBegin(GL_POINTS);
	Float mx = -1;
	for (int j = 0; j < GRIDY; j++) {
		for (int k = 0; k < GRIDZ; k++) {
			mx = max(mx, fabs(density.get(GRIDX / 2, j, k)));
		}
	}
	mx = max(mx, (Float)1.0);
    for (int is = 0; is < SHOW_SIZE_X; is++) {
        for (int js = 0; js < SHOW_SIZE_Y; js++) {
            int i = is*GRIDY / SHOW_SIZE_X;
            int j = js*GRIDZ / SHOW_SIZE_Y;
            Float x = is, y = js;
            ScreenCoor_to_ClipCoor(x, y, screenid);
			Float d = fabs(density(GRIDX / 2, i, j) / mx);
			//if (d != 0) printf("%f ", d);
            glColor3f(d, d, d);
            glVertex2f(x, y);
        }
    }
    glEnd();
    glFlush();
}

void Draw_Velocity_2d(const aryf &vxs, const aryf &vys,const aryf &vzs, const aryi &mask, SCREENID_T screenid) {
    LOGM("draw velocity\n");
	int i = GRIDX / 2;
	for (int j = 0; j < GRIDY; j += GRIDY / 64) {
		for (int k = 0; k < GRIDZ; k += GRIDZ / 64) {
			int j0 = j, k0 = k;
			Float vy = Interpolation_Water_Velocity(_Y, vys, i + 0.5, j0 + 0.5, k0 + 0.5, mask, false)*TIME_DELTA;
			Float vz = Interpolation_Water_Velocity(_Z, vzs, i + 0.5, j0 + 0.5, k0 + 0.5, mask, false)*TIME_DELTA;
			vy *= 1, vz *= 1;
			//if (j == 29 && vz != 0) printf("get velocity: %f %f %f %f\n", j0 + 0.5, k0 + 0.5, vy, vz);
			Float y0 = j0 + 0.5, z0 = k0 + 0.5;
			Float y1 = y0 + vy, z1 = z0 + vz;
			//printf("%f %f %f %f\n", y0, z0, y1, z1);
			GridCoor_to_ClipCoor(y0, z0, screenid);
			GridCoor_to_ClipCoor(y1, z1, screenid);
			glColor3f(1.0, 1.0, 1.0);
			glBegin(GL_LINES);
			glVertex2f(y0, z0);
			glVertex2f(y1, z1);
			glEnd();
			glColor3f(0, 1, 0);
			glBegin(GL_POINTS);
			glVertex2f(y0, z0);
			glEnd();
		}
	}
    glFlush();
}

void Draw_Nearest(int i, MarkerParticle * nearest[GRIDX][GRIDY][GRIDZ]){
	for (int j = 0; j < GRIDY; j++) {
		for (int k = 0; k < GRIDZ; k++) {
			if (!nearest[i][j][k]) continue;
			MarkerParticle p = *nearest[i][j][k];
			//cout << p.x << " " << p.y << endl;
			Float y0 = j, z0 = k;
			Float y1 = p.y, z1 = p.z;
			GridCoor_to_ClipCoor(y0, z0, THIRD_SCREEN);
			GridCoor_to_ClipCoor(y1, z1, THIRD_SCREEN);
			glColor3f(1.0, 1.0, 1.0);
			glBegin(GL_LINES);
			glVertex2f(y0, z0);
			glVertex2f(y1, z1);
			glEnd();
		}
	}
	glFlush();
}

#endif