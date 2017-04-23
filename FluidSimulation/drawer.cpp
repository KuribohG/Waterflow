#include "drawer.h"
#include "fluid_simulation.h"


static void ScreenCoor_to_ClipCoor(Float &x, Float &y, SCREENID_T screenid) {
    //x *= (SHOW_SIZE_X / GRID_SIZE_X);
    //y *= (SHOW_SIZE_Y / GRID_SIZE_Y);
    x = (x + screenid*SHOW_SIZE_X) * 2 / (SHOW_SIZE_X*TOTAL_SCREEN) - 1.0f;
    y = y * 2 / (SHOW_SIZE_Y) - 1.0f;
}
static void GridCoor_to_ClipCoor(Float &x, Float &y, SCREENID_T screenid) {
    x *= (SHOW_SIZE_X / GRID_SIZE_Y);
    y *= (SHOW_SIZE_Y / GRID_SIZE_Z);
    x = (x + screenid*SHOW_SIZE_X) * 2 / (SHOW_SIZE_X*TOTAL_SCREEN) - 1.0f;
    y = y * 2 / (SHOW_SIZE_Y)- 1.0f;
}

void Draw_Particle_2d(vector<MarkerParticle> &particles, SCREENID_T screenid) {
	LOGM("draw particles\n");
	glBegin(GL_POINTS);
	glColor3f(1.0, 1.0, 1.0);
	for (MarkerParticle &p : particles) {
		int ix = round(p.x);
		Float y = p.y, z = p.z;
		if (ix == 2) {
			GridCoor_to_ClipCoor(y, z, screenid);
			glVertex2f(y, z);
		}
	}
	glEnd();
	glFlush();
}

void Draw_Density_3d(Float *density, SCREENID_T screen,
                     const vector<pair<int, Float> > (*lightPath)[SHOW_SIZE_Y])
{
    glBegin(GL_POINTS);
    for (int i = 0; i < SHOW_SIZE_X; ++i)
        for (int j = 0; j < SHOW_SIZE_Y; ++j)
        {
            Float rgb(0), a(0);
            for (auto x : lightPath[i][j])
            {
                if (a > 0.995) break;
                Float rd = density[x.first] * x.second;
                rgb += rd * rd * (1 - a);
                a += (1 - a) * rd;
            }
            glColor3f(rgb, rgb, rgb);
            Float x = i, y = j;
            ScreenCoor_to_ClipCoor(x, y, screen);
            glVertex2f(x, y);
        }
    glEnd();
    glFlush();
}

void Draw_Mask_2d(int *mask, SCREENID_T screenid) {
	LOGM("draw density\n");
	assert(screenid < 5);
	glBegin(GL_POINTS);
	for (int is = 0; is < SHOW_SIZE_X; is++) {
		for (int js = 0; js < SHOW_SIZE_Y; js++) {
			int i = is*GRID_SIZE_Y / SHOW_SIZE_X;
			int j = js*GRID_SIZE_Z / SHOW_SIZE_Y;
			//LOGM("%d %d\n", i, j);
			//for (int k = 0; k < 3; k++) {
			//assert(0 <= density[ID(2, i, j)] && density[ID(2, i, j)] <= 1.0);
			Float x = is, y = js;
			ScreenCoor_to_ClipCoor(x, y, screenid);
			int d = mask[ID(2, i, j)];
			if (d == WATER) glColor3f(1, 1, 1);
			else if (d == AIR) glColor3f(0, 0, 0);
			else if (d == SOLID) glColor3f(0, 0, 1);
			//Float d = density[ID(2, i, j)];
			//glColor3f(d, d, d);
			glVertex2f(x, y);
			//pixels[p++] = BYT(density[ID(2, i, j)] * 255);
			//}
			//pixels[p++] = (BYT)255;
		}
	}
	glEnd();
	glFlush();
	//LOGM("pixel: %d\n", pixels[0]);
	//glDrawPixels(SHOW_SIZE_X, SHOW_SIZE_Y, GL_RGBA, GL_UNSIGNED_BYTE, pixels);
}

void Draw_Density_2d(Float *density, SCREENID_T screenid) {
    LOGM("draw density\n");
    assert(screenid < 5);
    glBegin(GL_POINTS);
    for (int is = 0; is < SHOW_SIZE_X; is++) {
        for (int js = 0; js < SHOW_SIZE_Y; js++) {
            int i = is*GRID_SIZE_Y / SHOW_SIZE_X;
            int j = js*GRID_SIZE_Z / SHOW_SIZE_Y;
            //LOGM("%d %d\n", i, j);
            //for (int k = 0; k < 3; k++) {
            assert(0 <= density[ID(2, i, j)] && density[ID(2, i, j)] <= 1.0);
            Float x = is, y = js;
            ScreenCoor_to_ClipCoor(x, y, screenid);
			Float d = abs(density[ID(2, i, j)]);
            glColor3f(d, d, d);
            glVertex2f(x, y);
            //pixels[p++] = BYT(density[ID(2, i, j)] * 255);
            //}
            //pixels[p++] = (BYT)255;
        }
    }
    glEnd();
    glFlush();
    //LOGM("pixel: %d\n", pixels[0]);
    //glDrawPixels(SHOW_SIZE_X, SHOW_SIZE_Y, GL_RGBA, GL_UNSIGNED_BYTE, pixels);
}

void Draw_Velocity_2d(int i, Float*, Float *vys, Float *vzs, SCREENID_T screenid) {
    LOGM("draw velocity\n");
    for (int j = 0; j < GRID_SIZE_Y; j ++) {
        for (int k = 0; k < GRID_SIZE_Z; k ++) {
            int j0 = j, k0 = k;
			Float vy = vys[ID(i, j0, k0)] * TIME_DELTA , vz = vzs[ID(i, j0, k0)] * TIME_DELTA;
            //LOGM("%f %f\n", vy, vz);
			Float y0 = j0 + 0.5, z0 = k0 + 0.5;
            Float y1 = y0 + vy, z1 = z0 + vz;
            //LOGM("%f %f %f %f\n", y0, z0, y1, z1);
            GridCoor_to_ClipCoor(y0, z0, screenid);
            GridCoor_to_ClipCoor(y1, z1, screenid);
            //y0 = y0 / SHOW_SIZE_X, y1 = y1 / SHOW_SIZE_X;
            //z0 = z0 / SHOW_SIZE_Y * 2 - 1, z1 = z1 / SHOW_SIZE_Y * 2 - 1;
            //y0 += SHOW_SIZE_X, y1 += SHOW_SIZE_X;
            glColor3f(1.0, 1.0, 1.0);
            glBegin(GL_LINES);
            //LOGM("lines\n");
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

