#include "surface.h"



void Mark_Water_By(vector<MarkerParticle> &particles, int *mask) {
	for (int i = 0; i < GRID_SIZE_X; i++) {
		for (int j = 0; j < GRID_SIZE_Y; j++) {
			for (int k = 0; k < GRID_SIZE_Z; k++) {
				if (mask[ID(i, j, k)] != SOLID) mask[ID(i, j, k)] = AIR;
			}
		}
	}
	for (MarkerParticle &p : particles) {
		int x = round(p.x), y = round(p.y), z = round(p.z);
		if (mask[ID(x, y, z)] != SOLID) mask[ID(x, y, z)] = WATER;
	}
}

void Place_Particles(vector<MarkerParticle> &particles, int *mask) {
	particles.clear();
	for (int i = 0; i < GRID_SIZE_X; i++) {
		for (int j = 0; j < GRID_SIZE_Y; j++) {
			for (int k = 0; k < GRID_SIZE_Z; k++) {
				if (mask[ID(i, j, k)] == WATER) {
					particles.push_back(MarkerParticle(i - 0.25, j - 0.25, k - 0.25));
					particles.push_back(MarkerParticle(i - 0.25, j - 0.25, k + 0.25));
					particles.push_back(MarkerParticle(i - 0.25, j + 0.25, k - 0.25));
					particles.push_back(MarkerParticle(i - 0.25, j + 0.25, k + 0.25));
					particles.push_back(MarkerParticle(i + 0.25, j - 0.25, k - 0.25));
					particles.push_back(MarkerParticle(i + 0.25, j - 0.25, k + 0.25));
					particles.push_back(MarkerParticle(i + 0.25, j + 0.25, k - 0.25));
					particles.push_back(MarkerParticle(i + 0.25, j + 0.25, k + 0.25));
				}
			}
		}
	}
	for (MarkerParticle &p : particles) {
		p.x += (randomF() - 0.5)*0.1;
		p.y += (randomF() - 0.5)*0.1;
		p.z += (randomF() - 0.5)*0.1;
	}
}

//todo: Runge_Kutta here
void Advect_Particles(vector<MarkerParticle> &particles, Float *vx, Float *vy, Float *vz, int *mask){
	for (MarkerParticle &p : particles) {
		int ix = round(p.x), iy = round(p.y), iz = round(p.z);
		int k = ID(ix, iy, iz);
		if (mask[k] != SOLID) {
			p.x += Interpolation_In_Water_3D(vx, p.x, p.y, p.z, mask)*TIME_DELTA;
			p.y += Interpolation_In_Water_3D(vy, p.x, p.y, p.z, mask)*TIME_DELTA;
			p.z += Interpolation_In_Water_3D(vz, p.x, p.y, p.z, mask)*TIME_DELTA;
		}
	}
}
