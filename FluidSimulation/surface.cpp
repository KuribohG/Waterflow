#include "surface.h"


void Mark_Water_By(vector<MarkerParticle> &particles, aryi &mask) {
	for (int i = 0; i < GRIDX; i++) {
		for (int j = 0; j < GRIDY; j++) {
			for (int k = 0; k < GRIDZ; k++) {
				if (mask(i, j, k) != SOLID) mask(i, j, k) = AIR;
			}
		}
	}
	for (MarkerParticle &p : particles) {
		int x = round(p.x), y = round(p.y), z = round(p.z);
		if (!mask.inside(x, y, z)) {
			LOGM("when marking marker particle flying outside: %f %f %f\n", p.x, p.y, p.z);
		}
		else if (mask(x, y, z) != SOLID) mask(x, y, z) = WATER;
	}
}

void Place_Particles(vector<MarkerParticle> &particles, aryi &mask) {
	particles.clear();
	for (int i = 0; i < GRIDX; i++) {
		for (int j = 0; j < GRIDY; j++) {
			for (int k = 0; k < GRIDZ; k++) {
				if (mask(i, j, k) == WATER) {
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
void Advect_Particles(vector<MarkerParticle> &particles, aryf &vx, aryf &vy, aryf &vz, aryi &mask){
	for (MarkerParticle &p : particles) {
		int ix = round(p.x), iy = round(p.y), iz = round(p.z);
		if (!mask.inside(ix, iy, iz)) {
			LOGM("when advecting marker particle flying outside: %f %f %f\n", p.x, p.y, p.z);
		}
		else if (mask(ix,iy,iz) != SOLID) {
			p.x += Interpolation_In_Water_3D(vx, p.x, p.y, p.z, mask)*TIME_DELTA;
			p.y += Interpolation_In_Water_3D(vy, p.x, p.y, p.z, mask)*TIME_DELTA;
			p.z += Interpolation_In_Water_3D(vz, p.x, p.y, p.z, mask)*TIME_DELTA;
		}
	}
}
