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
		int x = floor(p.x), y = floor(p.y), z = floor(p.z);
		if (!mask.inside(x, y, z)) {
			LOGM("when marking marker particle flying outside: %f %f %f\n", p.x, p.y, p.z);
		}
		else if (mask(x, y, z) != SOLID) {
			//LOGM("%f %f %f %d %d %d\n", p.x, p.y, p.z, x, y, z);
			mask(x, y, z) = WATER;
		}
	}
}

void Place_Particles(vector<MarkerParticle> &particles, aryi &mask) {
	const Float dxs[8] = { 0.25,0.25,0.25,0.25,0.75,0.75,0.75,0.75 };
	const Float dys[8] = { 0.25,0.25,0.75,0.75,0.25,0.25,0.75,0.75 };
	const Float dzs[8] = { 0.25,0.75,0.25,0.75,0.25,0.75,0.25,0.75 };
	particles.clear();
	for (int i = 0; i < GRIDX; i++) {
		for (int j = 0; j < GRIDY; j++) {
			for (int k = 0; k < GRIDZ; k++) {
				if (mask(i, j, k) == WATER) {
					for (int d = 0; d < 8; d++) {
						particles.push_back(MarkerParticle(i + dxs[d], j + dys[d], k + dzs[d]));
					}
				}
			}
		}
	}
	for (MarkerParticle &p : particles) {
		//p.x += (randomF() - 0.5)*0.1;
		//p.y += (randomF() - 0.5)*0.1;
		//p.z += (randomF() - 0.5)*0.1;
	}
	LOGM("marker particles set\n");
}

//todo: Runge_Kutta here
void Advect_Particles(vector<MarkerParticle> &particles, aryf &vx, aryf &vy, aryf &vz, aryi &mask){
	for (MarkerParticle &p : particles) {
		int ix = floor(p.x), iy = floor(p.y), iz = floor(p.z);
		if (!mask.inside(ix, iy, iz)) {
			//assert(false);
			printf("when advecting marker particle flying outside: %f %f %f\n", p.x, p.y, p.z);
		}
		else if (mask(ix,iy,iz) != SOLID) {
			Float z0 = p.z;
			Float pvx = Interpolation_Water_Velocity(_X, vx, p.x, p.y, p.z, mask, false);
			Float pvy = Interpolation_Water_Velocity(_Y, vy, p.x, p.y, p.z, mask, false);
			Float pvz = Interpolation_Water_Velocity(_Z, vz, p.x, p.y, p.z, mask, false);
			//if (ix == GRIDX / 2) LOGM("advect particle: %f %f %f %f\n", p.x, p.y, p.z, pvz);
			//assert(pvz <= 0);
			Float x1 = p.x + pvx*TIME_DELTA;
			Float y1 = p.y + pvy*TIME_DELTA;
			Float z1 = p.z + pvz*TIME_DELTA;
			if (mask.is(floor(x1), floor(y1), floor(z1), WATER) || mask.is(floor(x1), floor(y1), floor(z1), AIR)) {
				p.x = x1;
				p.y = y1;
				p.z = z1;
			}
			//p.x += pvx*TIME_DELTA;
			//p.y += pvy*TIME_DELTA;
			//p.z += pvz*TIME_DELTA;
		}
	}
}
