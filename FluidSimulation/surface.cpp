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

void Add_Particles_Single_Cell(vector<MarkerParticle> &particles, int i, int j, int k) {
	const Float dxs[8] = { 0.25,0.25,0.25,0.25,0.75,0.75,0.75,0.75 };
	const Float dys[8] = { 0.25,0.25,0.75,0.75,0.25,0.25,0.75,0.75 };
	const Float dzs[8] = { 0.25,0.75,0.25,0.75,0.25,0.75,0.25,0.75 };
	for (int d = 0; d < 8; d++) {
		particles.push_back(MarkerParticle(i + dxs[d], j + dys[d], k + dzs[d]));
	}
	/*p.x += (randomF() - 0.5)*0.1;
	p.y += (randomF() - 0.5)*0.1;
	p.z += (randomF() - 0.5)*0.1;*/
}

void Init_Particles(vector<MarkerParticle> &particles, aryi &mask) {//should only be called once for initialization
	particles.clear();
	for (int i = 0; i < GRIDX; i++) {
		for (int j = 0; j < GRIDY; j++) {
			for (int k = 0; k < GRIDZ; k++) {
				if (mask(i, j, k) == WATER) {
					Add_Particles_Single_Cell(particles, i, j, k);
				}
			}
		}
	}
	LOGM("initial marker particles set\n");
}

void Add_Single_Particle(vector<MarkerParticle> &particles, aryi &mask, int i, int j, int k) {
	if (mask.is(i, j, k, WATER)) {
		Add_Particles_Single_Cell(particles, i, j, k);
	}
}

void Add_Particles(vector<MarkerParticle> &particles, aryi &mask, int x0, int x1, int y0, int y1, int z0, int z1) {
	for (int i = x0; i <= x1; i++) {
		for (int j = y0; j <= y1; j++) {
			for (int k = z0; k <= z1; k++) {
				if (mask.is(i, j, k, WATER)) {
					Add_Particles_Single_Cell(particles, i, j, k);
				}
			}
		}
	}
	LOGM("add water source particles\n");
}

Float clip(Float x, Float min, Float max) {
	x = std::max(min, x);
	x = std::min(max, x);
	return x;
}

void RK3(MarkerParticle &p, aryf &vx, aryf &vy, aryf &vz, aryi &mask, int step = 3) {
	const Float step_size[3] = {0.0, 1.0 / 2, 3.0 / 4};
	const Float weight[3] = {2.0 / 9, 3.0 / 9, 4.0 / 9};
	Float time_delta = TIME_DELTA / step;
	for (int t = 0; t < step; t++) {
		Float pvx = Interpolation_Water_Velocity(_X, vx, p.x, p.y, p.z, mask);
		Float pvy = Interpolation_Water_Velocity(_Y, vy, p.x, p.y, p.z, mask);
		Float pvz = Interpolation_Water_Velocity(_Z, vz, p.x, p.y, p.z, mask);
		Float delta_x = 0, delta_y = 0, delta_z = 0;
		for (int s = 0; s < 3; s++) {
			Float px = p.x + step_size[s] * pvx * time_delta;
			Float py = p.y + step_size[s] * pvy * time_delta;
			Float pz = p.z + step_size[s] * pvz * time_delta;
			Float pvx1 = Interpolation_Water_Velocity(_X, vx, px, py, pz, mask);
			Float pvy1 = Interpolation_Water_Velocity(_Y, vy, px, py, pz, mask);
			Float pvz1 = Interpolation_Water_Velocity(_Z, vz, px, py, pz, mask);
			delta_x += pvx1 * time_delta * weight[s];
			delta_y += pvy1 * time_delta * weight[s];
			delta_z += pvz1 * time_delta * weight[s];
		}
		p.x += delta_x, p.y += delta_y, p.z += delta_z;
		//XXX: find the nearest grid which is not solid
		p.x = clip(p.x, 1 + EPS, GRIDX - 1 - EPS);
		p.y = clip(p.y, 1 + EPS, GRIDY - 1 - EPS);
		p.z = clip(p.z, 1 + EPS, GRIDZ - 1 - EPS);
	}
}

//todo: Runge_Kutta here
void Advect_Particles(vector<MarkerParticle> &particles, aryf &vx, aryf &vy, aryf &vz, aryi &mask){
	for (MarkerParticle &p : particles) {
		int ix = floor(p.x), iy = floor(p.y), iz = floor(p.z);
		if (!mask.inside(ix, iy, iz)) {
			assert(false);
			printf("when advecting marker particle flying outside: %f %f %f\n", p.x, p.y, p.z);
		}
		else if (mask(ix,iy,iz) != SOLID) {
			RK3(p, vx, vy, vz, mask);
		}
	}
}
