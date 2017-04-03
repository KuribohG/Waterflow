#include "simulation_cubic.h"

SimulationCubic::SimulationCubic(void){
	vx = (Float*) calloc(GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z, sizeof(Float));
	vy = (Float*) calloc(GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z, sizeof(Float));
	vz = (Float*) calloc(GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z, sizeof(Float));
	vx0 = (Float*)calloc(GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z, sizeof(Float));
	vy0 = (Float*)calloc(GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z, sizeof(Float));
	vz0 = (Float*)calloc(GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z, sizeof(Float));
	density = (Float*) calloc(GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z, sizeof(Float));
	s = (Float*)calloc(GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z, sizeof(Float));
	mask = (int*) calloc(GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z, sizeof(int));
	printf("allocated\n");
	for (int i = 0; i < GRID_SIZE_X; i++) {
		for (int j = 0; j < GRID_SIZE_Y; j++) {
			mask[ID(i, j, 0)] = mask[ID(i, j, GRID_SIZE_Z - 1)] = SOLID;
		}
	}
	for (int i = 0; i < GRID_SIZE_X; i++) {
		for (int k = 0; k < GRID_SIZE_Z; k++) {
			mask[ID(i, 0, k)] = mask[ID(i, GRID_SIZE_Y - 1, k)] = SOLID;
		}
	}
	for (int j = 0; j < GRID_SIZE_Y; j++) {
		for (int k = 0; k < GRID_SIZE_Z; k++) {
			mask[ID(0, j, k)] = mask[ID(GRID_SIZE_X - 1, j, k)] = SOLID;
		}
	}
	printf("masked\n");
	for (int i = 0; i < GRID_SIZE_X; i++) {
		//cout << i << endl;
		for (int j = 0; j < GRID_SIZE_Y; j++) {
			for (int k = 0; k < GRID_SIZE_Z; k++) {
				//cout << i << " " << j << " " << k << endl;
				int t = ID(i, j, k);
				//vx[t] = randomF();
				//vy[t] = randomF()*10;
				//vz[t] = randomF()*10;
				vx[t] = 0;
				//vy[t] = 0;
				//vz[t] = 0;
				if (20 <= j &&j<= 30 && 30 <= k &&k<= 40) {
					vy[t] = 100;
				}
				//if (20 <= j &&j<= 30 && 40 <= k &&k<= 50) {
				//	vy[t] = -10;
				//}
				//vz[t] = 0;
				//if (100 <= k&&k <= 400) {
				//vy[t] = 100;
				//vz[t] = 50;
					//vz[t] = 250 - k;
				//}
				//else vy[t] = 0;
				//vy[t] = (j + 0.0);
				//density[t] = randomF();
				if (30 <= j&&j <= 40 && 30 <= k&&k <= 40) density[t] = 1.0;
				else density[t] = 0;
				//density[t] = ((j > GRID_SIZE_Y / 2)*0.5);
			}
		}
	}
	printf("randomed\n");
}

SimulationCubic::~SimulationCubic() {
	free(vx);
	free(vy);
	free(vz);
	free(vx0);
	free(vy0);
	free(vz0);
	free(density);
	free(mask);
	free(s);
}

Float Clip(Float &x, Float _min, Float _max) {
	if (x < _min) x = _min;
	if (x > _max) x = _max;
	return x;
}

Float Interpolation_3D(Float *f, Float x, Float y,Float z) {
	int x0 = floor(x), x1 = x0 + 1;
	int y0 = floor(y), y1 = y0 + 1;
	int z0 = floor(z), z1 = z0 + 1;
	Float s0 = x - x0, s1 = x1 - x;
	Float u0 = y - y0, u1 = y1 - y;
	Float v0 = z - z0, v1 = z1 - z;
	Float v=0, sw=0;
	if (validX(x0) && validY(y0) && validZ(z0)) v += s1*u1*v1*f[ID(x0, y0, z0)], sw += s1*u1*v1;
	if (validX(x0) && validY(y0) && validZ(z1)) v += s1*u1*v0*f[ID(x0, y0, z1)], sw += s1*u1*v0;
	if (validX(x0) && validY(y1) && validZ(z0)) v += s1*u0*v1*f[ID(x0, y1, z0)], sw += s1*u0*v1;
	if (validX(x0) && validY(y1) && validZ(z1)) v += s1*u0*v0*f[ID(x0, y1, z1)], sw += s1*u0*v0;
	if (validX(x1) && validY(y0) && validZ(z0)) v += s0*u1*v1*f[ID(x1, y0, z0)], sw += s0*u1*v1;
	if (validX(x1) && validY(y0) && validZ(z1)) v += s0*u1*v0*f[ID(x1, y0, z1)], sw += s0*u1*v0;
	if (validX(x1) && validY(y1) && validZ(z0)) v += s0*u0*v1*f[ID(x1, y1, z0)], sw += s0*u0*v1;
	if (validX(x1) && validY(y1) && validZ(z1)) v += s0*u0*v0*f[ID(x1, y1, z1)], sw += s0*u0*v0;
	if (sw == 0) return 0;
	return v / sw;
}

Float Neighbor_Sum6(Float *x, int i, int j, int k) {
	Float tmp = 0;
	if (i - 1 >= 0) tmp += x[ID(i - 1, j, k)];
	if (i + 1 < GRID_SIZE_X) tmp += x[ID(i + 1, j, k)];
	if (j - 1 >= 0) tmp += x[ID(i, j - 1, k)];
	if (j + 1 < GRID_SIZE_Y) tmp += x[ID(i, j + 1, k)];
	if (k - 1 >= 0) tmp += x[ID(i, j, k - 1)];
	if (k + 1 < GRID_SIZE_Z) tmp += x[ID(i, j, k + 1)];
	return tmp;
}

void Bound_Solid(int axis, Float *x) {
	for (int j = 0; j < GRID_SIZE_Y; j++) {
		for (int k = 0; k < GRID_SIZE_Z; k++) {
			x[ID(0, j, k)] = (axis == 0) ? -x[ID(1, j, k)] : x[ID(1, j, k)];
			x[ID(GRID_SIZE_X - 1, j, k)] = (axis == 0) ? -x[ID(GRID_SIZE_X - 2, j, k)] : x[ID(GRID_SIZE_X - 2, j, k)];
		}
	}
	for (int i = 0; i < GRID_SIZE_X; i++) {
		for (int k = 0; k < GRID_SIZE_Z; k++) {
			x[ID(i, 0, k)] = (axis == 1) ? -x[ID(i, 1, k)] : x[ID(i, 1, k)];
			x[ID(i, GRID_SIZE_Y - 1, k)] = (axis == 1) ? -x[ID(i, GRID_SIZE_Y - 2, k)] : x[ID(i, GRID_SIZE_Y - 2, k)];
		}
	}
	for (int i = 0; i < GRID_SIZE_X; i++) {
		for (int j = 0; j < GRID_SIZE_Y; j++) {
			x[ID(i, j, 0)] = (axis == 2) ? -x[ID(i, j, 1)] : x[ID(i, j, 1)];
			x[ID(i, j, GRID_SIZE_Z - 1)] = (axis == 2) ? -x[ID(i, j, GRID_SIZE_Z - 2)] : x[ID(i, j, GRID_SIZE_Z - 2)];
		}
	}
}

void SimulationCubic::Linear_Solve(int axis, Float *x, Float *x0, Float a, Float wtsum, int iter) {
	//this is so-called "Semi-Lagrangian" Advectoin
	//we shall change it to MacCormack Advection in further development
	//"hyperbolic equation": du/dt + a*du/dx = 0
	memcpy(x, x0, sizeof(Float)*GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z);
	Float rcp = 1.0 / wtsum;
	for (int t = 0; t < iter; t++) {
		for (int i = 0; i < GRID_SIZE_X; i++) {
			for (int j = 0; j < GRID_SIZE_Y; j++) {
				for (int k = 0; k < GRID_SIZE_Z; k++) {
					if (mask[ID(i, j, k)] == WATER) {
						x[ID(i, j, k)] = (x0[ID(i, j, k)] + a*Neighbor_Sum6(x, i, j, k)) *rcp;
					}
				}
			}
		}
	}
	/*for (int j = 240; j < 255; j++) {
		for (int k = 240; k < 255; k++) {
			printf("%f ", x[ID(2, j, k)]);
		}
		printf("\n");
	}printf("\n");*/
	Bound_Solid(axis, x);
}

void SimulationCubic::Diffuse(int axis, Float *x, Float *x0, Float diff, Float dt, int iter) {//diffusion is blurring
	//solve x from x0
	Float a = dt*diff;
	printf("diffusion parameter: %d\n", iter);
	Linear_Solve(axis, x, x0, a, 6 * a + 1, iter);
}

void SimulationCubic::Project(Float *vx, Float *vy, Float *vz, Float *s, Float *div) {
	/*based on Helmholtz decomposition, a vector field can be resolved into the sum of
	a mass-conserving field and a gradient field, and we hope the velocity field is mass-conserving*/
	for (int i = 0; i < GRID_SIZE_X; i++) {
		for (int j = 0; j < GRID_SIZE_Y; j++) {
			for (int k = 0; k < GRID_SIZE_Z; k++) {
				if (mask[ID(i, j, k)] == WATER) {
					Float tmp = 0;
					if (i + 1 < GRID_SIZE_X) tmp += vx[ID(i + 1, j, k)];
					if (i - 1 >= 0) tmp -= vx[ID(i - 1, j, k)];
					if (j + 1 < GRID_SIZE_Y) tmp += vy[ID(i, j + 1, k)];
					if (j - 1 >= 0) tmp -= vy[ID(i, j - 1, k)];
					if (k + 1 < GRID_SIZE_Z) tmp += vz[ID(i, j, k + 1)];
					if (k - 1 >= 0) tmp -= vz[ID(i, j, k - 1)];
					div[ID(i, j, k)] = -0.5*tmp;//it's divergence of velocity field
					s[ID(i, j, k)] = 0;
				}
			}
		}
	}
	Bound_Solid(-1, div);
	Linear_Solve(-1, s, div, 1, 6+1, LINSOLVER_ITER);
	//basically, s is some "blurred" divergence
	//we substract the gradient of divergence from velocity field
	for (int i = 0; i < GRID_SIZE_X; i++) {
		for (int j = 0; j < GRID_SIZE_Y; j++) {
			for (int k = 0; k < GRID_SIZE_Z; k++) {
				if (mask[ID(i, j, k)] == WATER) {
					Float tmp = 0;
					if (i + 1 < GRID_SIZE_X) tmp += s[ID(i + 1, j, k)];
					if (i - 1 >= 0) tmp -= s[ID(i - 1, j, k)];
					vx[ID(i, j, k)] -= 0.5*tmp;
					tmp = 0;
					if (j + 1 < GRID_SIZE_Y) tmp += s[ID(i, j + 1, k)];
					if (j - 1 >= 0) tmp -= s[ID(i, j - 1, k)];
					vy[ID(i, j, k)] -= 0.5*tmp;
					tmp = 0;
					if (k + 1 < GRID_SIZE_Z) tmp += s[ID(i, j, k + 1)];
					if (k - 1 >= 0) tmp -= s[ID(i, j, k - 1)];
					vz[ID(i, j, k)] -= 0.5*tmp;
				}
			}
		}
	}
	Bound_Solid(0, vx);
	Bound_Solid(1, vy);
	Bound_Solid(2, vz);
}

void SimulationCubic::Advect(int axis, Float *density, Float *density0, Float *vx, Float *vy, Float *vz) {
	Float x, y, s0, s1, t0, t1;
	int i0, i1, j0, j1;
	for (int i = 0; i < GRID_SIZE_X; i++) {
		for (int j = 0; j < GRID_SIZE_Y; j++) {
			for (int k = 0; k < GRID_SIZE_Z; k++) {
				if (mask[ID(i, j, k)] == WATER) {
					Float x = i - TIME_DELTA*vx[ID(i, j, k)];
					//Float x = i;
					Float y = j - TIME_DELTA*vy[ID(i, j, k)];
					Float z = k - TIME_DELTA*vz[ID(i, j, k)];
					//printf("%f %f %f\n", TIME_DELTA*vx[ID(i, j, k)], TIME_DELTA*vy[ID(i, j, k)], TIME_DELTA*vz[ID(i, j, k)]);
					x = Clip(x, 0.5, GRID_SIZE_X + 0.5);
					y = Clip(y, 0.5, GRID_SIZE_Y + 0.5);
					z = Clip(z, 0.5, GRID_SIZE_Z + 0.5);
					density[ID(i, j, k)] = Interpolation_3D(density0, x, y, z);
				}
			}
		}
	}
	Bound_Solid(axis, density);
}


void SimulationCubic::Step_Time(void){
	printf("cubic step time \n");

	//velocity-evolution
	//Diffuse(0, vx0, vx, viscosity, TIME_DELTA, LINSOLVER_ITER);
	Diffuse(1, vy0, vy, viscosity, TIME_DELTA, LINSOLVER_ITER);
	Diffuse(2, vz0, vz, viscosity, TIME_DELTA, LINSOLVER_ITER);
	//now vx0, vy0, vz0 are "blurred" velocities
	Project(vx0, vy0, vz0, vx, vy);
	//swap(vx, vx0);
	//swap(vy, vy0);
	//swap(vz, vz0);
	cout << "advect x:\n";
	//Advect(0, vx, vx0, vx0, vy0, vz0);
	Advect(1, vy, vy0, vx0, vy0, vz0);
	Advect(2, vz, vz0, vx0, vy0, vz0);
	Project(vx, vy, vz, vx0, vy0);


	//diff = 1.0*GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z*200;
	//diff = 1e6;

	//density evolution
	Diffuse(-1, s, density, diff, TIME_DELTA, LINSOLVER_ITER);
	//now stored in s
	Float sm = 0;
	for (int j = 0; j < GRID_SIZE_Y; j++) {
		for (int k = 0;k<GRID_SIZE_Z; k++) {
			//printf("%f ", s[ID(2, j, k)]);
			sm += s[ID(2, j, k)];
		}
		//printf("\n");
	}printf("%f \n",sm);
	Advect(-1, density, s, vx, vy, vz);
	//now stored in density
	//memcpy(density, s, sizeof(Float)*GRID_SIZE_X*GRID_SIZE_Y*GRID_SIZE_Z);
	/*for (int j = 0; j < 10; j++) {
		for (int k = 0; k < 10; k++) {
			printf("%f ", density[ID(2, j, k)]);
		}
		printf("\n");
	}printf("\n");*/
	//printf("density: %lf\n", density[0]);
}
