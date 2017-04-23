#include "simulation_cubic.h"

//void Runge_Kutta(int i, int j, int k, Float delta, int iter,
//       Float vx[], Float vy[], Float vz[], Float &x, Float &y, Float &z);

SimulationCubic::SimulationCubic(void){
	vx.init(GRIDX + 1, GRIDY, GRIDZ); vx0.init(GRIDX + 1, GRIDY, GRIDZ);
	vy.init(GRIDX, GRIDY + 1, GRIDZ); vy0.init(GRIDX, GRIDY + 1, GRIDZ);
	vz.init(GRIDX, GRIDY, GRIDZ + 1); vz0.init(GRIDX, GRIDY, GRIDZ + 1);
	p.init(GRIDX, GRIDY, GRIDZ); p0.init(GRIDX, GRIDY, GRIDZ);
	mask.init(GRIDX, GRIDY, GRIDZ);
    for (int i = 0; i < GRIDX; i++) {
        for (int j = 0; j < GRIDY; j++) {
            mask(i, j, 0) = mask(i, j, GRIDZ - 1) = SOLID;
        }
    }
    for (int i = 0; i < GRIDX; i++) {
        for (int k = 0; k < GRIDZ; k++) {
            mask(i, 0, k) = mask(i, GRIDY - 1, k) = SOLID;
        }
    }
    for (int j = 0; j < GRIDY; j++) {
        for (int k = 0; k < GRIDZ; k++) {
            mask(0, j, k) = mask(GRIDX - 1, j, k) = SOLID;
        }
    }
    LOGM("solid mask set\n");
    for (int i = 0; i < GRIDX; i++) {
        for (int j = 0; j < GRIDY; j++) {
            for (int k = 0; k < GRIDZ; k++) {
				if ((k == GRIDZ - 2 || j >= 30) && mask(i, j, k) != SOLID) mask(i, j, k) = AIR;
            }
        }
    }
	LOGM("velocity set\n");
	Place_Particles(particles, mask);
	LOGM("marker particles set\n");
}

SimulationCubic::~SimulationCubic() {
	vx.clear(), vx0.clear(), vy.clear(), vy0.clear(), vz.clear(), vz0.clear();
	p.clear(), p0.clear();
	mask.clear();
}

Float Clip(Float &x, Float _min, Float _max) {
    if (x < _min) x = _min;
    if (x > _max) x = _max;
    return x;
}

void Bound_Solid(int axis, aryf &x) {
    for (int j = 0; j < x.m; j++) {
        for (int k = 0; k < x.w; k++) {
            x(0, j, k) = (axis == 0) ? -x(1, j, k) : x(1, j, k);
            x(x.n - 1, j, k) = (axis == 0) ? -x(x.n - 2, j, k) : x(x.n - 2, j, k);
        }
    }
    for (int i = 0; i < x.n; i++) {
        for (int k = 0; k < x.w; k++) {
            x(i, 0, k) = (axis == 1) ? -x(i, 1, k) : x(i, 1, k);
            x(i, x.m - 1, k) = (axis == 1) ? -x(i, x.m - 2, k) : x(i, x.m - 2, k);
        }
    }
    for (int i = 0; i < x.n; i++) {
        for (int j = 0; j < x.m; j++) {
            x(i, j, 0) = (axis == 2) ? -x(i, j, 1) : x(i, j, 1);
            x(i, j, x.w - 1) = (axis == 2) ? -x(i, j, x.w - 2) : x(i, j, x.w - 2);
        }
    }
}

void Bound_Surface(aryf &pressure, aryi &mask) {
	for (int i = 0; i < GRIDX; i++) {
		for (int j = 0; j < GRIDY; j++) {
			for (int k = 0; k < GRIDZ; k++) {
				if (mask(i,j,k) == AIR) pressure(i,j,k) = 0;
			}
		}
	}
}

void SimulationCubic::Apply_External_Forces(void) {
	for (int i = 0; i < GRIDX; i++) {
		for (int j = 0; j < GRIDY; j++) {
			for (int k = 0; k < GRIDZ; k++) {
				int t = mask(i, j, k);
				if (t == WATER) {
					vz(i, j, k) += -g*TIME_DELTA;
					vz(i, j, k + 1) += -g*TIME_DELTA;
				}
				if (i < GRIDX&&t == AIR&&mask(i + 1, j, k) == AIR) vx(i + 1, j, k) = 0;
				if (j < GRIDY&&t == AIR&&mask(i, j + 1, k) == AIR) vy(i, j + 1, k) = 0;
				if (k < GRIDZ&&t == AIR&&mask(i, j, k + 1) == AIR) vz(i, j, k + 1) = 0;
			}
		}
	}
}


void SimulationCubic::Linear_Solve(int axis, aryf &x, aryf &x0, Float a, Float wtsum, int iter) {
    //this is so-called "Semi-Lagrangian" Advectoin
    //we shall change it to MacCormack Advection in further development
    //"hyperbolic equation": du/dt + a*du/dx = 0
	x.copy(x0);
    Float rcp = 1.0 / wtsum;
    for (int t = 0; t < iter; t++) {
        for (int i = 0; i < GRIDX; i++) {
            for (int j = 0; j < GRIDY; j++) {
                for (int k = 0; k < GRIDZ; k++) {
                    if (mask(i, j, k) == WATER) {
                        x(i, j, k) = (x0(i, j, k) + a*Neighbor_Sum6(x, i, j, k)) *rcp;
                    } 
                }
            }
        }
    }
    /*for (int j = 240; j < 255; j++) {
        for (int k = 240; k < 255; k++) {
            LOGM("%f ", x[ID(2, j, k)]);
        }
        LOGM("\n");
    }LOGM("\n");*/
    Bound_Solid(axis, x);
}

/*void SimulationCubic::Diffuse(int axis, Float *x, Float *x0, Float diff, Float dt, int iter) {//diffusion is blurring
    //solve x from x0
    Float a = dt*diff;
    LOGM("diffusion parameter: %d\n", iter);
    Linear_Solve(axis, x, x0, a, 6 * a + 1, iter);
}*/

void SimulationCubic::Calc_Divergence(aryf &vx, aryf &vy, aryf &vz, aryf &div) {
	for (int i = 0; i < GRIDX; i++) {
		for (int j = 0; j < GRIDY; j++) {
			for (int k = 0; k < GRIDZ; k++) {
				if (mask(i, j, k) == WATER) {
					div(i, j, k) =
						+vx(i + 1, j, k) - vx(i, j, k)
						+ vy(i, j + 1, k) - vy(i, j, k)
						+ vz(i, j, k + 1) - vz(i, j, k);
				}
				/*if (mask[ID(i, j, k)] == WATER) {
					div[ID(i, j, k)] = -0.5*tmp;//it's divergence of velocity field
				}*/
			}
		}
	}
}

void SimulationCubic::Project(aryf &vx,aryf &vy,aryf &vz,aryf &p,aryf &div) {
    /*based on Helmholtz decomposition, a vector field can be resolved into the sum of
    a mass-conserving field and a gradient field, and we hope the velocity field is mass-conserving*/
    /*for (int i = 0; i < GRIDX; i++) {
        for (int j = 0; j < GRIDY; j++) {
            for (int k = 0; k < GRIDZ; k++) {
                if (mask[ID(i, j, k)] == WATER) {
                    Float tmp = 0;
                    if (i + 1 < GRIDX) tmp += vx[ID(i + 1, j, k)];
                    if (i - 1 >= 0) tmp -= vx[ID(i - 1, j, k)];
                    if (j + 1 < GRIDY) tmp += vy[ID(i, j + 1, k)];
                    if (j - 1 >= 0) tmp -= vy[ID(i, j - 1, k)];
                    if (k + 1 < GRIDZ) tmp += vz[ID(i, j, k + 1)];
                    if (k - 1 >= 0) tmp -= vz[ID(i, j, k - 1)];
                    div[ID(i, j, k)] = -0.5*tmp;//it's divergence of velocity field
                    s[ID(i, j, k)] = 0;
                }
            }
        }
    }*/
	Calc_Divergence(vx, vy, vz, div);
	p.set(0);
	//memset(s, 0, sizeof(s[0])*GRIDX*GRIDY*GRIDZ);
	//Bound_Solid(-1, div);
	
	//actually, p solved here is -deltaT*p, look at [Robert Bridson, p54]
	Linear_Solve(-1, p, div, 1, 6 + 1, LINSOLVER_ITER);

	Bound_Surface(p, mask);
	for (int i = 0; i < GRIDX; i++) {
		for (int j = 0; j < GRIDY; j++) {
			for (int k = 0; k < GRIDZ; k++) {
				if (i > 0) vx(i, j, k) += (p(i, j, k) - p(i - 1, j, k));
				if (j > 0) vy(i, j, k) += (p(i, j, k) - p(i, j - 1, k));
				if (k > 0) vz(i, j, k) += (p(i, j, k) - p(i, j, k - 1));
			}
		}
	}

	/*for (int i = 0; i < GRIDX*GRIDY*GRIDZ; i++) s[i] *= 2;
    //basically, s is some "blurred" divergence
    //we substract the gradient of divergence from velocity field
    for (int i = 0; i < GRIDX; i++) {
        for (int j = 0; j < GRIDY; j++) {
            for (int k = 0; k < GRIDZ; k++) {
                if (mask[ID(i, j, k)] == WATER) {
                    Float tmp = 0;
                    if (i + 1 < GRIDX) tmp += s[ID(i + 1, j, k)];
                    if (i - 1 >= 0) tmp -= s[ID(i - 1, j, k)];
                    vx[ID(i, j, k)] -= 0.5*tmp;
                    tmp = 0;
                    if (j + 1 < GRIDY) tmp += s[ID(i, j + 1, k)];
                    if (j - 1 >= 0) tmp -= s[ID(i, j - 1, k)];
                    vy[ID(i, j, k)] -= 0.5*tmp;
                    tmp = 0;
                    if (k + 1 < GRIDZ) tmp += s[ID(i, j, k + 1)];
                    if (k - 1 >= 0) tmp -= s[ID(i, j, k - 1)];
                    vz[ID(i, j, k)] -= 0.5*tmp;
                }
            }
        }
    }*/
    Bound_Solid(0, vx);
    Bound_Solid(1, vy);
    Bound_Solid(2, vz);
}

void SimulationCubic::Runge_Kutta(int i, int j, int k, Float delta, int iter,
	aryf &vx, aryf &vy, aryf &vz, Float &x, Float &y, Float &z) {
	assert(iter >= 1);

	delta /= iter;

	//int id = ID(i, j, k);
	x = i + 0.5 - delta * vx(i, j, k);
	y = j + 0.5 - delta * vy(i, j, k);
	z = k + 0.5 - delta * vz(i, j, k);

	for (int _ = 1; _ < iter; _++) {
		Float nx = x - delta * Interpolation_In_Water_3D(vx, x, y, z, mask);
		Float ny = y - delta * Interpolation_In_Water_3D(vy, x, y, z, mask);
		Float nz = z - delta * Interpolation_In_Water_3D(vz, x, y, z, mask);
		x = nx, y = ny, z = nz;
	}
}

void SimulationCubic::Advect(int axis, aryf &f, aryf &f0, aryf &vx, aryf &vy, aryf &vz) {
    Float x, y, s0, s1, t0, t1;
    int i0, i1, j0, j1;
    for (int i = 0; i < GRIDX; i++) {
        for (int j = 0; j < GRIDY; j++) {
            for (int k = 0; k < GRIDZ; k++) {
                if (mask(i, j, k) == WATER) {
                    Float x, y, z;
                    Runge_Kutta(i, j, k, TIME_DELTA, 2, vx, vy, vz, x, y, z);
                    //LOGM("%f %f %f\n", TIME_DELTA*vx[ID(i, j, k)], TIME_DELTA*vy[ID(i, j, k)], TIME_DELTA*vz[ID(i, j, k)]);
                    x = Clip(x, 0.5, GRIDX + 0.5);
                    y = Clip(y, 0.5, GRIDY + 0.5);
                    z = Clip(z, 0.5, GRIDZ + 0.5);
                    f(i, j, k) = Interpolation_In_Water_3D(f0, x, y, z,mask);
                }
            }
        }
    }
	Bound_Solid(axis, f);
}

void swap(aryf &a, aryf &b) {
	aryf c;
	c = a; a = b; b = c;
}

void SimulationCubic::Step_Time(void){
    LOGM("cubic step time \n");

    //velocity-evolution
	Apply_External_Forces();

    //Diffuse(0, vx0, vx, viscosity, TIME_DELTA, LINSOLVER_ITER);
    //Diffuse(1, vy0, vy, viscosity, TIME_DELTA, LINSOLVER_ITER);
    //Diffuse(2, vz0, vz, viscosity, TIME_DELTA, LINSOLVER_ITER);

	printf("%p %p %f %f\n", vx.f, vx0.f, vx.f[0], vx0.f[0]);
	swap(vx, vx0);
	swap(vy, vy0);
	swap(vz, vz0);
	printf("%p %p %f %f\n", vx.f, vx0.f, vx.f[0], vx0.f[0]);


    //now vx0, vy0, vz0 are "blurred" velocities
    //Project(vx0, vy0, vz0, vx, vy);

    //swap(vx, vx0);
    //swap(vy, vy0);
    //swap(vz, vz0);
    LOGM( "advect x:\n");
    Advect(0, vx, vx0, vx0, vy0, vz0);
    Advect(1, vy, vy0, vx0, vy0, vz0);
    Advect(2, vz, vz0, vx0, vy0, vz0);


	Project(vx, vy, vz, p, p0);
	//Calc_Divergence(vx, vy, vz, s);


    //diff = 1.0*GRIDX*GRIDY*GRIDZ*200;
    //diff = 1e6;

    //density evolution
    //Diffuse(-1, s, density, diff, TIME_DELTA, LINSOLVER_ITER);

	Advect_Particles(particles, vx, vy, vz, mask);
	Mark_Water_By(particles, mask);


	Calc_Divergence(vx, vy, vz, p);
	/*Calc_Divergence(vx, vy, vz, vx0);
	memset(s, 0, sizeof(s[0])*GRIDX*GRIDY*GRIDZ);
	Bound_Solid(-1, vx0);
	Linear_Solve(-1, s, vx0, 1, 6 + 1, LINSOLVER_ITER);
	Bound_Surface(s, mask);
	for (int i = 0; i < GRIDX*GRIDY*GRIDZ; i++) s[i] -= vx0[i];*/
}


