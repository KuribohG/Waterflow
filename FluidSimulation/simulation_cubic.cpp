#include "simulation_cubic.h"

void Print_Velocity(const aryf &vx, const aryf &vy, const aryf &vz, const aryi &mask) {
	printf("print velocity: \n");
	Float vz0 = INF;
	for (int i = GRIDX / 2; i <= GRIDX / 2; i++) {
		for (int j = 0; j < GRIDY; j++) {
			for (int k = 0; k < GRIDZ; k++) {
				if (mask.is(i, j, k, WATER)) {
					Float v1 = Interpolation_Water_Velocity(_Y, vy, i + 0.5, j + 0.5, k + 0.5, mask, false);
					Float v2 = Interpolation_Water_Velocity(_Z, vz, i + 0.5, j + 0.5, k + 0.5, mask, false);
					//printf("cell-center velocity: %f %f %f %f %f\n", i + 0.5, j + 0.5, k + 0.5, v1, v2);
					//if (vz0 > 100) vz0 = v2;
					//assert(vz0 == v2);
				}
				Float cvz = vz.get(i, j, k);
				if (cvz != 0) printf("on-surface particle: %d %d %d %f\n", i, j, k, cvz);
			}
		}
	}
}

/*void SimulationCubic::Extrapolate(aryf &f){
	for(int i=0;i<)
}*/

SimulationCubic::SimulationCubic(void){
	temp_dis.init(GRIDX + 1, GRIDY + 1, GRIDZ + 1);
	vx.init(GRIDX + 1, GRIDY, GRIDZ); vx0.init(GRIDX + 1, GRIDY, GRIDZ);
	vy.init(GRIDX, GRIDY + 1, GRIDZ); vy0.init(GRIDX, GRIDY + 1, GRIDZ);
	vz.init(GRIDX, GRIDY, GRIDZ + 1); vz0.init(GRIDX, GRIDY, GRIDZ + 1);
	p.init(GRIDX, GRIDY, GRIDZ); p0.init(GRIDX, GRIDY, GRIDZ);
	mask.init(GRIDX, GRIDY, GRIDZ);
	solver.init(GRIDX, GRIDY, GRIDZ);
    for (int i = 0; i < GRIDX; i++) {
        for (int j = 0; j < GRIDY; j++) {
			mask(i, j, 0) = SOLID; mask(i, j, GRIDZ - 1) = SOLID;
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
				if (mask.is(i, j, k, SOLID)) continue;
				mask(i, j, k) = AIR;
				//if (2 <= i&&i <= 8 && 20 <= j&&j <= 50 && !mask.is(i, j, k, SOLID)) mask(i, j, k) = WATER;
				//if (j == 30 && 40 <= k&&k <= 50) mask(i, j, k) = WATER;
				//if (20 <= j&&j <= 30 && 50 <= k&&k <= 60) mask(i, j, k) = WATER;
				if (29 <= j&&j <= 29 && 40 <= k&&k <= 50) mask(i, j, k) = WATER;
				//mask(i, 50, 50) = WATER;
				//if ((k == GRIDZ - 2 || j >= 40 || j <= 10) && mask(i, j, k) != SOLID) mask(i, j, k) = AIR;
				//if (!mask.is(i, j, k, SOLID)) mask(i, j, k) = AIR;
            }
        }
    }
	//mask(2, 30, 30) = WATER;
	LOGM("velocity set\n");
	Place_Particles(particles, mask);
}

SimulationCubic::~SimulationCubic() {
	temp_dis.clear();
	vx.clear(), vx0.clear(), vy.clear(), vy0.clear(), vz.clear(), vz0.clear();
	p.clear(), p0.clear();
	mask.clear();
}

Float Clip(Float &x, Float _min, Float _max) {
    if (x < _min) x = _min;
    if (x > _max) x = _max;
    return x;
}

void SimulationCubic::Bound_Solid(void) {
	for (int i = 0; i < GRIDX; i++) {
		for (int j = 0; j < GRIDY; j++) {
			for (int k = 0; k < GRIDZ; k++) {
				if (mask.is(i, j, k, SOLID)) {
					vx(i, j, k) = vx(i + 1, j, k) = 0;
					vy(i, j, k) = vy(i, j + 1, k) = 0;
					vz(i, j, k) = vz(i, j, k + 1) = 0;
				}
			}
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
				if (t == WATER || mask.is(i + 1, j, k, WATER));
				else vx(i, j, k) = 0;
				if (t == WATER || mask.is(i, j + 1, k, WATER));
				else vy(i, j, k) = 0;
				if (t == WATER || mask.is(i, j, k + 1, WATER)) {
					if(i==GRIDX/2) LOGM("add gravity: %d %d %d\n", i, j, k);
					vz(i, j, k) += -g*TIME_DELTA;
				}
				else vz(i, j, k) = 0;
				//if (i == 2 && (t == WATER || mask.is(i, j, k + 1, WATER))) LOGM("(%d %d)", j, k);
			}
		}
	}
}








void SimulationCubic::Calc_Divergence(aryf &vx, aryf &vy, aryf &vz, aryf &div) {
	div.set(0);
	for (int i = 0; i < GRIDX; i++) {
		for (int j = 0; j < GRIDY; j++) {
			for (int k = 0; k < GRIDZ; k++) {
				if (mask(i, j, k) == WATER) {
					div(i, j, k) =
						+vx(i + 1, j, k) - vx(i, j, k)
						+ vy(i, j + 1, k) - vy(i, j, k)
						+ vz(i, j, k + 1) - vz(i, j, k);
				}
			}
		}
	}
}

void SimulationCubic::Project(aryf &vx,aryf &vy,aryf &vz,aryf &p,aryf &div) {
	LOGM("project\n");
	solver.Solve_Pressure(vx, vy, vz, mask);
	solver.Send_Back_To(p);

	//actually, p solved here is -deltaT*p, look at [Robert Bridson, p54]
	//Linear_Solve(-1, p, div, 1, 6 + 1, LINSOLVER_ITER);

	//-+ //todo: apply bound condition before, not after linear_solve
	//-+ //todo: now this routine's solid bound fails, try to fix it
	double scale = TIME_DELTA / DENSITY;
	Bound_Surface(p, mask);
	for (int i = 0; i < GRIDX; i++) {
		for (int j = 0; j < GRIDY; j++) {
			for (int k = 0; k < GRIDZ; k++) {
				bool s0 = mask.is(i, j, k, SOLID);
				bool si = mask.is(i - 1, j, k, SOLID);
				bool sj = mask.is(i, j - 1, k, SOLID);
				bool sk = mask.is(i, j, k - 1, SOLID);
				if (i > 0 && !s0 && !mask.is(i - 1, j, k, SOLID)) vx(i, j, k) -= (p(i, j, k) - p(i - 1, j, k))*scale;
				else vx(i, j, k) = 0;
				if (j > 0 && !s0 && !mask.is(i, j - 1, k, SOLID)) vy(i, j, k) -= (p(i, j, k) - p(i, j - 1, k))*scale;
				else vy(i, j, k) = 0;
				if (k > 0 && !s0 && !mask.is(i, j, k - 1, SOLID)) vz(i, j, k) -= (p(i, j, k) - p(i, j, k - 1))*scale;
				else vz(i, j, k) = 0;
			}
		}
	}


}

void SimulationCubic::Runge_Kutta(int i, int j, int k, Float delta, int iter,
	const aryf &vx, const aryf &vy, const aryf &vz, Float &x, Float &y, Float &z) {
	assert(iter >= 1);

	delta /= iter;

	//int id = ID(i, j, k);
	x = i + 0.5 - delta * vx.get(i, j, k);
	y = j + 0.5 - delta * vy.get(i, j, k);
	z = k + 0.5 - delta * vz.get(i, j, k);

	for (int _ = 1; _ < iter; _++) {
		Float nx = x - delta*Interpolation_Water_Velocity(_X, vx, x, y, z, mask, true);
		Float ny = y - delta*Interpolation_Water_Velocity(_Y, vy, x, y, z, mask, true);
		Float nz = z - delta*Interpolation_Water_Velocity(_Z, vz, x, y, z, mask, true);
		x = nx, y = ny, z = nz;
	}
}

void SimulationCubic::Advect_Velocity(int axis, aryf &f, const aryf &f0, const aryf &vx, const aryf &vy, const aryf &vz) {
	LOGM("advect along: %d\n", axis);
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
					f(i, j, k) = Interpolation_Water_Velocity(axis, f0, x, y, z, mask, true);
                }
            }
        }
    }
	//Bound_Solid(axis, f);
}

void swap(aryf &a, aryf &b) {
	aryf c;
	c = a; a = b; b = c;
}

void SimulationCubic::Step_Time(void){
    printf("cubic step time \n");

    //velocity-evolution
	Apply_External_Forces();

    swap(vx, vx0);
    swap(vy, vy0);
    swap(vz, vz0);
    Advect_Velocity(0, vx, vx0, vx0, vy0, vz0);
	Advect_Velocity(1, vy, vy0, vx0, vy0, vz0);
	Advect_Velocity(2, vz, vz0, vx0, vy0, vz0);

	printf("after self-advection: \n"); Print_Velocity(vx, vy, vz, mask);

	Bound_Solid();

	Advect_Particles(particles, vx, vy, vz, mask);
	Mark_Water_By(particles, mask);

	//Project(vx, vy, vz, p, p0);
	//Calc_Divergence(vx, vy, vz, s);
	//LOGM("velocity: %f\n", Interpolation_Water_Velocity(2, vz, 2.5, 30.5, 30.5, mask));




	//Calc_Divergence(vx, vy, vz, p);
}
