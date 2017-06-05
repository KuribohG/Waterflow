#include "simulation_cubic.h"
#include<fstream>
#include<sstream>

SimulationCubic::SimulationCubic(){}

void SimulationCubic::Init(int GRIDX, int GRIDY, int GRIDZ){
	printf("cubic init: %d %d %d\n", GRIDX, GRIDY, GRIDZ);
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
            }
        }
    }
}

SimulationCubic::~SimulationCubic() {
	temp_dis.clear();
	vx.clear(), vx0.clear(), vy.clear(), vy0.clear(), vz.clear(), vz0.clear();
	p.clear(), p0.clear();
	mask.clear();
}

void SimulationCubic::Mark_Single_Water(int i, int j, int k) {
	if (mask.is(i, j, k, AIR)) {
		mask(i, j, k) = WATER;
	}
}

void SimulationCubic::Mark_Water(int x0, int x1, int y0, int y1, int z0, int z1) {
	for (int i = x0; i <= x1; i++) {
		for (int j = y0; j <= y1; j++) {
			for (int k = z0; k <= z1; k++) {
				Mark_Single_Water(i, j, k);
			}
		}
	}
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
	printf("apply external forces\n");
	for (int i = 0; i < vz.n; i++) {
		for (int j = 0; j < vz.m; j++) {
			for (int k = 0; k < vz.w; k++) {
				vz(i, j, k) += -g*TIME_DELTA;
			}
		}
	}
}

void SimulationCubic::Cancel_Air_Velocity(void) {
	for (int i = 0; i < GRIDX; i++) {
		for (int j = 0; j < GRIDY; j++) {
			for (int k = 0; k < GRIDZ; k++) {
				if (mask(i, j, k) != WATER && i - 1 >= 0 && mask(i - 1, j, k) != WATER) {
					vx(i, j, k) = 0;
				}
				if (mask(i, j, k) != WATER && j - 1 >= 0 && mask(i, j - 1, k) != WATER) {
					vy(i, j, k) = 0;
				}
				if (mask(i, j, k) != WATER && k - 1 >= 0 && mask(i, j, k - 1) != WATER) {
					vz(i, j, k) = 0;
				}
			}
		}
	}
}

void SimulationCubic::Calc_Divergence(aryf &vx, aryf &vy, aryf &vz, aryf &div) {
	printf("calc divergence:\n");
	div.set(0);
	for (int i = 0; i < GRIDX; i++) {
		for (int j = 0; j < GRIDY; j++) {
			for (int k = 0; k < GRIDZ; k++) {
				if (mask(i, j, k) == WATER) {
					div(i, j, k) =
						+vx(i + 1, j, k) - vx(i, j, k)
						+ vy(i, j + 1, k) - vy(i, j, k)
						+ vz(i, j, k + 1) - vz(i, j, k);
					printf("water: %d %d %d, vz: %f %f, div: %.15f\n", i, j, k, vz.get(i, j, k), vz.get(i, j, k + 1), div.get(i, j, k));
					if (mask.is(i - 1, j, k, SOLID)) div(i, j, k) -= vx.get(i, j, k);
					if (mask.is(i + 1, j, k, SOLID)) div(i, j, k) += vx.get(i + 1, j, k);
					if (mask.is(i, j - 1, k, SOLID)) div(i, j, k) -= vy.get(i, j, k);
					if (mask.is(i, j + 1, k, SOLID)) div(i, j, k) += vy.get(i, j + 1, k);
					if (mask.is(i, j, k - 1, SOLID)) div(i, j, k) -= vz.get(i, j, k);
					if (mask.is(i, j, k + 1, SOLID)) div(i, j, k) += vz.get(i, j, k + 1);
				}
			}
		}
	}
}

void SimulationCubic::Project(aryf &vx,aryf &vy,aryf &vz,aryf &p,aryf &div) {
	LOGM("project\n");
	//Cancel_Air_Velocity();
	//Print_Velocity(vx, vy, vz, mask);
	//Calc_Divergence(vx, vy, vz, p0);
	solver.Solve_Pressure(vx, vy, vz, mask);
	solver.Send_Back_To(p);
	//printf("pressure density:\n");Print_Density(p);
	//return;
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
				if (i > 0 && !s0 && !mask.is(i - 1, j, k, SOLID)) {
					vx(i, j, k) -= (p(i, j, k) - p(i - 1, j, k))*scale;
				}
				else {
					//if (vx.get(i, j, k) != 0) printf("X solid boundary: %d %d %d %f\n", i, j, k, vx.get(i, j, k));
					vx(i, j, k) = 0;
				}
				if (j > 0 && !s0 && !mask.is(i, j - 1, k, SOLID)) vy(i, j, k) -= (p(i, j, k) - p(i, j - 1, k))*scale;
				else {
					//if (vy.get(i, j, k) != 0) printf("Y solid boundary: %d %d %d %f\n", i, j, k, vy.get(i, j, k));
					vy(i, j, k) = 0;
				}
				if (k > 0 && !s0 && !mask.is(i, j, k - 1, SOLID)) {
					vz(i, j, k) -= (p(i, j, k) - p(i, j, k - 1))*scale;
					if (p(i, j, k) - p(i, j, k - 1) != 0) {
						//printf("update vz: %d %d %d %f\n", i, j, k, (p(i, j, k) - p(i, j, k - 1))*scale);
					}
				}
				else {
					//printf("%d %d %d,s0: %d, k-1: %d\n", i, j, k, s0, mask.is(i, j, k - 1, SOLID));
					//if (vz.get(i, j, k) != 0) printf("Z solid boundary: %d %d %d %f\n", i, j, k, vz.get(i, j, k));
					vz(i, j, k) = 0;
				}
			}
		}
	}
}

void Bound_Add(Float &x, Float a, Float bound_low, Float bound_high) {
	if (bound_low <= x + a&&x + a <= bound_high) x += a;
}

void SimulationCubic::Runge_Kutta(int axis, int i, int j, int k, Float delta, int iter,
	const aryf &vx, const aryf &vy, const aryf &vz, Float &x, Float &y, Float &z) {
	assert(iter >= 1);

	delta /= iter;
	x = i, y = j, z = k;
	if (axis == _X) y += 0.5, z += 0.5;
	else if (axis == _Y) x += 0.5, z += 0.5;
	else if (axis == _Z) x += 0.5, y += 0.5;
	//int id = ID(i, j, k);
	Bound_Add(x, -vx.get(i, j, k), 0, GRIDX);
	Bound_Add(y, -vy.get(i, j, k), 0, GRIDY);
	Bound_Add(z, -vz.get(i, j, k), 0, GRIDZ);
	/*x -= vx.get(i, j, k);
	y -= vy.get(i, j, k);
	z -= vz.get(i, j, k);*/
	/*x = i + 0.5 - delta * vx.get(i, j, k);
	y = j + 0.5 - delta * vy.get(i, j, k);
	z = k + 0.5 - delta * vz.get(i, j, k);*/
	//NOTE: x,y,z should be bounded here
	for (int _ = 1; _ < iter; _++) {
		Bound_Add(x, -delta*Interpolation_Water_Velocity(_X, vx, x, y, z, mask), 0, GRIDX);
		Bound_Add(y, -delta*Interpolation_Water_Velocity(_Y, vy, x, y, z, mask), 0, GRIDY);
		Bound_Add(z, -delta*Interpolation_Water_Velocity(_Z, vz, x, y, z, mask), 0, GRIDZ);
		/*Float nx = x - delta*Interpolation_Water_Velocity(_X, vx, x, y, z, mask);
		Float ny = y - delta*Interpolation_Water_Velocity(_Y, vy, x, y, z, mask);
		Float nz = z - delta*Interpolation_Water_Velocity(_Z, vz, x, y, z, mask);
		x = nx, y = ny, z = nz;*/
	}
}

void SimulationCubic::Advect_Velocity(int axis, aryf &f, const aryf &f0, const aryf &vx, const aryf &vy, const aryf &vz) {
	LOGM("advect along: %d\n", axis);
    for (int i = 0; i < f.n; i++) {
        for (int j = 0; j < f.m; j++) {
            for (int k = 0; k < f.w; k++) {
				f(i, j, k) = f0.get(i, j, k);
				bool flag = false;
				if (axis == _X) {
					if (mask.is(i, j, k, WATER) || mask.is(i - 1, j, k, WATER)) flag = true;
				}
				else if (axis == _Y) {
					if (mask.is(i, j, k, WATER) || mask.is(i, j - 1, k, WATER)) flag = true;
				}
				else if (axis == _Z) {
					if (mask.is(i, j, k, WATER) || mask.is(i, j, k - 1, WATER)) flag = true;
				}
				//flag = false;
				if (flag) {
					Float x, y, z;
					Runge_Kutta(axis, i, j, k, TIME_DELTA, 2, vx, vy, vz, x, y, z);
					//LOGM("%f %f %f\n", TIME_DELTA*vx[ID(i, j, k)], TIME_DELTA*vy[ID(i, j, k)], TIME_DELTA*vz[ID(i, j, k)]);
					x = Clip(x, 0.5, GRIDX + 0.5);
					y = Clip(y, 0.5, GRIDY + 0.5);
					z = Clip(z, 0.5, GRIDZ + 0.5);
					f(i, j, k) = Interpolation_Water_Velocity(axis, f0, x, y, z, mask);
				}
				/*int i1 = i, j1 = j, k1 = k;
				if (axis == _X) {
					if (mask.is(i, j, k, WATER) || mask.is(i + 1, j, k, WATER)) {
						Float x, y, z;
					}
				}
				if (mask.is(i, j, k, WATER) || mask.is(i + 1, j, k, WATER)) {
					Float x, y, z;

				}
                if (mask(i, j, k) == WATER) {
                    Float x, y, z;
					Runge_Kutta(axis, i, j, k, TIME_DELTA, 2, vx, vy, vz, x, y, z);
                    //LOGM("%f %f %f\n", TIME_DELTA*vx[ID(i, j, k)], TIME_DELTA*vy[ID(i, j, k)], TIME_DELTA*vz[ID(i, j, k)]);
                    x = Clip(x, 0.5, GRIDX + 0.5);
                    y = Clip(y, 0.5, GRIDY + 0.5);
                    z = Clip(z, 0.5, GRIDZ + 0.5);
					f(i, j, k) = Interpolation_Water_Velocity(axis, f0, x, y, z, mask);
                }*/
            }
        }
    }
	//Bound_Solid(axis, f);
}

Float Kernel(Float dx, Float dy, Float dz) {
	Float res = 1;
	if (dx >= 0 && dx <= 1) {
		res *= 1 - dx;
	} else if (dx >= -1 && dx <= 0) {
		res *= 1 + dx;
	} else return 0;
	if (dy >= 0 && dy <= 1) {
		res *= 1 - dy;
	} else if (dy >= -1 && dy <= 0) {
		res *= 1 + dy;
	} else return 0;
	if (dz >= 0 && dz <= 1) {
		res *= 1 - dz;
	} else if (dz >= -1 && dz <= 0) {
		res *= 1 + dz;
	} else return 0;
	return res;
}

void SimulationCubic::Advect_PIC_Preprocess() {
	printf("advect PIC preprocess\n");
	for (int i = 0; i < GRIDX; i++) {
		for (int j = 0; j < GRIDY; j++) {
			for (int k = 0; k < GRIDZ; k++) {
				v[i][j][k].clear();
			}
		}
	}
	for (MarkerParticle &p : particles) {
		int x = (int)std::floor(p.x);
		int y = (int)std::floor(p.y);
		int z = (int)std::floor(p.z);
		if (x < 0 || x >= GRIDX) continue;
		if (y < 0 || y >= GRIDY) continue;
		if (z < 0 || z >= GRIDZ) continue;
		v[x][y][z].emplace_back(&p);
	}
}

void SimulationCubic::Advect_PIC(int axis, aryf &f, const aryf &f0, const aryf &vx, const aryf &vy, const aryf &vz) {
	printf("advect PIC along axis %d\n", axis);
	if (axis == _X) {
		for (int i = 0; i <= GRIDX; i++) {
			for (int j = 0; j < GRIDY; j++) {
				for (int k = 0; k < GRIDZ; k++) {
					Float u = 0, d = 0;
					for (int dx = -1; dx <= 0; dx++) {
						for (int dy = -1; dy <= 1; dy++) {
							for (int dz = -1; dz <= 1; dz++) {
								if (i + dx < 0 || i + dx >= GRIDX) continue;
								if (j + dy < 0 || j + dy >= GRIDY) continue;
								if (k + dz < 0 || k + dz >= GRIDZ) continue;
								for (MarkerParticle *ptr : v[i + dx][j + dy][k + dz]) {
									MarkerParticle &p = *ptr;
									Float ker = Kernel(i - p.x, j + 0.5 - p.y, k + 0.5 - p.z);
									d += ker, u += ker * p.vx;
								}
							}
						}
					}
					if (fabs(d) < EPS) f(i, j, k) = 0;
					else f(i, j, k) = u / d;
				}
			}
		}
	}
	if (axis == _Y) {
		for (int i = 0; i < GRIDX; i++) {
			for (int j = 0; j <= GRIDY; j++) {
				for (int k = 0; k < GRIDZ; k++) {
					Float u = 0, d = 0;
					for (int dx = -1; dx <= 1; dx++) {
						for (int dy = -1; dy <= 0; dy++) {
							for (int dz = -1; dz <= 1; dz++) {
								if (i + dx < 0 || i + dx >= GRIDX) continue;
								if (j + dy < 0 || j + dy >= GRIDY) continue;
								if (k + dz < 0 || k + dz >= GRIDZ) continue;
								for (MarkerParticle *ptr : v[i + dx][j + dy][k + dz]) {
									MarkerParticle &p = *ptr;
									Float ker = Kernel(i + 0.5 - p.x, j - p.y, k + 0.5 - p.z);
									d += ker, u += ker * p.vy;
								}
							}
						}
					}
					if (fabs(d) < EPS) f(i, j, k) = 0;
					else f(i, j, k) = u / d;
				}
			}
		}
	}
	if (axis == _Z) {
		for (int i = 0; i < GRIDX; i++) {
			for (int j = 0; j < GRIDY; j++) {
				for (int k = 0; k <= GRIDZ; k++) {
					Float u = 0, d = 0;
					for (int dx = -1; dx <= 1; dx++) {
						for (int dy = -1; dy <= 1; dy++) {
							for (int dz = -1; dz <= 0; dz++) {
								if (i + dx < 0 || i + dx >= GRIDX) continue;
								if (j + dy < 0 || j + dy >= GRIDY) continue;
								if (k + dz < 0 || k + dz >= GRIDZ) continue;
								for (MarkerParticle *ptr : v[i + dx][j + dy][k + dz]) {
									MarkerParticle &p = *ptr;
									Float ker = Kernel(i + 0.5 - p.x, j + 0.5 - p.y, k - p.z);
									d += ker, u += ker * p.vz;
								}
							}
						}
					}
					if (fabs(d) < EPS) f(i, j, k) = 0;
					else f(i, j, k) = u / d;
				}
			}
		}
	}
}

void swap(aryf &a, aryf &b) {
	aryf c;
	c = a; a = b; b = c;
}


void SimulationCubic::Pour_Source(int framenum, vector<WaterSource> &sources) {
	printf("pour sources\n");
	for (WaterSource &s : sources) {
		if (framenum > s.pourend) continue;
		for (int i = s.x0; i <= s.x1; i++) {
			for (int j = s.y0; j <= s.y1; j++) {
				for (int k = s.z0; k <= s.z1; k++) {
					vx.set(i, j, k, s.init_vx), vx.set(i + 1, j, k, s.init_vx);
					vy.set(i, j, k, s.init_vy), vy.set(i, j + 1, k, s.init_vy);
					vz.set(i, j, k, s.init_vz), vz.set(i, j, k + 1, s.init_vz);
					if (mask.is(i, j, k, AIR)) {
						Mark_Single_Water(i, j, k);
						//Add_Single_Particle(particles, mask, i, j, k, s.init_vx, s.init_vy, s.init_vz);
						Add_Single_Particle(particles, mask, i, j, k);
						//printf("%f %f %f\n", s.init_vx, s.init_vy, s.init_vz);
					}
				}
			}
		}
	}
}

void SimulationCubic::Step_Time(int framenum, vector<WaterSource> &sources){
    printf("cubic step time \n");
	Float t0 = omp_get_wtime(), t;
	static int tot = 0;
    //velocity-evolution
	Apply_External_Forces();
	Pour_Source(framenum, sources);
	printf("particle num: %d\n", particles.size());
	t = omp_get_wtime(); printf("apply external forces&pour source time cost: %.2fs\n", (t - t0 + 0.0)); t0 = t;

	Get_Particles_Velocity(particles, vx, vy, vz, mask);
	Advect_Particles(particles, vx, vy, vz, mask);
	Advect_PIC_Preprocess();
	t = omp_get_wtime(); printf("advect particles time cost: %.2fs\n", (t - t0 + 0.0)); t0 = t;
	
    swap(vx, vx0);
	swap(vy, vy0);
	swap(vz, vz0);
	Advect_PIC(0, vx, vx0, vx0, vy0, vz0);
	Advect_PIC(1, vy, vy0, vx0, vy0, vz0);
	Advect_PIC(2, vz, vz0, vx0, vy0, vz0);

	t = omp_get_wtime(); printf("PIC advection time cost: %.2fs\n", (t - t0 + 0.0)); t0 = t;

	Get_Particles_Velocity(particles, vx, vy, vz, mask);
	Mark_Water_By(particles, mask);

	t = omp_get_wtime(); printf("update water mask time cost: %.2fs\n", (t - t0 + 0.0)); t0 = t;
	
	Project(vx, vy, vz, p, p0);

	t = omp_get_wtime(); printf("project time cost: %.2fs\n", (t - t0 + 0.0)); t0 = t;
}
