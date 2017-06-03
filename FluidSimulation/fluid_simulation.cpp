#include "fluid_simulation.h"
#include <cmath>
#include <queue>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <string>
#include <cstdlib>
using namespace std;

Float Random(void) {
	return rand() / (RAND_MAX - 1.0);
}

WaterSource::WaterSource(int _x0, int _x1, int _y0, int _y1, int _z0, int _z1, Float _gen_rate, Float _init_vx, Float _init_vy, Float _init_vz) :
	x0(_x0), x1(_x1), y0(_y0), y1(_y1), z0(_z0), z1(_z1), gen_rate(_gen_rate / FPS), init_vx(_init_vx), init_vy(_init_vy), init_vz(_init_vz)
{
}

FluidSimulation::FluidSimulation():cubic(){
	string filename;
	cout << "please enter scene file name: ";
	//cin >> filename;
	filename = "scenes/test.box";
	cout << filename << endl;
	Read_Scene_File(filename.c_str());
	//mask(2, 30, 30) = WATER;
	LOGM("velocity set\n");
	Init_Particles(cubic.particles, cubic.mask);
    signed_dis.init(GRIDX, GRIDY, GRIDZ);
    init_Density_3d();
}

#ifdef OPENGL

//#define D3
void FluidSimulation::Draw_On_Screen(void){
	printf("draw on screen: \n");
#ifdef D3
    Draw_Density_3d(cubic.density, RIGHT_SCREEN, lightPath);
#else
	//Draw_Nearest(GRIDX / 2, nearest);
	//Draw_Mask_2d(cubic.mask, THIRD_SCREEN);
    Draw_Density_2d(cubic.p, THIRD_SCREEN);
#endif
	Draw_Velocity_2d(cubic.vx, cubic.vy, cubic.vz, cubic.mask, LEFT_SCREEN);
	Draw_Particle_2d(cubic.particles, RIGHT_SCREEN);
}

#endif

void FluidSimulation::Read_Scene_File(const char * filename) {
	cerr << "read scene file: " << filename << endl;
	ifstream fin;
	fin.open(filename);
	string buff, cmd;
	while (getline(fin, buff)) {
		stringstream sin(buff);
		sin >> cmd;
		//cout << buff << endl;
		if (buff[0] == '#');
		else if (cmd == "box") {
			int x0, x1, y0, y1, z0, z1;
			sin >> x0 >> x1 >> y0 >> y1 >> z0 >> z1;
			cubic.Mark_Water(x0, x1, y0, y1, z0, z1);
		}
		else if (cmd == "source") {
			int x0, x1, y0, y1, z0, z1;
			Float rate, vx, vy, vz;
			sin >> x0 >> x1 >> y0 >> y1 >> z0 >> z1 >> rate >> vx >> vy >> vz;
			sources.emplace_back(x0, x1, y0, y1, z0, z1, rate, vx, vy, vz);
		}
	}
	fin.close();
}

void FluidSimulation::Pour_Source(void) {
	for (WaterSource &s : sources) {
		for (int i = s.x0; i <= s.x1; i++) {
			for (int j = s.y0; j <= s.y1; j++) {
				for (int k = s.z0; k <= s.z1; k++) {
					Float t = Random();
					if (t <= s.gen_rate) {
						cubic.Mark_Single_Water(i, j, k);
						Add_Single_Particle(cubic.particles, cubic.mask, i, j, k);
					}
				}
			}
		}
	}
}

Float Kernel_Func(Float s_square) {
    if (s_square >= 1) {
        return 0;
    }
    Float x = 1 - s_square;
    return x * x * x;
}



Float Square_Dis(Float x, Float y, Float z, Float xx, Float yy, Float zz) {
    return (x - xx) * (x - xx) + (y - yy) * (y - yy) + (z - zz) * (z - zz);
}

void FluidSimulation::Calculate_Signed_Distance() {
    const Float h = 4.0;
    const Float r = 1.0;
    for (int i = 0; i < GRIDX; i++) {
        for (int j = 0; j < GRIDY; j++) {
            for (int k = 0; k < GRIDZ; k++) {
                cubic.v[i][j][k].clear();
            }
        }
    }
    for (MarkerParticle &p : cubic.particles) {
        int x = (int)std::floor(p.x);
        int y = (int)std::floor(p.y);
        int z = (int)std::floor(p.z);
		if (x < 0 || x >= GRIDX) continue;
		if (y < 0 || y >= GRIDY) continue;
		if (z < 0 || z >= GRIDZ) continue;
        cubic.v[x][y][z].emplace_back(&p);
    }
    for (int i = 0; i < GRIDX; i++) {
        for (int j = 0; j < GRIDY; j++) {
            for (int k = 0; k < GRIDZ; k++) {
                Float X = 0, Y = 0, Z = 0, d = 0;
                Float x = i + 0.5, y = j + 0.5, z = k + 0.5;
                for (int dx = -2; dx <= 2; dx++) {
                    for (int dy = -2; dy <= 2; dy++) {
                        for (int dz = -2; dz <= 2; dz++) {
                            if (std::abs(dx) + std::abs(dy) + std::abs(dz) >= 5) continue;
                            if (i + dx < 0 || i + dx >= GRIDX) continue;
                            if (j + dy < 0 || j + dy >= GRIDY) continue;
                            if (k + dz < 0 || k + dz >= GRIDZ) continue;
                            for (MarkerParticle *ptr : cubic.v[i + dx][j + dy][k + dz]) {
                                MarkerParticle &p = *ptr;
                                Float dis_square = Square_Dis(p.x, p.y, p.z, x, y, z);
                                Float weight = Kernel_Func(dis_square / h);
                                d += weight, X += weight * p.x, Y += weight * p.y, Z += weight * p.z;
                            }
                        }
                    }
                }
                X /= d, Y /= d, Z /= d;
                Float norm = sqrt(Square_Dis(x, y, z, X, Y, Z));
                signed_dis(i, j, k) = norm - r;
            }
        }
    }
	for (int i = 0; i < GRIDX; i++) {
		for (int j = 0; j < GRIDY; j++) {
			for (int k = 0; k < GRIDZ; k++) {
				Float w = signed_dis.get(i, j, k);
				if (!finite(w) || isnan(w)) {
					signed_dis(i, j, k) = 1;
				}
                // hack on signed_dis, to make mesh close
                if (cubic.mask(i, j, k) == SOLID && signed_dis(i, j, k) < 1e-5) {
                    signed_dis(i, j, k) = 0.0001;
                } else if (cubic.mask(i, j, k) == WATER && signed_dis(i, j, k) > -1e-5) {
                    signed_dis(i, j, k) = -10;
                }
				//else if (fabs(w) < TSDF_EPS) signed_dis(i, j, k) = 0;
			}
		}
	}
}

struct Position {
    int x, y, z;
    Position() {}
    Position(int _x, int _y, int _z)
            : x(_x), y(_y), z(_z) {}
};

std::queue<Position> known;

void FluidSimulation::Calculate_Nearest_Particle() {
    const Float INF = 1e20;
    for (int i = 0; i < GRIDX; i++) {
        for (int j = 0; j < GRIDY; j++) {
            for(int k = 0; k < GRIDZ; k++) {
                dis[i][j][k] = INF;
                nearest[i][j][k] = nullptr;
            }
        }
    }
    for (int i = 0; i < GRIDX; i++) {
        for (int j = 0; j < GRIDY; j++) {
            for (int k = 0; k < GRIDZ; k++) {
                if (cubic.mask(i, j, k) == WATER) continue;
                for (int dx = -1; dx <= 1; dx++) {
                    for (int dy = -1; dy <= 1; dy++) {
                        for (int dz = -1; dz <= 1; dz++) {
                            if (i + dx < 0 || i + dx >= GRIDX) continue;
                            if (j + dy < 0 || j + dy >= GRIDY) continue;
                            if (k + dz < 0 || k + dz >= GRIDZ) continue;
                            vector<MarkerParticle *> &particles = cubic.v[i + dx][j + dy][k + dz];
                            Float X = i + 0.5, Y = j + 0.5, Z = k + 0.5;
                            for (MarkerParticle *particle : particles) {
                                Float d = Square_Dis(particle->x, particle->y, particle->z, X, Y, Z);
                                if (!nearest[i][j][k] || d < dis[i][j][k]) {
                                    dis[i][j][k] = d;
                                    nearest[i][j][k] = particle;
                                    known.push(Position(i, j, k));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    while (!known.empty()) {
        Position pos = known.front();
        known.pop();
        int i = pos.x, j = pos.y, k = pos.z;
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                for (int dz = -1; dz <= 1; dz++) {
                    if (i + dx < 0 || i + dx >= GRIDX) continue;
                    if (j + dy < 0 || j + dy >= GRIDY) continue;
                    if (k + dz < 0 || k + dz >= GRIDZ) continue;
                    if (cubic.mask(i + dx, j + dy, k + dz) == WATER) continue;
                    Float x = nearest[i][j][k]->x, y = nearest[i][j][k]->y, z = nearest[i][j][k]->z;
                    Float d = Square_Dis(x, y, z, i + dx + 0.5, j + dy + 0.5, k + dz + 0.5);
                    if (!nearest[i + dx][j + dy][k + dz] || d < dis[i + dx][j + dy][k + dz]) {
                        dis[i + dx][j + dy][k + dz] = d;
                        nearest[i + dx][j + dy][k + dz] = nearest[i][j][k];
                        known.push(Position(i + dx, j + dy, k + dz));
                    }
                }
            }
        }
    }
    for (int i = 0; i < GRIDX; i++) {
        for (int j = 0; j < GRIDY; j++) {
            for (int k = 0; k < GRIDZ; k++) {
                if (nearest[i][j][k]) {
                    signed_dis(i, j, k) = dis[i][j][k];
                }
            }
        }
    }
}

void FluidSimulation::Get_Full_Velocity() {
	printf("extrapolating: \n");
    for (int i = 0; i < GRIDX; i++) {
        for (int j = 0; j < GRIDY; j++) {
            for (int k = 0; k < GRIDZ; k++) {
				if (cubic.mask(i, j, k) != WATER && i - 1 >= 0 && cubic.mask(i - 1, j, k) != WATER) {
					Float x = nearest[i][j][k]->vx;
					Float xx = nearest[i - 1][j][k]->vx;
					cubic.vx(i, j, k) = (x + xx) / 2;
				}
                if (cubic.mask(i, j, k) != WATER && j - 1 >= 0 && cubic.mask(i, j - 1, k) != WATER) {
                    Float y = nearest[i][j][k]->vy;
                    Float yy = nearest[i][j - 1][k]->vy;
                    cubic.vy(i, j, k) = (y + yy) / 2;
                }
				if (cubic.mask(i, j, k) != WATER && k - 1 >= 0 && cubic.mask(i, j, k - 1) != WATER) {
					Float z = nearest[i][j][k]->vz;
					Float zz = nearest[i][j][k - 1]->vz;
                    cubic.vz(i, j, k) = (z + zz) / 2;
				}
				//if (i == GRIDX / 2 && j == 2) printf("%d %d %d %f\n", i, j, k, cubic.vz(i, j, k));
            }
        }
    }
}

void FluidSimulation::Step_Time(void){
	printf("==========================================================================================\n");
	static int T0 = clock();
	int tstep = clock();
	framenum++;
	printf("start to step frame %d\n", framenum);
	Pour_Source();
    cubic.Step_Time();
	int t0 = clock();
    Calculate_Signed_Distance();
	int t1 = clock(); printf("calc signed distance time cost: %.2fs\n", (t1 - t0 + 0.0) / CLOCKS_PER_SEC);
    Calculate_Nearest_Particle();
	int t2 = clock(); printf("calc nearest particle time cost: %.2fs\n", (t2 - t1 + 0.0) / CLOCKS_PER_SEC);
	//printf("before extrapolation: \n"); Print_Velocity(cubic.vx, cubic.vy, cubic.vz, cubic.mask);
	Get_Full_Velocity();
	int t3 = clock(); printf("calc nearest particle time cost: %.2fs\n", (t2 - t1 + 0.0) / CLOCKS_PER_SEC);
	//printf("after extrapolation: \n"); Print_Velocity(cubic.vx, cubic.vy, cubic.vz, cubic.mask);
    
	if (framenum % 1 == 0) {

		char name[50];
		sprintf(name, "objs/meshs.%04d.obj", framenum);
		meshcubes.Reconstruct(signed_dis, 0.0);
		printf("dump to: %s\n", name);
		char pngname[100];
		sprintf(pngname, "%04d.png", framenum);
		
		meshcubes.Dump_Obj(name);
		//meshcubes.Dump_GOC(name, pngname, 1200, 900);
		//getchar();
	}
	//if (framenum >= 0) { printf("input: \n"); getchar(); }
	//LOGM("continue\n");
	int t4 = clock();
	printf("frame %d step done, step time cost: %.2lfs, all time cost: %.2lfs\n", framenum, (t4-tstep+0.0)/CLOCKS_PER_SEC, (t4 - T0 + 0.0) / CLOCKS_PER_SEC);
	printf("==========================================================================================\n");
	//if (framenum >= 200) exit(0);
}

struct P_3d
{
    Float x, y, z;
    P_3d(Float _x, Float _y, Float _z) : x(_x), y(_y), z(_z) {}
    P_3d operator*(const Float a) const { return P_3d(x * a, y * a, z * a);  }
    P_3d operator+(const P_3d& o) { return P_3d(x + o.x, y + o.y, z + o.z); }
    P_3d operator-(const P_3d& o) { return P_3d(x - o.x, y - o.y, z - o.z); }
    Float abs() const { return sqrt(x * x + y * y + z * z); }
    P_3d cross(P_3d o) const { return P_3d(y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x); }
    tuple<int, int, int> to_int()  { return make_tuple(int(x), int(y), int(z)); };
};
 struct Cuboid
{
    Float xm, xM, ym, yM, zm, zM;
    Cuboid (P_3d p1, P_3d p2) :
            xm(min(p1.x, p2.x)), xM(max(p1.x, p2.x)),
            ym(min(p1.y, p2.y)), yM(max(p1.y, p2.y)),
            zm(min(p1.z, p2.z)), zM(max(p1.z, p2.z)) { }
};
 struct Segment
{
    P_3d p1, p2;
    Segment(P_3d _p1, P_3d _p2) : p1(_p1), p2(_p2) {}
    Float len() { return (p1 - p2).abs(); }
    bool Clip(Cuboid View)
    {
        Float umin = 0, umax = 1,
                x1 = p1.x, dx = p2.x - x1,
                y1 = p1.y, dy = p2.y - y1,
                z1 = p1.z, dz = p2.z - z1 ;
#define CL(l) if (fabs(d##l) > EPS) {\
         Float a =  (View.l##M - l##1) / d##l, b =  (View.l##m - l##1) / d##l;\
         if (a > b) swap(a, b); umin = max(umin, a); umax = min(umax, b); }

        CL(x); CL(y); CL(z);
#undef CL
        if (umin >= umax - EPS) return false;
        p1 = P_3d(x1 + dx * umin, y1 + dy * umin, z1 + dz * umin);
        p2 = P_3d(x1 + dx * umax, y1 + dy * umax, z1 + dz * umax);
        return true;
    }
    P_3d atx(Float x) throw(int)
    {
        if (fabs(p2.x - p1.x) < EPS)
            throw 233;
        Float inc = (x - p1.x) / (p2.x - p1.x);
        return P_3d(x, p1.y + (p2.y - p1.y) * inc, p1.z + (p2.z - p1.z) * inc);

    }
};
const Float pixel_size = 1.0f / 8.0f;
using p3df = P_3d;
using sef = Segment;
using cuf = Cuboid;
void FluidSimulation::init_Density_3d()
{
  /*double ang = -0;
    p3df vrp(30, 32, 0), vpn(-cos(ang), sin(ang), 0), uvp(0, 0, 1), vvp(vpn.cross(uvp)), prp(50, 32, 32);
    cuf jar(p3df(0, 0, 0), p3df(GRIDX, GRIDY, GRIDZ));
    for (int i = 0; i < SHOW_SIZE_X; ++i)
        for (int j = 0; j < SHOW_SIZE_Y; ++j)
            try
            {
                sef ray(prp, vrp + (vvp * i + uvp * j) * pixel_size);

                ray.p2 = ray.atx(0); //保证光线不能与OXY平面平行
                if (!ray.Clip(jar)) continue;
                while (true)
                {
                    sef thu = ray;
                    if (thu.len() < EPS) break;
                    p3df pPos = thu.p1 + (thu.p2 - thu.p1) * (0.05 / thu.len());

                    if (pPos.to_int() == thu.p2.to_int()) break;
                    int x, y, z; tie(x, y, z) = pPos.to_int();

                    thu.Clip(cuf(p3df(x, y, z), p3df(x + 1, y + 1, z + 1)));

                    lightPath[i][j].push_back(make_pair(ID(x, y, z), thu.len()));
                    ray.p1 = thu.p2;
                }
            } catch(int err) { LOGM("err: %d\n", err); }*/
}

