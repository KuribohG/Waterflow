#include "gridmath.h"
void Print_Velocity(const aryf &vx, const aryf &vy, const aryf &vz, const aryi &mask) {
	printf("print velocity: \n");
	Float vz0 = INF;
	for (int i = GRIDX / 2; i <= GRIDX / 2; i++) {
		for (int j = 0; j < GRIDY; j++) {
			for (int k = 0; k < GRIDZ; k++) {
				if (j == 29) printf("%d %d %f\n", j, k, vz.get(i, j, k));
			}
		}
	}
}


Float Neighbor_Sum6(aryf &x, int i, int j, int k) {
	Float tmp = 0;
	if (i - 1 >= 0) tmp += x(i - 1, j, k);
	if (i + 1 < x.n) tmp += x(i + 1, j, k);
	if (j - 1 >= 0) tmp += x(i, j - 1, k);
	if (j + 1 < x.m) tmp += x(i, j + 1, k);
	if (k - 1 >= 0) tmp += x(i, j, k - 1);
	if (k + 1 < x.w) tmp += x(i, j, k + 1);
	return tmp;
}

bool Valid_Water(int x, int y, int z, const aryi &mask) {
	return mask.inside(x, y, z) && mask.get(x, y, z) == WATER;
}

Float Extrapolation_Water_Velocity(const aryf &f, Float x, Float y, Float z, const aryi &mask) {
	int x0 = floor(x), y0 = floor(y), z0 = floor(z);
	int high = 7;//magic number!
	Float mnsqr = INF, v;
	for (int s = 0; s <= high; s++) {
		for (int di = -s; di <= s; di++) {
			int r = s - abs(di);
			for (int dj = -r; dj <= r; dj++) {
				for (int sgn = -1; sgn <= 1; sgn += 2) {
					int dk = (r - abs(dj))*sgn;
					if (mask.is(x0 + di, y0 + dj, z0 + dk, WATER)) {
						high = s;
						Float dissqr = (x0 + di - x)*(x0 + di - x) + (y0 + dj - y)*(y0 + dj - y) + (z0 + dk - z)*(z0 + dk - z);
						if (dissqr < mnsqr) {
							mnsqr = dissqr;
							v = f.get(x0 + di, y0 + dj, z0 + dk);
						}
					}
				}
			}
		}
	}
	return v;
}

Float Interpolation_Water_Velocity(int axis, const aryf &f, Float x, Float y, Float z, const aryi &mask, bool extrapolation) {
	if (axis == 0) y -= 0.5, z -= 0.5;
	else if (axis == 1) x -= 0.5, z -= 0.5;
	else if (axis == 2) x -= 0.5, y -= 0.5;
	int x0 = floor(x), x1 = x0 + 1;
	int y0 = floor(y), y1 = y0 + 1;
	int z0 = floor(z), z1 = z0 + 1;
	Float s0 = x - x0, s1 = x1 - x;
	Float u0 = y - y0, u1 = y1 - y;
	Float v0 = z - z0, v1 = z1 - z;
	Float v = 0, sw = 0;
	if (f.inside(x0, y0, z0)) v += s1*u1*v1*f.get(x0, y0, z0), sw += s1*u1*v1;
	if (f.inside(x0, y0, z1))v += s1*u1*v0*f.get(x0, y0, z1), sw += s1*u1*v0;
	if (f.inside(x0, y1, z0))v += s1*u0*v1*f.get(x0, y1, z0), sw += s1*u0*v1;
	if (f.inside(x0, y1, z1))v += s1*u0*v0*f.get(x0, y1, z1), sw += s1*u0*v0;
	if (f.inside(x1, y0, z0))v += s0*u1*v1*f.get(x1, y0, z0), sw += s0*u1*v1;
	if (f.inside(x1, y0, z1))v += s0*u1*v0*f.get(x1, y0, z1), sw += s0*u1*v0;
	if (f.inside(x1, y1, z0))v += s0*u0*v1*f.get(x1, y1, z0), sw += s0*u0*v1;
	if (f.inside(x1, y1, z1))v += s0*u0*v0*f.get(x1, y1, z1), sw += s0*u0*v0;
	/*
	if (Valid_Water(x0, y0, z0, mask)) v += s1*u1*v1*f.get(x0, y0, z0), sw += s1*u1*v1;
	if (Valid_Water(x0, y0, z1, mask)) v += s1*u1*v0*f.get(x0, y0, z1), sw += s1*u1*v0;
	if (Valid_Water(x0, y1, z0, mask)) v += s1*u0*v1*f.get(x0, y1, z0), sw += s1*u0*v1;
	if (Valid_Water(x0, y1, z1, mask)) v += s1*u0*v0*f.get(x0, y1, z1), sw += s1*u0*v0;
	if (Valid_Water(x1, y0, z0, mask)) v += s0*u1*v1*f.get(x1, y0, z0), sw += s0*u1*v1;
	if (Valid_Water(x1, y0, z1, mask)) v += s0*u1*v0*f.get(x1, y0, z1), sw += s0*u1*v0;
	if (Valid_Water(x1, y1, z0, mask)) v += s0*u0*v1*f.get(x1, y1, z0), sw += s0*u0*v1;
	if (Valid_Water(x1, y1, z1, mask)) v += s0*u0*v0*f.get(x1, y1, z1), sw += s0*u0*v0;
	//if(abs()) printf("%d %f %f %f %f %f\n", axis, x, y, z, v, sw);
	if (sw == 0) {
		if (!extrapolation) return 0;
		else return Extrapolation_Water_Velocity(f, x, y, z, mask);
	}
	*/
	//if (fabs(sw) < 1e-10) printf("++++++++++++++++++++++++++++++%f %f %f\n", x, y, z);
	//assert(fabs(sw) > 1e-10);
	assert(sw != 0);
	return v / sw;
}
