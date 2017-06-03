#pragma once
#include "shared.h"

template<typename T>
class Array3D {
public:
	T *f = NULL;
	int n, m, w;
	Array3D() {}
	~Array3D(void) {}
	void clear(void) {
		if (f != NULL) free(f);
	}
	void operator = (Array3D<T> &b) {
		//printf("before: %p %p\n", f, b.f);
		f = b.f;
		n = b.n;
		m = b.m;
		w = b.w;
		//printf("after: %p %p\n", f, b.f);
	}
	void init(int _n, int _m, int _w) {
		n = _n;
		m = _m;
		w = _w;
		f = (T*)malloc(n*m*w * sizeof(T));
		//printf("%d %d %d %p\n", n, m, w, f);
		memset(f, 0, sizeof(T)*n*m*w);
		//printf("%p\n", f);
	}
	int ID(int i, int j, int k)const {
		return i*m*w + j*w + k;
	}
	T get(int i, int j, int k)const {
		return f[ID(i, j, k)];
	}
	T& operator () (int i, int j, int k) {
		assert(f != NULL);
		assert(0 <= i&&i < n && 0 <= j&&j < m && 0 <= k&&k < w);
		//assert(fabs(f[i*m*w + j*w + k]) < 1000);
		return f[i*m*w + j*w + k];
	}
	bool inside(int i, int j, int k)const {
		return 0 <= i&&i < n && 0 <= j&&j < m && 0 <= k&&k < w;
	}
	void copy(Array3D<T> &b) {
		assert(n == b.n);
		assert(m == b.m);
		assert(w == b.w);
		memcpy(f, b.f, sizeof(T)*n*m*w);
	}
	bool is(int i, int j, int k, T typ)const {
		return inside(i, j, k) && get(i, j, k) == typ;
	}
	void set(int a) {
		memset(f, a, sizeof(T)*n*m*w);
	}
};
template<typename T>
class Vec1D {
public:
	T *f = NULL;
	int n;
	Vec1D() {}
	~Vec1D(void) {}
	void reset(void) {
		n = 0;
	}
	void clear(void) {
		if (f != NULL) free(f);
	}
	void init(int _n) {
		n = _n;
		f = (T*)malloc(n * sizeof(T));
		memset(f, 0, sizeof(T)*n);
	}
	T get(int i)const {
		return f[i];
	}
	T& operator [] (int i) {
		assert(0 <= i&&i < n);
		return f[i];
	}
};
typedef Array3D<Float> aryf;
typedef Array3D<int> aryi;
typedef Vec1D<double> vecd;

Float Neighbor_Sum6(aryf &x, int i, int j, int k);
bool Valid_Water(int x, int y, int z, aryi &mask);

Float Interpolation_Water_Velocity(int axis, const aryf &f, Float x, Float y, Float z, const aryi &mask, bool extrapolation = false);
//Float Interpolation_In_Water_3D(aryf &f, Float x, Float y, Float z, aryi &mask);

void Print_Velocity(const aryf &vx, const aryf &vy, const aryf &vz, const aryi &mask);

void Print_Density(const aryf &vx);