#pragma once
#include "shared.hpp"

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
		printf("%d %d %d %p\n", n, m, w, f);
		memset(f, 0, sizeof(T)*n*m*w);
		//printf("%p\n", f);
	}
	T& operator () (int i, int j, int k);
	//T& operator [] (int i);
	inline int ID(int i, int j, int k);
	inline bool inside(int i, int j, int k);
	void copy(Array3D<T> &b) {
		assert(n == b.n);
		assert(m == b.m);
		assert(w == b.w);
		memcpy(f, b.f, sizeof(T)*n*m*w);
	}
	void set(int a) {
		memset(f, a, sizeof(T)*n*m*w);
	}
};
typedef Array3D<Float> aryf;
typedef Array3D<int> aryi;

Float Neighbor_Sum6(aryf &x, int i, int j, int k);
bool Valid_Water(int x, int y, int z, aryi &mask);

Float Interpolation_In_Water_3D(aryf &f, Float x, Float y, Float z, aryi &mask);