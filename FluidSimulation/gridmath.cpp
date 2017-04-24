#include "gridmath.h"



/*template<typename T>
T& Array3D<T>::operator[](int i){
	return f[i];
}*/

template<typename T>
inline int Array3D<T>::ID(int i, int j, int k)const{
	return i*m*w + j*w + k;
}

template<typename T>
inline bool Array3D<T>::inside(int i, int j, int k)const{
	return 0 <= i&&i < n && 0 <= j&&j < m && 0 <= k&&k < w;
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

Float Interpolation_Water_Velocity(int axis, const aryf &f, Float x, Float y, Float z, const aryi &mask) {
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
	if (Valid_Water(x0, y0, z0, mask)) v += s1*u1*v1*f.get(x0, y0, z0), sw += s1*u1*v1;
	if (Valid_Water(x0, y0, z0, mask)) v += s1*u1*v0*f.get(x0, y0, z1), sw += s1*u1*v0;
	if (Valid_Water(x0, y0, z0, mask)) v += s1*u0*v1*f.get(x0, y1, z0), sw += s1*u0*v1;
	if (Valid_Water(x0, y0, z0, mask)) v += s1*u0*v0*f.get(x0, y1, z1), sw += s1*u0*v0;
	if (Valid_Water(x0, y0, z0, mask)) v += s0*u1*v1*f.get(x1, y0, z0), sw += s0*u1*v1;
	if (Valid_Water(x0, y0, z0, mask)) v += s0*u1*v0*f.get(x1, y0, z1), sw += s0*u1*v0;
	if (Valid_Water(x0, y0, z0, mask)) v += s0*u0*v1*f.get(x1, y1, z0), sw += s0*u0*v1;
	if (Valid_Water(x0, y0, z0, mask)) v += s0*u0*v0*f.get(x1, y1, z1), sw += s0*u0*v0;
	if (sw == 0) return 0;
	return v / sw;
}

/*Float Interpolation_In_Water_3D(aryf &f, Float x, Float y, Float z, aryi &mask) {
	// todo: somehow re-write all interpolation scenes to take 0.5 offset into consideration
	//x -= 0.5, y -= 0.5, z -= 0.5;
	int x0 = floor(x), x1 = x0 + 1;
	int y0 = floor(y), y1 = y0 + 1;
	int z0 = floor(z), z1 = z0 + 1;
	Float s0 = x - x0, s1 = x1 - x;
	Float u0 = y - y0, u1 = y1 - y;
	Float v0 = z - z0, v1 = z1 - z;
	Float v = 0, sw = 0;
	if (Valid_Water(x0, y0, z0, mask)) v += s1*u1*v1*f(x0, y0, z0), sw += s1*u1*v1;
	if (Valid_Water(x0, y0, z0, mask)) v += s1*u1*v0*f(x0, y0, z1), sw += s1*u1*v0;
	if (Valid_Water(x0, y0, z0, mask)) v += s1*u0*v1*f(x0, y1, z0), sw += s1*u0*v1;
	if (Valid_Water(x0, y0, z0, mask)) v += s1*u0*v0*f(x0, y1, z1), sw += s1*u0*v0;
	if (Valid_Water(x0, y0, z0, mask)) v += s0*u1*v1*f(x1, y0, z0), sw += s0*u1*v1;
	if (Valid_Water(x0, y0, z0, mask)) v += s0*u1*v0*f(x1, y0, z1), sw += s0*u1*v0;
	if (Valid_Water(x0, y0, z0, mask)) v += s0*u0*v1*f(x1, y1, z0), sw += s0*u0*v1;
	if (Valid_Water(x0, y0, z0, mask)) v += s0*u0*v0*f(x1, y1, z1), sw += s0*u0*v0;
	if (sw == 0) return 0;
	return v / sw;
}*/