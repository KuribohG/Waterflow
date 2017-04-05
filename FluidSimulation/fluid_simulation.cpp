#include "fluid_simulation.h"


FluidSimulation::FluidSimulation()
{
    init_Density_3d();
}
#define D3
void FluidSimulation::Draw_On_Screen(void){
#ifdef D3
    Draw_Density_3d(cubic.density, RIGHT_SCREEN, lightPath);
#else
    Draw_Density_2d(cubic.density, RIGHT_SCREEN);
#endif
    Draw_Velocity_2d(2, cubic.vx, cubic.vy, cubic.vz, LEFT_SCREEN);
}

void FluidSimulation::Step_Time(void){
    cubic.Step_Time();

}

template<class T> struct P_3d
{
    T x, y, z;
    P_3d(T _x, T _y, T _z) : x(_x), y(_y), z(_z) {}
    template<class T2> P_3d(P_3d<T2> u) : x(u.x), y(u.y), z(u.z) {}
    P_3d operator*(const T a) const { return P_3d(x * a, y * a, z * a);  }
    P_3d operator+(const P_3d& o) { return P_3d(x + o.x, y + o.y, z + o.z); }
    P_3d operator-(const P_3d& o) { return P_3d(x - o.x, y - o.y, z - o.z); }
    T abs() const { return sqrt(x * x + y * y + z * z); }
    P_3d cross(P_3d o) const { return P_3d(y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x); }
    tuple<int, int, int> to_int()  { return make_tuple(int(x), int(y), int(z)); };
};
template<class T> struct Cuboid
{
    T xm, xM, ym, yM, zm, zM;
    Cuboid (P_3d<T> p1, P_3d<T> p2) :
            xm(min(p1.x, p2.x)), xM(max(p1.x, p2.x)),
            ym(min(p1.y, p2.y)), yM(max(p1.y, p2.y)),
            zm(min(p1.z, p2.z)), zM(max(p1.z, p2.z)) { }
};
template<class T> struct Segment
{
    P_3d<T> p1, p2;
    Segment(P_3d<T> _p1, P_3d<T> _p2) : p1(_p1), p2(_p2) {}
    T len() { return (p1 - p2).abs(); }
    bool Clip(Cuboid<T> View)
    {
        T umin = 0, umax = 1,
                x1 = p1.x, dx = p2.x - x1,
                y1 = p1.y, dy = p2.y - y1,
                z1 = p1.z, dz = p2.z - z1 ;
#define CL(l) if (fabs(d##l) > EPS<T>) {\
         T a =  (View.l##M - l##1) / d##l, b =  (View.l##m - l##1) / d##l;\
         if (a > b) swap(a, b); umin = max(umin, a); umax = min(umax, b); }

        CL(x); CL(y); CL(z);
#undef CL
        if (umin >= umax - EPS<T>) return false;
        p1 = P_3d<T>(x1 + dx * umin, y1 + dy * umin, z1 + dz * umin);
        p2 = P_3d<T>(x1 + dx * umax, y1 + dy * umax, z1 + dz * umax);
        return true;
    }
    P_3d<T> atx(T x) throw(int)
    {
        if (fabs(p2.x - p1.x) < EPS<T>)
            throw 233;
        T inc = (x - p1.x) / (p2.x - p1.x);
        return P_3d<T>(x, p1.y + (p2.y - p1.y) * inc, p1.z + (p2.z - p1.z) * inc);

    }
};
const Float pixel_size = 1.0f / 8.0f;
using p3df = P_3d<Float>;
using sef = Segment<Float>;
using cuf = Cuboid<Float>;
void FluidSimulation::init_Density_3d()
{
    p3df vrp(30, 32, 0), vpn(-1, 0, 0), uvp(0, 0, 1), vvp(vpn.cross(uvp)), prp(50, 32, 32);
    cuf jar(p3df(0, 0, 0), p3df(GRID_SIZE_X, GRID_SIZE_Y, GRID_SIZE_Z));
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
                    thu.p1 = thu.p1 + (thu.p2 - thu.p1) * (0.1 / thu.len());

                    if (thu.p1.to_int() == thu.p2.to_int()) break;
                    int x, y, z; tie(x, y, z) = thu.p1.to_int();

                    thu.Clip(cuf(p3df(x, y, z), p3df(x + 1, y + 1, z + 1)));

                    lightPath[i][j].push_back(make_pair(ID(x, y, z), thu.len()));
                    ray.p2 = thu.p1;
                }
            } catch(int err) { fprintf(stderr, "err: %d\n", err); }
}

