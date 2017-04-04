#include "drawer.h"


static void ScreenCoor_to_ClipCoor(Float &x, Float &y, SCREENID_T screenid) {
    //x *= (SHOW_SIZE_X / GRID_SIZE_X);
    //y *= (SHOW_SIZE_Y / GRID_SIZE_Y);
    x = (x + screenid*SHOW_SIZE_X) * 2 / (SHOW_SIZE_X*TOTAL_SCREEN) - 1.0f;
    y = y * 2 / (SHOW_SIZE_Y) - 1.0f;
}
static void GridCoor_to_ClipCoor(Float &x, Float &y, SCREENID_T screenid) {
    x *= (SHOW_SIZE_X / GRID_SIZE_Y);
    y *= (SHOW_SIZE_Y / GRID_SIZE_Z);
    x = (x + screenid*SHOW_SIZE_X) * 2 / (SHOW_SIZE_X*TOTAL_SCREEN) - 1.0f;
    y = y * 2 / (SHOW_SIZE_Y)- 1.0f;
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
    P_3d<T> atz(T z) throw(int)
    {
        if (fabs(p2.z - p1.z) < EPS<T>) throw 233;
        T inc = (z - p1.z) / (p2.z - p1.z);
        return P_3d<T>(p1.x + (p2.x - p1.x) * inc, p1.y + (p2.y - p1.y) * inc, z);
    }
};
const Float pixel_size = 1.0f / 8.0f;
using p3df = P_3d<Float>;
using sef = Segment<Float>;
using cuf = Cuboid<Float>;
void init_Density_3d()
{
    p3df vrp(20, 32, 32), vpn(0, 0, -pixel_size), uvp(0, pixel_size, 0), vvp(vpn.cross(uvp)), prp(30, 32, 32);
    cuf jar(p3df(0, 0, 0), p3df(GRID_SIZE_X, GRID_SIZE_Y, GRID_SIZE_Z));
    for (int i = 0; i < SHOW_SIZE_X; ++i)
        for (int j = 0; j < SHOW_SIZE_Y; ++j)
            try
            {
                sef ray(prp, vrp + uvp * i + vvp * j);
                ray.p2 = ray.atz(0); //保证光线不能与OXY平面平行
                if (!ray.Clip(jar)) throw 2333;
                while (ray.len() > EPS<Float>)
                {
                    sef thu = ray;
                    thu.p1 = thu.p1 + (thu.p2 - thu.p1) * (0.1 / thu.len());

                    int x = (thu.p1.x), y = thu.p1.y, z = thu.p1.z;
                    thu.Clip(cuf(p3df(x, y, z), p3df(x + 1, y + 1, z + 1)));
                    ray.p2 = thu.p1;
                }
            } catch(int err) { fprintf(stderr, "err: %d\n", err); }


}
vector<pair<int, Float> > lightPath[SHOW_SIZE_X][SHOW_SIZE_Y];


void Draw_Density_3d(Float *density, SCREENID_T screen)
{
    glBegin(GL_POINTS);
    for (int i = 0; i < SHOW_SIZE_X; ++i)
        for (int j = 0; j < SHOW_SIZE_Y; ++j)
        {
            Float rgb(0), a(0);
            for (auto x : lightPath[i][j])
            {
                if (a > 0.995) break;
                Float rd = density[x.first] * x.second;
                rgb += rd * rd * (1 - a);
                a += (1 - a) * rd;
            }
            glColor3f(rgb, rgb, rgb);
            glVertex2f(i, j);
        }
    glEnd();
    glFlush();
}

void Draw_Density_2d(Float *density, SCREENID_T screenid) {
    printf("draw density\n");
    assert(screenid < 5);
    glBegin(GL_POINTS);
    for (int is = 0; is < SHOW_SIZE_X; is++) {
        for (int js = 0; js < SHOW_SIZE_Y; js++) {
            int i = is*GRID_SIZE_Y / SHOW_SIZE_X;
            int j = js*GRID_SIZE_Z / SHOW_SIZE_Y;
            //printf("%d %d\n", i, j);
            //for (int k = 0; k < 3; k++) {
            assert(0 <= density[ID(2, i, j)] && density[ID(2, i, j)] <= 1.0);
            Float x = is, y = js;
            ScreenCoor_to_ClipCoor(x, y, screenid);
            Float d = density[ID(2, i, j)];
            glColor3f(d, d, d);
            glVertex2f(x, y);
            //pixels[p++] = BYT(density[ID(2, i, j)] * 255);
            //}
            //pixels[p++] = (BYT)255;
        }
    }
    glEnd();
    glFlush();
    //printf("pixel: %d\n", pixels[0]);
    //glDrawPixels(SHOW_SIZE_X, SHOW_SIZE_Y, GL_RGBA, GL_UNSIGNED_BYTE, pixels);
}

void Draw_Velocity_2d(int i, Float*, Float *vys, Float *vzs, SCREENID_T screenid) {
    printf("draw velocity\n");
    for (int j = 0; j < GRID_SIZE_Y; j ++) {
        for (int k = 0; k < GRID_SIZE_Z; k ++) {
            int j0 = j, k0 = k;
            Float vy = vys[ID(i, j0, k0)] / 10.0f, vz = vzs[ID(i, j0, k0)] / 10.0f;
            //printf("%f %f\n", vy, vz);
            Float y0 = j0, z0 = k0;
            Float y1 = j0 + vy, z1 = k0 + vz;
            //printf("%f %f %f %f\n", y0, z0, y1, z1);
            GridCoor_to_ClipCoor(y0, z0, screenid);
            GridCoor_to_ClipCoor(y1, z1, screenid);
            //y0 = y0 / SHOW_SIZE_X, y1 = y1 / SHOW_SIZE_X;
            //z0 = z0 / SHOW_SIZE_Y * 2 - 1, z1 = z1 / SHOW_SIZE_Y * 2 - 1;
            //y0 += SHOW_SIZE_X, y1 += SHOW_SIZE_X;
            glColor3f(1.0, 1.0, 1.0);
            glBegin(GL_LINES);
            //printf("lines\n");
            glVertex2f(y0, z0);
            glVertex2f(y1, z1);
            glEnd();
            glColor3f(0, 1, 0);
            glBegin(GL_POINTS);
            glVertex2f(y0, z0);
            glEnd();
        }
    }
    glFlush();
}

