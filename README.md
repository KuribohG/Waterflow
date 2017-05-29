场景文件格式：
#...：注释，只支持整行注释，不支持行内注释。
box x0 x1 y0 y1 z0 z1 #在初始时候增加一个装满水的长方体，坐标为闭区间[x0,x1],[y0,y1],[z0,z1]。它不会覆盖固体格子和网格外的地方。水的初始速度为0。
source x0 x1 y0 y1 z0 z1 rate vx vy vz #增加一个水源，同上，其中每个格子，每一秒平均生成rate次水。vx vy vz是它吐出水的初始速度（暂且不支持）。

现在todo（星号为不重要）：
1. 渲染（wukan）。可选方案：blender API，自写光线追踪，借用wkw组photon mapping。
2. 全面切换到PIC advection方法，不再采用grid advection。
3. 修extrapolation的bug：x,y,z=0时不计算。
4. （或许在PIC上线之后）设法支持水源的初始速度。
5*. 监控projection的表现，之前提到的“在解方程之前apply boundary conditions”是否有必要。
6*. 支持出surface mesh：mesh->TSDF，TSDF场的advection
7*. 把marker particle增加随机扰动（这个有优化，但不用急着做，因为现在对齐的marker particle更利于debug）
8*. Marching Cubes支持normal，材料：http://www.angelfire.com/linux/myp/MCAdvanced/MCImproved.html

路线图：2D(dummies)---extended 2D(v1.0)---better extended 2D---naive 3D---better and faster 3D

get started:
v1.0（盛满水的鱼缸里的染料）
1. 论文：Real-Time Fluid Dynamics for Games，Stam
2. 一个网页：https://mikeash.com/pyblog/fluid-simulation-for-dummies.html
3. 网页里的代码：https://mikeash.com/pyblog/blog/images/fluid.c.file
v2.0（空气和水）
1. 书本：fluid simulation for computer graphics robert bridson
2. 一个项目：https://github.com/rlguy/GridFluidSim3D

2017.5.5
现在状态：除了mesh和extrapolation之外，实现了gridfluidsim3D的全部主要要素。（飞沫等不在计划中，略去）

2017.5.29
状态：实现了MAC方法的全要素（除飞沫等），但MAC advection出现问题，准备切换到PIC advection。

有趣的材料： [http://gamedev.stackexchange.com/questions/177/what-is-some-good-examples-about-creating-2d-fluids](http://gamedev.stackexchange.com/questions/177/what-is-some-good-examples-about-creating-2d-fluids)