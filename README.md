FRESH POTS(TO BE CONTINUE...) （排名不分先后）

路线图：extended 2D(v1.0)---better extended 2D---naive 3D---better and faster 3D

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

现在todo：
1. 渲染。可选方案：blender API，自写raytracing，wkw组photon mapping。
2. 调研第10章，研究hybrid方法
3. 支持出surface mesh：mesh->TSDF，TSDF场的advection
4. 把marker particle增加随机扰动（这个有优化，但不用急着做，因为现在对齐的marker particle更利于debug）
5. 搞懂blender怎么渲染实心玻璃材质。
6. Marching Cubes支持normal，材料：http://www.angelfire.com/linux/myp/MCAdvanced/MCImproved.html



有趣的材料： [http://gamedev.stackexchange.com/questions/177/what-is-some-good-examples-about-creating-2d-fluids](http://gamedev.stackexchange.com/questions/177/what-is-some-good-examples-about-creating-2d-fluids)