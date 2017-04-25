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

2017.4.25
现在状态：
advect：Runge_Kutta
Project：书上的解线性方程组

现在todo：
1. 支持出surface mesh

现在问题：
1. “贴壁”：贴着壁的速度总是0，这个问题是project之前的，换言之，是apply_external_forces（或者advect，但更可能是前者）出了问题
2. “悬浮”：有一些marker particle悬浮在半空中，它们的速度为0，因此就一直那么悬在空中了


有趣的材料： [http://gamedev.stackexchange.com/questions/177/what-is-some-good-examples-about-creating-2d-fluids](http://gamedev.stackexchange.com/questions/177/what-is-some-good-examples-about-creating-2d-fluids)