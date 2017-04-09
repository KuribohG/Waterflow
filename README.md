FRESH POTS(TO BE CONTINUE...) （排名不分先后）

路线图：extended 2D(v1.0)---better extended 2D---naive 3D---better and faster 3D

1. 把Semi-Lagrangian换成MacCormack
2. 支持水源和外力
3. 搞清楚[http://cowboyprogramming.com/2008/04/01/practical-fluid-mechanics/](http://cowboyprogramming.com/2008/04/01/practical-fluid-mechanics/) （一个晚于GDC03的系统）干了什么，有没有可能加入我们的系统 catch：它没有强制散度为零，而是用特殊的advection技巧，以及一种根据压强算速度的方法把散度变得像零。这个是“主体思想”，是假的社会主义理想，不予考虑
4. 比较"Real-Time Eulerian Water Simulation Using a Restricted Tall Cell Grid"（共产主义理想）比我们多了什么
5. 研究怎么支持3D情况，搞清楚[https://github.com/rlguy/GridFluidSim3D](https://github.com/rlguy/GridFluidSim3D) （社会主义理想）他们干了什么 catch：这和任务7相同
6. 支持鼠标交互，彻底搞懂v1.0的逻辑
7. 泛读fluid simulation for computer graphics(Robert Bridson)，对每一个子任务找出我们我们可能的优化方向 catch：
   * ~~advection：可以运用Runge-Kutta法改进semi-laglarian法~~
   * p53,lin_solve的理论依据，这一过程可以写成一个真正的线性方程组，用conjugate gradient算法求解
8. 用假的光线追踪（固定相机视角，具体方法见GPU Gem3）支持三维渲染
9. 支持free surface和重力
10. 更新advection算法
11. 阅读论文：Fast Grid-Free Surface Tracking，研究如何解决surface-tracking问题（着重寻找此领域的开山鼻祖paper）

有趣的材料： [http://gamedev.stackexchange.com/questions/177/what-is-some-good-examples-about-creating-2d-fluids](http://gamedev.stackexchange.com/questions/177/what-is-some-good-examples-about-creating-2d-fluids)