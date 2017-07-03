场景文件格式：
#...：注释。一个被注释的行必须以'#'开头。
gridsize GRIDX GRIDY GRIDZ # 指定模拟网格的大小。如果没有，则默认为256*256*256. 以下所有参数均应以256*256*256的网格为基准，程序会按照真实大小进行自动缩放。
end n #指定运算n帧后终止运行。如果没有此指令，则程序不会停止。
box x0 x1 y0 y1 z0 z1 #在初始时候增加一个装满水的长方体，坐标为闭区间[x0,x1],[y0,y1],[z0,z1]。它不会覆盖固体格子和网格外的地方。水的初始速度为0。
source x0 x1 y0 y1 z0 z1 vx vy vz end#增加一个水源，同上。vx vy vz是它吐出水的初始速度，它在第end+1帧将不再刷出水。
periodbox x0 x1 y0 y1 z0 z1 semi_period vx0 vy0 vx1 vy1 vz0 vz1 end # 增加一块速度随时间周期振荡的区域。这块区域的位置由x0等指定，半周期由semi_period指定，单位为帧，含义是：在第0,2,4...个半周期内会设定vx0,vy0,vz0的速度，第1,3,5...个半周期则设置vx1,vy1,vz1的速度。在时间end之后不再生效。（此功能暂不支持）

现在todo（星号为不重要）：
1. 监控mask内部出气泡的问题。
2. 渲染（wukan）。用blender python API出精美demo，自己实现一份光线追踪，作为代码结构上的补齐。如果wkw组photon mapping能出成果，亦可借用之。
3. 支持水源的初始速度：支持从水源particle的速度直接出网格速度。调研为什么此事挂掉了。
4. 调研GridFluidSim，出飞沫等额外特效。
5. better method for sampling the marker particles
6*. Marching Cubes支持normal，材料：http://www.angelfire.com/linux/myp/MCAdvanced/MCImproved.html
7. zyh: Catmull-Rom interpolation in advection
8. zyh: FLIP method, tune parameters
9. 更好的pressure solve？
10. 解决side.box的interpolation goes wrong问题。

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

2017.5.30
实现了MAC+PIC方法。修正了projection和advection的bug，修正了模型内部气泡问题。接下来首要事项是demo艺术。

2017.5.31
发现bug：advection自更新导致的速度耗散，在projection后产生了pressure爆炸。此外，把固体内部面的速度改成0会影响PIC advection。

2017.6.2
修正了extrapolation的自更新问题。因此修正了已知的所有bug。

2017.6.3
更新了场景文件格式。将GRIDXYZ改为可读取。场景进一步和代码解耦。现在运行时间较长，需要优化。

2017.6.4
修改了水源函数，改成正确的“如果是空气则刷新”。运用了OpenMP进行优化。

2017.6.22
今晚跑freefall.aligned的渲染，在demo.freefall.aligned里

有趣的材料： [http://gamedev.stackexchange.com/questions/177/what-is-some-good-examples-about-creating-2d-fluids](http://gamedev.stackexchange.com/questions/177/what-is-some-good-examples-about-creating-2d-fluids)
