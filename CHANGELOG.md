# Change logs

## v0.0.1

- 添加 fpm 和 visual studio 支持
- 支持核函数:
  - cubic spline kernel;
  - gauss kernel;
  - quintic kernel.
- 支持空间维度: 1~3D

## v0.0.2

- 重构接口 (Refactor the Interface)
  - 所有函数参数接口
- 将许可证改为 BSD 3-Clause
- 添加 vtk 输出，支持 paraview 可视化
- 使用自制子组件库: easy_string, master_time, progress_bar

## TO DO

- 树形搜索法 (2D) (WIP)
- 添加进度条 (WIP)
- 测试
- 具有材料强度的动力学
- OpenMP 并行
- 前后处理模块增强
- 三维算例计算, 更好地支持 3D
- 添加读取 BDF 能力 (https://gitee.com/zoziha/set-particles)
