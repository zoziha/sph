---
title: 开发者日志
---

## v0.0.1 - 0.1.2

- 添加 fpm 和 visual studio 支持
- 支持核函数:
  - cubic spline kernel;
  - gauss kernel;
  - quintic kernel.
- 支持空间维度: 1~3D
- 重构接口 (Refactor the Interface)
  - 所有函数参数接口
- 将许可证改为 BSD 3-Clause
- 添加 vtk 输出，支持 paraview 可视化
- 使用自制子组件库: easy_string, master_time, progress_bar

## v0.1.2 -

- 对非Windows-Intel-OneAPI环境采用部分esc码彩色样式;
- 拟支持toml-f解析配置文件;
- 移除原始output子例程;
- 使用error_stop(base), info;
- 数型搜索法(2D);
- 添加进度条.

## TO DO

- 树形搜索法 (1, 3D)
- 支持配置文件toml (WIP)
- 测试
- 具有材料强度的动力学
- OpenMP 并行
- 前后处理模块增强
- 三维算例计算, 更好地支持 3D
- 添加读取 BDF 能力 (https://gitee.com/zoziha/set-particles)
