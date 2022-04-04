# 开发工具

本项目目前在 Windows10 操作系统上使用 vs code + msys2 gfortran 和 visual-studio + ifort 进行开发，保证程序的兼容性。

但相信可以方便迁移至 Linux 环境。其中 GFORTRAN 环境易于迁移，INTEL 环境需要注意以下三点:

- 预处理设置;
- 堆栈设置;
- 数组边界检测.
