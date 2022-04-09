# 开发工具

本项目目前在 Windows10 操作系统上使用 `vs code` + `msys2 gfortran` 和 `visual-studio` + `ifort` 进行开发，
保证程序的兼容性。

但相信可以方便迁移至 Linux 环境。其中 GFORTRAN 环境易于迁移，INTEL 环境需要注意以下三点:

- 预处理设置;
- 堆栈设置;
- 数组边界检测.

## 构建系统

目前，`sph` 项目支持以下构建系统：

- `fortran-lang/fpm`;
- `visual-studio`.

并且在未来有限的时间内，不打算支持更多的构建系统。

## 代码文档生成器

目前，`sph` 项目使用 [ford](https://github.com/Fortran-FOSS-Programmers/ford) 构建代码文档。

```sh
ford FORD-project-file.md
```

最新的 `ford` 生成的静态网页详见[此链接](https://zoziha-sph-api-docs.netlify.app/)。

## 规范与教程

`ford` 主要针对代码文档，面对开发文档，尚属欠缺；同时为了方便根据开发 markdown 文档生成 pdf 文档，
`sph` 项目也提供了 mdbook 来构建《sph 规范与教程》，最新的静态网页详见
[此链接](https://zoziha-sph-specs-and-tutorial.netlify.app/)。

```sh
cd doc && mdbook build
```

## 参考链接

- [ford](https://github.com/Fortran-FOSS-Programmers/ford)
- [mdbook](https://github.com/rust-lang/mdBook)
