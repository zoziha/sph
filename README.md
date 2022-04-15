# 光滑粒子流体动力学（SPH）([English](./README_EN.md))

一份社区驱动的开源光滑粒子流体动力学（SPH）代码，起始代码版本源自[课本《光滑粒子流体动力学--一种无网格粒子法》](doc/books/光滑粒子流体动力学：一种无网格粒子法.pdf)。

[![wakatime](https://wakatime.com/badge/user/ca8e3153-da86-47e8-ba89-1fac0c842c19.svg)](https://wakatime.com/@ca8e3153-da86-47e8-ba89-1fac0c842c19)
[![Netlify Status](https://api.netlify.com/api/v1/badges/49b0928e-4dc4-4c1c-9aa1-b3f0095fb152/deploy-status)](https://app.netlify.com/sites/zoziha-sph-api-docs/deploys)
[![Netlify Status](https://api.netlify.com/api/v1/badges/ba328ec0-0f67-4349-86b1-2d5f5847c98a/deploy-status)](https://app.netlify.com/sites/zoziha-sph-specs-and-tutorial/deploys)

## 开始

### 依赖

- Lua链接库 >= 5.3

如果是 visual-studio 用户，需要将[编译的 Lua 官方库](https://gitee.com/zoziha/sph/issues/I5138J#note_9613327_link)，放置在 `solution` 文件夹下。

### 获取代码

```sh
git clone https://gitee.com/zoziha/SPH.git
cd SPH
```

### 使用[fortran-lang/fpm](https://github.com/fortran-lang/fpm)构建代码

FPM是社区驱动的Fortran语言的包管理器和代码构建工具，适用于C/C++/Fortran代码的构建。  
你可以通过提供的`fpm.toml`构建代码：

```sh
fpm build
# 运行SPH主程序
fpm run sph --profile release
# 运行后处理程序，生成vtk给ParaView使用
fpm run vtk
fpm test
```

### 其它构建系统

除了fpm，本项目将有可能支持`visual-studio / code::blocks`进行构建，具体的解决方案见`./solution`文件夹。

```sh
fpm update      # 仍然需要使用fpm拉取上游依赖库!
F5 / CTRL + F5  # 运行与调试
```

然后，在visual-studio中检查和选中`src`和`build/dependencies/*`中的源码文件，并构建和运行。

### [FORD](https://github.com/Fortran-FOSS-Programmers/ford) API 文档

```sh
ford FORD-project-file.md
```

最新生成、发布的API文档详见[此链接](https://zoziha-sph-api-docs.netlify.app/)。

### 规范与教程

[mdbook](https://github.com/rust-lang/mdBook) 是一个 markdown 文档生成器。

```sh
cd doc && mdbook build
```

最新生成、发布的规范、教程文档详见[此链接](https://zoziha-sph-specs-and-tutorial.netlify.app/)。

## 链接

- [spheric/SPH Codes](https://spheric-sph.org/sph-projects-and-codes)
- [zoziha/SPH-homework](https://github.com/zoziha/SPH-homework)
