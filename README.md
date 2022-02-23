# 光滑粒子流体动力学（SPH）([English README](./README_EN.md))

[![fpm](https://github.com/zoziha/SPH/workflows/fpm/badge.svg)](https://github.com/zoziha/SPH/actions)
[![msys2-fpm](https://github.com/zoziha/SPH/workflows/msys2-fpm/badge.svg)](https://github.com/zoziha/SPH/actions)

一份社区驱动的开源光滑粒子流体动力学（SPH）代码，起始代码版本源自课本《光滑粒子流体动力学--一种无网格粒子法》。

| 项目 | 描述 |
| :-: | :-: |
| 版本 | 0.0.3 |
| 许可证 | BSD 3-Clause |
| 版权 | Copyright (c) 2021~2022 SPH 贡献者 |

## 开始

### 获取代码

```sh
git clone https://github.com/zoziha/SPH.git
cd SPH
```

### 使用[fortran-lang/fpm](https://github.com/fortran-lang/fpm)构建代码

FPM是社区驱动的Fortran语言的包管理器和代码构建工具，适用于c/c++/fortran代码的构建。  
你可以通过提供的`fpm.toml`构建代码：

```sh
fpm build
# 运行SPH主程序
fpm run sph --profile release
# 运行后处理程序，生成vtk给ParaView使用
fpm run vtk
```

### 其它构建系统

除了fpm，本项目将有可能支持visual studio进行构建。

## 链接

+ [spheric/SPH Codes](https://spheric-sph.org/sph-projects-and-codes)
