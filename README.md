# 光滑粒子流体动力学（SPH）

一份社区驱动的开源光滑粒子流体动力学（SPH）代码，起始代码版本源自课本《光滑粒子流体动力学--一种无网格粒子法》。

| 项目 | 描述 |
| --- | --- |
| 版本 | 0.0.1 |
| 许可证 | Public Domain |
| 版权 | Copyright (c) 2021 SPH 贡献者 |

## 开始

### 获取代码

```sh
git clone https://github.com/zoziha/SPH.git
cd SPH
```

### 使用Make构建代码

本项目支持Make工具构建代码：

```sh
make
```

### 使用[fortran-lang/fpm](https://github.com/fortran-lang/fpm)构建代码

FPM是社区驱动的Fortran语言的包管理器和代码构建工具，适用于c/c++/fortran代码的构建。  
你可以通过提供的`fpm.toml`构建代码：

```sh
fpm build
fpm run
```

### 其它构建系统

除了make和fpm，本项目还支持cmake和visual studio进行构建。

## 链接

+ [spheric/SPH Codes](https://spheric-sph.org/sph-projects-and-codes)
