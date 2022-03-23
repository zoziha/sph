# 光滑粒子流体动力学（SPH）([English](./README_EN.md))

一份社区驱动的开源光滑粒子流体动力学（SPH）代码，起始代码版本源自[课本《光滑粒子流体动力学--一种无网格粒子法》](doc/books/光滑粒子流体动力学：一种无网格粒子法.pdf)。

## 开始

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
```

### 其它构建系统

除了fpm，本项目将有可能支持visual-studio进行构建。

```sh
fpm update # 仍然需要使用fpm拉取上游依赖库!
```

然后，在visual studio中检查和选中`src`和`build/dependencies/*`中的源码文件，并构建和运行。

## 链接

+ [spheric/SPH Codes](https://spheric-sph.org/sph-projects-and-codes)
