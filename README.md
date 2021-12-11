# Smoothed Particle Hydrodynamics（SPH）([中文文档](./README_CN.md))

A community-driven open source smoothed particle hydrodynamics (SPH) code, the starting code version is derived from the textbook "Smooth Particle hydrodynamics - A Meshfree Particle Method".

| Item | Info |
| --- | --- |
| Version | 0.0.1 |
| License | Public Domain |
| Copyright | Copyright (c) 2021 SPH Contributors |

## Start

### Get the Code

```sh
git clone https://github.com/zoziha/SPH.git
cd SPH
```

### Build with CMake

This project supports CMake tool to build code:

```sh
mkdir build && cd build
cmake -G "MSYS Makefiles" ..
make
```

### Build with [fortran-lang/fpm](https://github.com/fortran-lang/fpm)

FPM is a community-driven Fortran language package manager and code build system, suitable for c/c++/fortran code construction.
You can build the code with the provided `fpm.toml`:

```sh
fpm build
fpm run
```

## Links

+ [spheric/SPH Codes](https://spheric-sph.org/sph-projects-and-codes)
