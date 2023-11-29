# Smoothed Particle Hydrodynamics（SPH）([中文文档](./README_CN.md))

[![fpm](https://github.com/zoziha/SPH/workflows/fpm/badge.svg)](https://github.com/zoziha/SPH/actions)
[![msys2-fpm](https://github.com/zoziha/SPH/workflows/msys2-fpm/badge.svg)](https://github.com/zoziha/SPH/actions)

A community-driven open source smoothed particle hydrodynamics (SPH) code, the starting code version is derived from the textbook "Smooth Particle hydrodynamics - A Meshfree Particle Method".

| Item | Info |
| :-: | :-: |
| Version | 0.0.2 |
| License | BSD 3-Clause  |
| Copyright | Copyright (c) 2021 SPH Contributors |

## Get Started

### Get the Code

```sh
git clone https://github.com/zoziha/SPH.git
cd SPH
```

### Build with [fortran-lang/fpm](https://github.com/fortran-lang/fpm)

FPM is a community-driven Fortran language package manager and code build system, suitable for c/c++/fortran code construction.
You can build the code with the provided `fpm.toml`:

```sh
cd data && mkdir all && mkdir paraview
fpm build
# SPH main program
fpm run sph --profile release
# SPH post-processing program for ParaView
fpm run vtk
```

## Links

+ [spheric/SPH Codes](https://spheric-sph.org/sph-projects-and-codes)
