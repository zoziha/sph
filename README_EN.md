# Smoothed Particle Hydrodynamics（SPH）([中文文档](./README.md))

A community-driven open source smoothed particle hydrodynamics (SPH) code, the starting code version is derived from the textbook "Smooth Particle hydrodynamics - A Meshfree Particle Method".

## Getting Started

### Get the Code

```sh
git clone https://gitee.com/zoziha/SPH.git
cd SPH
```

### Build with [fortran-lang/fpm](https://github.com/fortran-lang/fpm)

FPM is a community-driven Fortran language package manager and code build system, suitable for c/c++/fortran code construction.
You can build the code with the provided `fpm.toml`:

```sh
fpm build
# SPH main program
fpm run sph --profile release
# SPH post-processing program for ParaView
fpm run vtk
```

## Links

- [spheric/SPH Codes](https://spheric-sph.org/sph-projects-and-codes)
- [zoziha/SPH-homework](https://github.com/zoziha/SPH-homework)
