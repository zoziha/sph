# output_m

该模块是拓展原始代码的第一个模块，包含了一些程序输出的例程。

## 单元例程说明

### `output_all`

#### 描述

输出保存每个保存时间步的求解信息（拓展）。

#### 语法

`call [[output_md(module):output_all(subroutine)]](x, vx, mass, rho, p, u, c, itype, hsml, ntotal, n)`

#### 参数

`x/vx`: 浮点型的rank-2数组。
该参数是`intent(in)`的。<br>
粒子的坐标/速度。

`mass/rho/p/u/c`: 浮点型的rank-1数组。
该参数是`intent(in)`的。<br>
粒子的质量/密度/压强/内能/声速。

`itype`: 整型的标量。
该参数是`intent(in)`的。<br>
粒子的类型。

`hsml`: 浮点型的rank-1数组。
该参数是`intent(in)`的。<br>
粒子的光滑长度。

`ntotal`: 整型的标量。
该参数是`intent(in)`的。
粒子的总数。

`n`: 整型的标量。
该参数是`intent(in)`的。
打印的时间戳数量。

#### 输出

在目标输出路径的`all`文件夹输出若干时间步的`dat`文件，包含：粒子的位置、速度、质量、密度、压强、内能、粒子类型、光滑长度。