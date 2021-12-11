# gprof

gprof可以用来分析程序的性能。

## 简单使用教程

[参考教程](https://zhuanlan.zhihu.com/p/385842627)

1. 使用编译参数`-pg`；
2. 执行可执行程序；
3. 查看gmon.out。

![SPH运行效率](images/gprof10-2.png)

显然，我们能看出来，粒子搜索法占用了非常多的运行时间。