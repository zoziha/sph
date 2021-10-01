# 使用GNUPlot作为可视化工具

GNUPlot是一款免费、开源、交互型的数据可视化软件。<br>
我们可以使用它简单地进行数据的可视化。

## 绘制粒子位置

```gnuplot
set title "Particles Location Info"
set xlabel "X"
set ylabel "Y"
unset key

# 打印单帧粒子位置：初始化位置
plot "ini_xv.dat" using 2:3 w p pt 7

# 循环打印粒子位置
do for [i=0:20] { plot sprintf('f_%dxv.dat', i) using 2:3 w p pt 6; pause 0.5 }
```

![初始化粒子位置信息](images/ini_x.PNG)