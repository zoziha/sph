#! /mingw64/bin/lua

-- 根据CFL条件，计算时间步
-- 目前，`sph` 项目暂不支持自适应时间步长，语言预先定义了一个时间步长

io.write("这是计算时间步长的程序\n")
io.write("程序采用SI国际单位, 即秒为单位\n")
io.write("作者: 左志华\n")
io.write("版本: v1.0\n\n")

io.write("请输入光滑长度 (0.001?): ")
local l0 = io.read("n")

io.write("通常为保证流体密度的相对变化率不大于 1%, 声速应大于流场最大速度的 10 倍\n")
io.write("请输入声速 (1480? 80? 340?): ")
local c = io.read("n")

io.write("根据CFL条件 (CFL = t*c/l0 <= 0.2?), 时间步长上限为:")
io.write(0.2 * l0 / c, "\n")
