# 日志器

`sph` 项目使用 `stdlib_logger` 来记录日志，在不出现循环 bug 的情况下，增加日志器
并不会减慢程序的高效运行，尤其是避免在核心计算部分的循环内部记录日志。

同样的，`toml` 配置文、`lua` 前处理脚本都在一定程度上不影响核心计算效率，原因是它们的工作量
属于比较小的量级，但又十分重要，特别是面向用户体验。