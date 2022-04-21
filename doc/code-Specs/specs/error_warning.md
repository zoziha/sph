# 错误与警告

目前“SPH”没有日志器功能，采用简单的示错预警功能。

“SPH”的错误提示代码样式为：

```fortran
use config_m, only: stderr
write (stderr, fmt='(a, g0.4)') "错误提示字符", error_value
error stop "*<module::procedure>*"
```
