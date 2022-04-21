# toml 配置文件

配置文件是控制台程序重要的人机交互方式之一。

`sph.toml` 文件可配置的参数信息为：

参数 `parameter`:

- `visc_artificial`: 人工粘性，默认为 F;
- `heat_artificial`: 人工热量，默认为 F;
- `visc`: 粘性，默认为 T;
- `self_gravity`: 自重，默认为 F;
- `nnps`: 粒子搜索算法，默认为 1-直接搜索;
- `skf`: 光滑核函数，默认为 1-三次样条核函数;
- `name`: 项目名称, 默认为 untitled;

前处理 `pre-process`:

- `dofile`: 自定义粒子初始化 Lua 脚本，默认为 F;

详见 `./src/config/toml_info.f90` 文件。
