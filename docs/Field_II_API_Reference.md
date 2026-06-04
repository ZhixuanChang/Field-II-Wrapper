# Field-II API Reference

> 基于 Field II User Guide (Release 3.30, April 5, 2021) 整理
> 用于 Field-II-Wrapper 项目开发参考

---

## 目录

1. [程序架构概述](#1-程序架构概述)
2. [通用命令 (General Commands)](#2-通用命令-general-commands)
3. [场参数设置 (Field Configuration)](#3-场参数设置-field-configuration)
4. [换能器定义 (Transducer Definition)](#4-换能器定义-transducer-definition)
5. [换能器操控 (Transducer Manipulation)](#5-换能器操控-transducer-manipulation)
6. [阵元操控 (Element Manipulation)](#6-阵元操控-element-manipulation)
7. [场计算 (Field Calculation)](#7-场计算-field-calculation)
8. [典型仿真流程](#8-典型仿真流程)
9. [Wrapper 项目中的映射关系](#9-wrapper-项目中的映射关系)
10. [关键数据结构](#10-关键数据结构)

---

## 1. 程序架构概述

Field-II 由 C 语言核心程序和 Matlab M-函数接口组成：

- **计算引擎**: C 程序执行所有声场计算，管理所有数据
- **M-函数接口**: 三类命名约定
  - `field_*` — 初始化和系统命令
  - `xdc_*` — 换能器定义和操控 (trans**d**u**c**er)
  - `calc_*` — 声场计算

**核心概念**:
- **空间脉冲响应 (Spatial Impulse Response)**: 换能器被 Dirac 脉冲激励时，空间某点的声场随时间的变化
- **Born 近似**: 忽略多次散射，仅计算单次散射
- **物理阵元 vs 数学阵元**: 物理阵元被细分为更小的数学阵元以提高计算精度
- **时间线 (Time Line)**: 聚焦和变迹通过时间线动态控制

**初始默认值** (field_init 后):

| 参数 | 值 | 说明 |
|------|------|------|
| `c` | 1540 m/s | 声速 |
| `fs` | 100 MHz | 采样频率 |
| `use_att` | 0 | 不使用衰减 |
| `att` | 0.0 dB/m | 频率无关衰减 |
| `freq_att` | 0.0 dB/[m·Hz] | 频率相关衰减 |
| `att_f0` | 0.0 Hz | 衰减中心频率 |
| `use_rectangles` | 1 | 使用矩形描述孔径 |
| `use_triangles` | 0 | 使用三角形描述孔径 |
| `use_lines` | 0 | 使用线段描述孔径 |
| `fast_integration` | 0 | 使用 Romberg 积分 |

---

## 2. 通用命令 (General Commands)

### 2.1 field_init — 初始化 Field-II

```
field_init(suppress)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `suppress` | int (可选) | 0: 不显示启动画面; -1: 无 ASCII 输出 |

**说明**: 必须在使用任何其他 Field-II 函数之前调用。设置默认参数值。

---

### 2.2 field_end — 终止 Field-II

```
field_end
```

**说明**: 释放 Field-II 占用的所有存储空间。

---

### 2.3 field_debug — 调试开关

```
field_debug(state)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `state` | int | 1: 开启调试输出; 0: 关闭 |

---

### 2.4 field_info — 显示系统信息

```
field_info
```

打印当前 Field-II 配置状态（版本、孔径数、声速、采样频率、衰减参数等）。

---

### 2.5 field_guide — 显示用户手册

```
field_guide
```

调用 Acrobat Reader 显示 Field-II 用户手册 PDF。

---

## 3. 场参数设置 (Field Configuration)

### 3.1 set_field — 设置场参数

```
set_field(option_name, value)
```

| 选项名 | 值类型 | 说明 |
|--------|--------|------|
| `'use_att'` | int | 是否启用衰减 (≠0 为启用) |
| `'att'` | float | 频率无关衰减 [dB/m] |
| `'freq_att'` | float | 频率相关衰减 [dB/(m·Hz)]，以 att_f0 为中心 |
| `'att_f0'` | float | 衰减中心频率 [Hz] |
| `'debug'` | int | 调试信息 (1=开启) |
| `'c'` | float | 声速 [m/s] |
| `'fs'` | float | 采样频率 [Hz] |
| `'show_time'` | int | 显示计算耗时 (>2 为打印间隔秒数) |
| `'use_rectangles'` | int | 使用矩形孔径 (1=是) |
| `'use_triangles'` | int | 使用三角形孔径 (1=是) |
| `'use_lines'` | int | 使用线段孔径 (1=是) |
| `'fast_integration'` | int | 快速梯形积分 (1) vs Romberg 积分 (0) |

**示例**:
```matlab
set_field('att', 1.5 * 100);           % 1.5 dB/cm
set_field('freq_att', 0.5 * 100 / 1e6); % 0.5 dB/[MHz·cm]
set_field('att_f0', 3e6);              % 3 MHz 中心
set_field('use_att', 1);               % 启用衰减
set_field('fs', 100e6);                % 100 MHz 采样
set_field('c', 1540);                  % 1540 m/s 声速
```

**注意**: 频率无关衰减应等于 `freq_att * att_f0`，否则深部组织衰减可能不正确。

---

### 3.2 set_sampling — 设置采样频率 (已过时)

```
set_sampling(fs)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `fs` | float | 采样频率 [Hz] |

**注意**: 此函数已被 `set_field('fs', fs)` 取代。

---

## 4. 换能器定义 (Transducer Definition)

### 4.1 低级定义函数

#### xdc_rectangles — 矩形孔径定义

```
Th = xdc_rectangles(rect, center, focus)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `rect` | matrix | 矩形信息，每行一个矩形。列：`[no, x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4, apo, width, height, cx,cy,cz]` |
| `center` | matrix | 物理阵元中心坐标，每行一个 `[x, y, z]` |
| `focus` | vector | 固定焦点 `[x, y, z]` |
| `Th` | handle | 孔径句柄 |

**rect 矩阵列索引**:

| 列 | 变量 | 说明 |
|-----|------|------|
| 1 | no | 物理阵元编号 (从1开始) |
| 2-4 | x1,y1,z1 | 第一个角点 (左下) |
| 5-7 | x2,y2,z2 | 第二个角点 |
| 8-10 | x3,y3,z3 | 第三个角点 |
| 11-13 | x4,y4,z4 | 第四个角点 |
| 14 | apo | 该数学阵元的固定变迹值 |
| 15 | width | 数学阵元宽度 (x方向) |
| 16 | height | 数学阵元高度 (y方向) |
| 17-19 | cx,cy,cz | 矩形中心 |

**注意**: 角点必须按顺时针排列。当前实现要求矩形与 x/y 轴对齐 (不支持旋转)。

---

#### xdc_triangles — 三角形孔径定义

```
Th = xdc_triangles(data, center, focus)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `data` | matrix | 三角形信息，每行 `[no, x1,y1,z1, x2,y2,z2, x3,y3,z3, apo]` |
| `center` | matrix | 物理阵元中心坐标 |
| `focus` | vector | 固定焦点 `[x, y, z]` |

---

#### xdc_lines — 线段边界孔径定义

```
Th = xdc_lines(lines, center, focus)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `lines` | matrix | 线段信息，每行 `[no_phys, no_mat, slope, infinity, intersect, above]` |
| `center` | matrix | 物理阵元中心坐标 |
| `focus` | vector | 固定焦点 |

**限制**: 仅适用于 x-y 平面内的平面阵元。

---

### 4.2 单阵元定义

#### xdc_piston — 圆形平面换能器

```
Th = xdc_piston(radius, ele_size)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `radius` | float | 物理阵元半径 [m] |
| `ele_size` | float | 数学阵元尺寸 [m] |
| `Th` | handle | 孔径句柄 |

---

#### xdc_concave — 凹面换能器

```
Th = xdc_concave(radius, focal_radius, ele_size)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `radius` | float | 物理阵元半径 [m] |
| `focal_radius` | float | 几何焦距 [m] |
| `ele_size` | float | 数学阵元尺寸 [m] |

---

### 4.3 一维阵列定义

#### xdc_linear_array — 线性阵列

```
Th = xdc_linear_array(no_elements, width, height, kerf, no_sub_x, no_sub_y, focus)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `no_elements` | int | 物理阵元数 |
| `width` | float | 阵元宽度 (x方向) [m] |
| `height` | float | 阵元高度 (y方向) [m] |
| `kerf` | float | 阵元间距 (x方向) [m] |
| `no_sub_x` | int | x方向数学细分 |
| `no_sub_y` | int | y方向数学细分 |
| `focus` | vector | 固定焦点 `[x, y, z]` |

---

#### xdc_convex_array — 凸面阵列

```
Th = xdc_convex_array(no_elements, width, height, kerf, Rconvex, no_sub_x, no_sub_y, focus)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `Rconvex` | float | 凸面半径 [m] |

---

#### xdc_focused_array — 仰角聚焦线性阵列

```
Th = xdc_focused_array(no_elements, width, height, kerf, Rfocus, no_sub_x, no_sub_y, focus)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `Rfocus` | float | 仰角 (elevation) 聚焦半径 [m] |

---

#### xdc_convex_focused_array — 仰角聚焦凸面阵列

```
Th = xdc_convex_focused_array(no_elements, width, height, kerf, Rconvex, Rfocus, no_sub_x, no_sub_y, focus)
```

**限制**: `π*Rconvex >= (kerf*(no_elements-1) + width*no_elements)`。所有物理尺寸必须为正。

---

### 4.4 二维 / 多排阵列定义

#### xdc_linear_multirow — 多排线性阵列

```
Th = xdc_linear_multirow(no_elem_x, width, no_elem_y, heights, kerf_x, kerf_y, no_sub_x, no_sub_y, focus)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `no_elem_x` | int | x方向物理阵元数 |
| `width` | float | x方向阵元宽度 [m] |
| `no_elem_y` | int | y方向物理阵元数 (行数) |
| `heights` | vector | 各行的高度 [m]，长度为 no_elem_y |
| `kerf_x` | float | x方向阵元间距 [m] |
| `kerf_y` | float | y方向阵元间距 [m] |

**物理阵元编号规则**: 从最负 x, y 坐标开始，先增加 x，再增加 y。

---

#### xdc_focused_multirow — 仰角聚焦多排阵列 (1.5D)

```
Th = xdc_focused_multirow(no_elem_x, width, no_elem_y, heights, kerf_x, kerf_y, Rfocus, no_sub_x, no_sub_y, focus)
```

---

#### xdc_convex_focused_multirow — 仰角聚焦凸面多排阵列

```
Th = xdc_convex_focused_multirow(no_elements, width, heights, kerf, Rconvex, Rfocus, no_sub_x, no_sub_y, focus)
```

**限制**:
- `π*Rconvex >= (kerf*(no_elements-1) + width*no_elements)`
- `sum(heights) + (no_elem_y-1)*kerf_y > 2*Rfocus`

---

#### xdc_2d_array — 二维 (稀疏) 阵列

```
Th = xdc_2d_array(no_ele_x, no_ele_y, width, height, kerf_x, kerf_y, enabled, no_sub_x, no_sub_y, focus)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `no_ele_x` | int | x方向阵元数 |
| `no_ele_y` | int | y方向阵元数 |
| `enabled` | matrix | `(no_ele_x × no_ele_y)` 矩阵，1 表示启用该阵元 |
| | | `enabled(1,1)` 对应 x-y 平面左下角阵元 |

**稀疏阵列编号规则**: 从最负 xy 的启用阵元开始，先增加 x 再增加 y。

---

### 4.5 格式转换

#### xdc_convert — 矩形转三角形描述

```
xdc_convert(Th)
```

数学阵元数量会翻倍 (一个矩形 = 两个三角形)。

---

## 5. 换能器操控 (Transducer Manipulation)

### 5.1 xdc_apodization — 变迹时间线

```
xdc_apodization(Th, times, values)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `Th` | handle | 孔径句柄 |
| `times` | vector | 各变迹区生效时间 [s] |
| `values` | matrix | 变迹值矩阵 `[时间点数 × 物理阵元数]`，每行至少一个非零值 |

---

### 5.2 xdc_focus — 聚焦时间线

```
xdc_focus(Th, times, points)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `Th` | handle | 孔径句柄 |
| `times` | vector | 各聚焦区生效时间 [s] |
| `points` | matrix | 焦点坐标 `[N×3]`，每行 `[x, y, z]` |

**延迟计算原理**:
```
ti = (sqrt((xf-xc)^2 + (yf-yc)^2 + (zf-zc)^2) - sqrt((xi-xc)^2 + (yi-yc)^2 + (zi-zc)^2)) / c
```
其中 (xf,yf,zf) 为焦点，(xc,yc,zc) 为参考中心，(xi,yi,zi) 为第 i 个物理阵元中心。

---

### 5.3 xdc_focus_times / xdc_times_focus — 用户自定义延迟

```
xdc_focus_times(Th, times, delays)
xdc_times_focus(Th, times, delays)   % 推荐使用此版本
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `Th` | handle | 孔径句柄 |
| `times` | vector | 各延迟区生效时间 [s] |
| `delays` | matrix | 延迟值矩阵 `[时间点数 × 物理阵元数]` [s] |

两个函数功能相同，`xdc_times_focus` 为兼容性新增。

---

### 5.4 xdc_dynamic_focus — 动态聚焦

```
xdc_dynamic_focus(Th, time, dir_zx, dir_zy)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `Th` | handle | 孔径句柄 |
| `time` | float | 动态聚焦生效时间 [s] |
| `dir_zx` | float | z-x 平面内的方向角 [rad] |
| `dir_zy` | float | z-y 平面内的方向角 [rad] |

方向角以焦距中心为参考点。

---

### 5.5 xdc_center_focus — 设置聚焦参考中心

```
xdc_center_focus(Th, point)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `Th` | handle | 孔径句柄 |
| `point` | vector | 聚焦参考中心点 `[x, y, z]` |

用于线性阵列成像中移动发射/接收波束原点。

---

### 5.6 xdc_excitation — 设置激励脉冲

```
xdc_excitation(Th, pulse)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `Th` | handle | 孔径句柄 |
| `pulse` | row vector | 激励脉冲采样序列 |

---

### 5.7 xdc_impulse — 设置脉冲响应

```
xdc_impulse(Th, pulse)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `Th` | handle | 孔径句柄 |
| `pulse` | row vector | 脉冲响应采样序列 |

**注意**: 改变采样频率后必须重新设置所有孔径的脉冲。

---

### 5.8 xdc_baffle — 设置障板条件

```
xdc_baffle(Th, soft_baffle)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `Th` | handle | 孔径句柄 |
| `soft_baffle` | int | 1: 软障板 (Rayleigh-Sommerfeld 积分); 0: 刚性障板 (默认) |

软障板: `h_soft(t) = h_rigid(t) * zp / (c*t)`，其中 zp 为场点的 z 坐标。

---

### 5.9 xdc_quantization — 延迟量化

```
xdc_quantization(Th, min_delay)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `Th` | handle | 孔径句柄 |
| `min_delay` | float | 最小延迟步长 [s]; 0 = 不量化 |

**注意**: 设置后需要重新设定聚焦时间线才能生效。不影响用户自定义延迟。

---

### 5.10 xdc_get — 获取孔径数据

```
data = xdc_get(Th, info_type)
```

| info_type | 返回内容 |
|-----------|---------|
| `'rect'` | 矩形数学阵元信息 (见下方结构) |
| `'tri'` | 三角形数学阵元信息 |
| `'focus'` | 聚焦时间线矩阵 |
| `'apo'` | 变迹时间线矩阵 |

**rect 返回矩阵 (每列一个数学阵元)**:

| 行 | 信息 |
|-----|------|
| 1 | 物理阵元编号 |
| 2 | 该物理阵元内的数学阵元编号 |
| 3 | 数学阵元宽度 [m] |
| 4 | 数学阵元高度 [m] |
| 5 | 数学阵元静态变迹值 |
| 6 | xz 平面夹角正切 |
| 7 | yz 平面夹角正切 |
| 8-10 | 数学阵元中心 (x, y, z) |
| 11-22 | 四个角点 (x, y, z) |
| 23 | 数学阵元延迟值 [s] |
| 24-26 | 物理阵元中心 (x, y, z) |

**tri 返回矩阵 (每列一个数学阵元)**:

| 行 | 信息 |
|-----|------|
| 1 | 物理阵元编号 |
| 2 | 数学阵元编号 |
| 3 | 变迹值 |
| 4-6 | 数学阵元中心 |
| 7-15 | 三个角点 |

---

### 5.11 xdc_show — 显示孔径信息

```
xdc_show(Th, info_type)
```

| info_type | 显示内容 |
|-----------|---------|
| `'elements'` | 阵元信息 |
| `'focus'` | 聚焦时间线 |
| `'apo'` | 变迹时间线 |
| `'all'` | 全部信息 (默认) |

参数可选，默认显示全部。

---

### 5.12 xdc_free — 释放孔径

```
xdc_free(Th)
```

释放孔径占用的存储空间。

---

## 6. 阵元操控 (Element Manipulation)

这些函数对单个数学阵元或物理阵元进行细粒度操控。

### 6.1 ele_apodization — 数学阵元变迹

```
ele_apodization(Th, element_no, apo)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `Th` | handle | 孔径句柄 |
| `element_no` | column vector | 要设置变迹的物理阵元编号列表 |
| `apo` | matrix | `[物理阵元数 × 数学阵元总数]` 变迹值矩阵 |

**说明**: 此变迹值会乘以空间脉冲响应，独立于物理阵元的动态变迹。

---

### 6.2 ele_delay — 数学阵元延迟

```
ele_delay(Th, element_no, delays)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `Th` | handle | 孔径句柄 |
| `element_no` | column vector | 物理阵元编号列表 |
| `delays` | matrix | `[物理阵元数 × 数学阵元总数]` 延迟值矩阵 [s] |

**典型用途**: 模拟固定声透镜 (如仰角聚焦透镜)。

---

### 6.3 ele_waveform — 物理阵元波形

```
ele_waveform(Th, element_no, samples)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `Th` | handle | 孔径句柄 |
| `element_no` | column vector | 物理阵元编号列表 |
| `samples` | matrix | `[物理阵元数 × 采样点数]` 波形采样值 |

**典型用途**: 模拟各阵元使用不同激励波形。

---

## 7. 场计算 (Field Calculation)

### 7.1 calc_h — 空间脉冲响应

```
[h, start_time] = calc_h(Th, points)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `Th` | handle | 发射孔径句柄 |
| `points` | matrix | 场点坐标 `[N×3]`，每行 `[x, y, z]` |
| `h` | vector | 空间脉冲响应 [m/s] |
| `start_time` | float | h 第一个样本对应的时间 [s] |

---

### 7.2 calc_hp — 发射声场

```
[hp, start_time] = calc_hp(Th, points)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `Th` | handle | 发射孔径句柄 |
| `points` | matrix | 场点坐标 `[N×3]` |
| `hp` | vector | 发射声压场 |
| `start_time` | float | 起始时间 [s] |

**计算过程**: 空间脉冲响应 卷积 激励波形 卷积 脉冲响应。

---

### 7.3 calc_hhp — 脉冲回波场

```
[hhp, start_time] = calc_hhp(Th1, Th2, points)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `Th1` | handle | 发射孔径句柄 |
| `Th2` | handle | 接收孔径句柄 |
| `points` | matrix | 场点坐标 `[N×3]` |
| `hhp` | vector | 接收电压信号 |
| `start_time` | float | 起始时间 [s] |

---

### 7.4 calc_scat — 散射体回波

```
[scat, start_time] = calc_scat(Th1, Th2, points, amplitudes)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `Th1` | handle | 发射孔径 |
| `Th2` | handle | 接收孔径 |
| `points` | matrix | 散射体坐标 `[N×3]` |
| `amplitudes` | row vector | 散射体幅度，每个散射体一个值 |
| `scat` | vector | 接收电压信号 (所有阵元叠加) |
| `start_time` | float | 起始时间 [s] |

---

### 7.5 calc_scat_multi — 多通道散射回波 (Wrapper 主要使用)

```
[scat, start_time] = calc_scat_multi(Th1, Th2, points, amplitudes)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `Th1` | handle | 发射孔径 |
| `Th2` | handle | 接收孔径 |
| `points` | matrix | 散射体坐标 `[N×3]` |
| `amplitudes` | row vector | 散射体幅度 |
| `scat` | matrix | `[采样点数 × 接收阵元数]` 每列一个接收通道 |
| `start_time` | float | 起始时间 [s] |

**与 calc_scat 区别**: 返回每个接收阵元的独立信号，而非叠加信号。

---

### 7.6 calc_scat_all — 全合成孔径散射回波

```
[scat, start_time] = calc_scat_all(Th1, Th2, points, amplitudes, dec_factor)
```

| 参数 | 类型 | 说明 |
|------|------|------|
| `Th1` | handle | 发射孔径 |
| `Th2` | handle | 接收孔径 |
| `points` | matrix | 散射体坐标 `[N×3]` |
| `amplitudes` | row vector | 散射体幅度 |
| `dec_factor` | int | 输出降采样因子 (fs_out = fs / dec_factor) |
| `scat` | matrix | `[采样点数 × (发射阵元数×接收阵元数)]` |
| `start_time` | float | 起始时间 [s] |

**数据组织**: 信号依次为: 发射1-接收1, 发射1-接收2, ..., 发射2-接收1, ...

**注意**: 32阵元将产生 1024 个信号，数据量很大。降采样需确保降采样后的采样频率不引起混叠。

---

## 8. 典型仿真流程

### 8.1 基本仿真流程

```matlab
% 1. 初始化
field_init;
set_field('fs', fs);
set_field('c', c);

% 2. 创建换能器
tx_aper = xdc_linear_array(N, width, height, kerf, n_sub_x, n_sub_y, focus);
rx_aper = xdc_linear_array(N, width, height, kerf, n_sub_x, n_sub_y, focus);

% 3. 设置脉冲
impulse = sin(2*pi*f0*(0:1/fs:2/f0)) .* hanning(length(t))';
xdc_impulse(tx_aper, impulse);
xdc_impulse(rx_aper, impulse);

excitation = sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation(tx_aper, excitation);

% 4. 仿真
for i = 1:N_emissions
    xdc_apodization(tx_aper, 0, tx_apod(i,:));
    xdc_times_focus(tx_aper, 0, tx_delay(i,:));
    [rf_cur, t_start(i)] = calc_scat_multi(tx_aper, rx_aper, scat_pos, scat_amp);
    rf_data{i} = rf_cur;
end

% 5. 数据整合
[rfdata, t_start] = f_rf_comb(rf_data, t_start, fs);

% 6. 清理
xdc_free(tx_aper);
xdc_free(rx_aper);
```

### 8.2 Wrapper 项目的 calc_rf 高层封装

```matlab
[rfdata, t_start] = calc_rf(fc, fs, tx_aper, rx_aper, tx_delay, tx_apod, scat_pos, scat_amp);
```

`calc_rf` 自动完成步骤 3-5。

### 8.3 线性阵列扫描模式

使用 `xdc_center_focus` 移动波束原点，配合 `xdc_apodization` 控制活动孔径:

```matlab
for i = 1:no_lines
    x = (i - 1 - no_lines/2) * dx;
    xdc_center_focus(emit_aperture, [x 0 0]);
    xdc_focus(emit_aperture, 0, [x 0 z_focus]);
    apo = [zeros(1, i-1) hamming(N_active)' zeros(1, N_elements-N_active-i+1)];
    xdc_apodization(emit_aperture, 0, apo);
    [v, t1] = calc_scat(emit_aperture, receive_aperture, positions, amplitudes);
end
```

### 8.4 相控阵列扫描模式

通过改变聚焦方向实现扇扫:

```matlab
theta = -sector/2;
for i = 1:no_lines
    xdc_focus(emit_aperture, 0, [70*sin(theta) 0 70*cos(theta)]/1000);
    [v, t1] = calc_scat(emit_aperture, receive_aperture, point_position, 1);
    theta = theta + d_theta;
end
```

---

## 9. Wrapper 项目中的映射关系

### 9.1 Wrapper 函数 → Field-II 底层 API

| Wrapper 函数 | 底层调用的 Field-II API | 用途 |
|-------------|------------------------|------|
| `xdc_single_rect` | `xdc_rectangles` | 创建单个矩形阵元 |
| `xdc_concatenate` | `xdc_rectangles` + `xdc_get` | 合并多个孔径 |
| `xdc_matrix_array` | `xdc_linear_array` + `xdc_linear_multirow` (等效) | 二维矩阵阵列 |
| `xdc_sparse_array_rect` | `xdc_single_rect` + `xdc_concatenate` | 稀疏阵列 |
| `xdc_rect_array` | `xdc_single_rect` + `xdc_concatenate` | 通用矩形阵列 |
| `xdc_linear_array` | `xdc_single_rect` + `xdc_concatenate` | 线性阵列封装 |
| `xdc_mills_cross` | `xdc_single_rect` + `xdc_concatenate` | Mills Cross 阵列 |
| `xdc_display_aper` | `xdc_get` | 孔径可视化 |
| `calc_rf` | `xdc_impulse`, `xdc_excitation`, `xdc_apodization`, `xdc_times_focus`, `calc_scat_multi` | 完整仿真流程 |
| `f_rf_comb` | (数据处理) | RF 数据对齐合并 |
| `gen_sphere_solid` | (几何生成) | 3D 球体散射体集 |

### 9.2 项目已使用的 Field-II API 总览

**初始化类**:
- `field_init`
- `set_field('fs', ...)`, `set_field('c', ...)`

**换能器定义类**:
- `xdc_rectangles` — 底层矩形定义 (被 xdc_single_rect 调用)
- `xdc_linear_array` — 线性阵列
- `xdc_piston` — 活塞换能器

**换能器操控类**:
- `xdc_impulse` — 设置脉冲响应
- `xdc_excitation` — 设置激励
- `xdc_apodization` — 设置变迹
- `xdc_times_focus` — 设置延迟
- `xdc_get` — 获取孔径数据
- `xdc_show` — 显示孔径信息 (通过 xdc_display_aper)
- `xdc_free` — 释放孔径

**场计算类**:
- `calc_scat_multi` — 多通道散射计算 (主要使用)
- `calc_scat` — 单通道散射计算

**未在项目中直接使用但可用的 Field-II API**:
- `xdc_triangles` — 三角形定义 (旋转阵列需要)
- `xdc_lines` — 线段边界
- `xdc_convex_array` / `xdc_convex_focused_array` — 凸面阵列
- `xdc_focused_array` / `xdc_focused_multirow` — 仰角聚焦阵列
- `xdc_convex_focused_multirow` — 凸面聚焦多排阵列
- `xdc_2d_array` — 二维稀疏阵列
- `xdc_concave` — 凹面换能器
- `xdc_baffle` — 障板条件
- `xdc_dynamic_focus` — 动态聚焦
- `xdc_center_focus` — 聚焦参考中心
- `xdc_convert` — 矩形转三角形
- `xdc_quantization` — 延迟量化
- `ele_apodization` / `ele_delay` / `ele_waveform` — 阵元级操控
- `calc_h` / `calc_hp` / `calc_hhp` — 声场计算
- `calc_scat_all` — 全合成孔径
- `field_info` / `field_debug` / `field_end` — 系统命令

---

## 10. 关键数据结构

### 10.1 坐标系

- **X轴**: 阵元宽度方向 (线性阵列的横向)
- **Y轴**: 阵元高度方向 (仰角方向)
- **Z轴**: 声传播方向 (垂直于换能器面)
- 阵列默认以原点为中心，除非指定 `center` 参数

### 10.2 物理阵元 vs 数学阵元

```
物理阵元 (Physical Element): 真实的换能器阵元 (如 64 个)
    └── 数学阵元 (Mathematical Element): 计算细分 (如 5×10 个子阵元)
```

- 数学细分由 `n_sub_x` 和 `n_sub_y` 控制
- 物理阵元编号: 从最负 xy 开始，先 x 后 y
- 数学阵元编号: 在每个物理阵元内部，同样从最负 xy 开始

### 10.3 聚焦延迟计算

Field-II 内部使用以下公式 (设置 `xdc_center_focus` 后):

```
ti = (|(xf,yf,zf) - (xc,yc,zc)| - |(xi,yi,zi) - (xc,yc,zc)|) / c
```

- `(xf, yf, zf)`: 焦点位置
- `(xc, yc, zc)`: 参考中心点 (由 xdc_center_focus 设定，初始为 (0,0,0))
- `(xi, yi, zi)`: 第 i 个物理阵元的中心
- `c`: 声速
- `ti`: 第 i 个阵元的延迟时间

### 10.4 衰减模型

- 频率相关衰减: 以 `att_f0` 为中心线性化，使该频率处衰减为 0 dB
- 频率无关衰减: 在孔径上的路径差异中体现
- 频率相关衰减: 使用平均传播距离近似
- 假设最小相位

### 10.5 矩形角点排序规则

`xdc_rectangles` 要求四个角点按顺时针排列 (1→2→3→4)，且物理阵元编号必须递增。

### 10.6 变迹/延迟矩阵格式

用于 `xdc_apodization` 和 `xdc_times_focus` 的矩阵为 `[时间点数 × 物理阵元数]`:
- 每行对应一个时间区间的变迹/延迟配置
- 列数 = 物理阵元总数
- 对于多排阵列: 先填充 x 方向再填充 y 方向

```matlab
% 多排阵列的 Hanning 变迹示例
apo = reshape(hanning(no_elem_x) * ones(1, no_elem_y), 1, no_elem_x * no_elem_y);
xdc_apodization(Th, 0, apo);
```

---

## 参考文献

1. J. A. Jensen and N. B. Svendsen. Calculation of pressure fields from arbitrarily shaped, apodized, and excited ultrasound transducers. *IEEE Trans. Ultrason., Ferroelec., Freq. Contr.*, 39:262–267, 1992.
2. J. A. Jensen. Field: A program for simulating ultrasound systems. *Med. Biol. Eng. Comp.*, 10th Nordic-Baltic Conference on Biomedical Imaging, Vol. 4, Supplement 1, Part 1:351–353, 1996.
3. J. A. Jensen. A model for the propagation and scattering of ultrasound in tissue. *J. Acoust. Soc. Am.*, 89:182–191, 1991.

---

> **文档版本**: 1.0 | **生成日期**: 2026-06-04 | **基于**: Field II User Guide Release 3.30
