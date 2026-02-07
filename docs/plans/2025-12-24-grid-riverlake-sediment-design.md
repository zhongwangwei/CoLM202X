# GridRiverLakeSediment 模块设计文档

**日期：** 2025-12-24
**目标：** 将 CaMa-Flood 泥沙模块移植到 GridRiverLakeFlow 框架中

---

## 1. 概述

### 1.1 背景

泥沙模块已在 CaMa-Flood (`extends/CaMa/src/cmf_ctrl_sed_mod.F90`) 中实现，需要移植到 GridRiverLakeFlow 框架中，使其能够在定义 GridRiverLakeFlow 时使用并输出相关结果。

### 1.2 数据来源

静态数据已保存在 `DEF_UnitCatchment_file` (如 `grid_routing_data.nc`) 中：
- `sed_frc(nseqmax, sed_n)` - 泥沙粒径分布比例，3个粒径级别
- `sed_slope(nseqmax, slope_layers)` - 坡面坡度，10层

### 1.3 设计决策

| 决策项 | 选择 |
|--------|------|
| 移植范围 | 完整移植（悬沙、推移质、侵蚀产沙、沉积层交换） |
| 产沙驱动 | 降水驱动（类似 CaMa-Flood） |
| 时间步长 | 独立子循环（与水流时间步分离） |
| 输出变量 | 完整变量（包括诊断变量） |
| 代码组织 | 独立模块 `MOD_Grid_RiverLakeSediment.F90` |

---

## 2. 模块架构

### 2.1 新建文件

**`main/HYDRO/MOD_Grid_RiverLakeSediment.F90`**

```
MOD_Grid_RiverLakeSediment
├── 参数常量
│   ├── nsed (粒径级数，从NC文件读取)
│   ├── totlyrnum (沉积层数)
│   ├── 物理常数 (lambda, psedD, pwatD, visKin, vonKar等)
│   └── 产沙参数 (pyld, pyldc, pyldpc, dsylunit)
│
├── 静态数据 (从 DEF_UnitCatchment_file 读取)
│   ├── sed_frc(:,:)     - 泥沙粒径分布比例 [numucat, nsed]
│   ├── sed_slope(:,:)   - 坡面坡度 [numucat, nlfp]
│   └── sDiam(:)         - 各粒径直径 [nsed]
│
├── 状态变量
│   ├── sedcon(:,:)      - 悬沙浓度 [numucat, nsed]
│   ├── layer(:,:)       - 活动层泥沙量 [numucat, nsed]
│   └── seddep(:,:,:)    - 沉积层泥沙量 [numucat, totlyrnum, nsed]
│
├── 诊断变量
│   ├── sedout(:,:)      - 悬沙通量
│   ├── bedout(:,:)      - 推移质通量
│   ├── netflw(:,:)      - 净交换通量
│   ├── sedinp(:,:)      - 侵蚀产沙输入
│   └── shearvel(:), critshearvel(:,:), setvel(:)等
│
└── 累积输出变量 (用于历史输出)
```

### 2.2 修改文件清单

| 文件 | 修改内容 |
|------|----------|
| `include/define.h` | 添加 `#define GridRiverLakeSediment` 预处理开关 |
| `main/MOD_Namelist.F90` | 添加泥沙相关 namelist 参数 |
| `main/HYDRO/MOD_Grid_RiverLakeNetwork.F90` | 添加 sed_frc, sed_slope 数据读取 |
| `main/HYDRO/MOD_Grid_RiverLakeFlow.F90` | 添加泥沙模块调用点 |
| `main/HYDRO/MOD_Grid_RiverLakeTimeVars.F90` | 添加泥沙状态变量的 restart 读写 |
| `main/HYDRO/MOD_Grid_RiverLakeHist.F90` | 添加泥沙输出变量与写出逻辑 |
| `Makefile` 或 CMake 配置 | 添加新模块编译 |

---

## 3. 初始化流程

```fortran
SUBROUTINE grid_sediment_init()
   ! 1. 读取泥沙参数维度
   !    - 从 DEF_UnitCatchment_file 读取 nsed (sed_n 维度)
   !    - 从 DEF_UnitCatchment_file 读取 nlfp (slope_layers 维度)

   ! 2. 读取静态数据 (使用 readin_riverlake_parameter 模式)
   !    - sed_frc(numucat, nsed) ← 'sed_frc'
   !    - sed_slope(numucat, nlfp) ← 'sed_slope'

   ! 3. 设置粒径直径 (通过 namelist 或硬编码)
   !    - sDiam(:) = [0.0002, 0.002, 0.02] ! 示例：粉砂、细砂、粗砂

   ! 4. 计算沉降速度
   !    - setvel(:) = calc_settling_velocity(sDiam, visKin, psedD, pwatD)

   ! 5. 分配状态变量数组
   !    - sedcon(numucat, nsed)
   !    - layer(numucat, nsed)
   !    - seddep(numucat, totlyrnum, nsed)

   ! 6. 初始化状态变量
   !    - 如有 restart 文件则读取
   !    - 否则根据 sed_frc 和活动层厚度初始化：
   !      layer(i,:) = lyrdph * topo_rivwth(i) * topo_rivlen(i) * sed_frc(i,:)
END SUBROUTINE
```

### 3.1 Namelist 配置

```fortran
! 在 MOD_Namelist 中添加
logical :: DEF_USE_SEDIMENT = .false.    ! 是否启用泥沙模块
real(r8) :: DEF_SED_LAMBDA  = 0.4        ! 孔隙率
real(r8) :: DEF_SED_LYRDPH  = 0.00005    ! 活动层厚度 [m]
real(r8) :: DEF_SED_DENSITY = 2.65       ! 泥沙密度 [g/cm³]
integer  :: DEF_SED_TOTLYRNUM = 5        ! 沉积层数
character(len=256) :: DEF_SED_DIAMETER = "0.0002,0.002,0.02"  ! 粒径 [m]

! 产沙参数
real(r8) :: DEF_SED_PYLD    = 0.01       ! 产沙系数
real(r8) :: DEF_SED_PYLDC   = 2.0        ! 坡度指数
real(r8) :: DEF_SED_PYLDPC  = 2.0        ! 降水指数
```

---

## 4. 核心计算流程

### 4.1 主计算子程序

```fortran
SUBROUTINE grid_sediment_calc(deltime)
   ! 输入: deltime - 水流计算的总时间步长

   ! 1. 累积水流诊断变量
   !    - 在每个水流子步调用 sediment_diag_accumulate()
   !    - 累积 veloc_riv, wdsrf_ucat 用于后续平均

   ! 2. 泥沙独立子循环
   DO WHILE (sed_time_remaining > 0)

      dt_sed = min(sed_time_remaining, DEF_SED_DT_MAX)

      ! 2.1 计算平均水流变量
      avg_veloc = acc_veloc / acc_time
      avg_wdsrf = acc_wdsrf / acc_time

      ! 2.2 计算泥沙输运参数
      CALL calc_sediment_params(avg_veloc, avg_wdsrf, &
           shearvel, critshearvel, susvel)

      ! 2.3 悬沙与推移质平流输运
      CALL calc_sediment_advection(dt_sed, &
           sedcon, layer, sedout, bedout)

      ! 2.4 侵蚀产沙 (降水驱动)
      CALL calc_sediment_yield(precip, sed_slope, sedinp)

      ! 2.5 悬浮-沉积交换
      CALL calc_sediment_exchange(dt_sed, &
           sedcon, layer, netflw)

      ! 2.6 沉积层重分布
      CALL calc_layer_redistribution(layer, seddep)

      ! 2.7 更新累积输出变量
      CALL sediment_diag_accumulate_output(dt_sed)

      sed_time_remaining = sed_time_remaining - dt_sed
   ENDDO
END SUBROUTINE
```

### 4.2 关键物理计算函数

| 函数 | 功能 | 公式来源 |
|------|------|----------|
| `calc_settling_velocity` | 沉降速度 | Stokes-Rubey 公式 |
| `calc_shear_velocity` | 剪切流速 | Manning 公式 |
| `calc_critical_shear` | 临界剪切力 | Shields / Egiazoroff 方程 |
| `calc_suspend_velocity` | 悬浮速度 | Uchida & Fukuoka (2019) Eq.44 |

### 4.3 上下游通信

使用现有的 worker_push_data 机制：
- `worker_push_data(push_next2ucat, ...)` - 获取下游泥沙浓度
- `worker_push_data(push_ups2ucat, ...)` - 汇总上游泥沙通量

---

## 5. 模块集成

### 5.1 在 MOD_Grid_RiverLakeFlow.F90 中的调用点

```fortran
SUBROUTINE grid_riverlake_flow(year, deltime)
   ...
   ! 在水流计算主循环中添加泥沙诊断累积
   DO WHILE (any(dt_res > 0))
      ...
      ! 现有水流计算
      ...

#ifdef GridRiverLakeSediment
      ! 每个水流子步累积泥沙所需的水流变量
      IF (DEF_USE_SEDIMENT) THEN
         CALL sediment_diag_accumulate(dt_all, veloc_riv, wdsrf_ucat)
      ENDIF
#endif
      ...
   ENDDO

#ifdef GridRiverLakeSediment
   ! 水流计算完成后，执行泥沙计算
   IF (DEF_USE_SEDIMENT) THEN
      CALL grid_sediment_calc(acctime_rnof)
   ENDIF
#endif
   ...
END SUBROUTINE
```

### 5.2 降水数据传递

```fortran
SUBROUTINE grid_sediment_forcing_put(precip_uc)
   ! precip_uc: 单元集水区的降水量 [numucat]
   sed_precip(:) = precip_uc(:)
END SUBROUTINE
```

### 5.3 静态数据读取

在 `MOD_Grid_RiverLakeNetwork.F90` 的 `build_riverlake_network()` 中添加：

```fortran
#ifdef GridRiverLakeSediment
   IF (DEF_USE_SEDIMENT) THEN
      CALL readin_riverlake_parameter(parafile, 'sed_frc',   rdata2d=sed_frc)
      CALL readin_riverlake_parameter(parafile, 'sed_slope', rdata2d=sed_slope)
   ENDIF
#endif
```

### 5.4 Restart 支持

在 `MOD_Grid_RiverLakeTimeVars.F90` 中添加：

```fortran
! READ
CALL vector_read_and_scatter(file_restart, sedcon, numucat, 'sedcon', ucat_data_address)
CALL vector_read_and_scatter(file_restart, layer,  numucat, 'layer',  ucat_data_address)
! seddep 需要特殊处理 (3D数组)

! WRITE
CALL vector_gather_and_write(sedcon, numucat, totalnumucat, ucat_data_address, &
   file_restart, 'sedcon', 'ucatch')
```

---

## 6. 历史输出系统

### 6.1 累积变量

```fortran
! 在 MOD_Grid_RiverLakeHist.F90 中添加
real(r8), allocatable :: a_sedcon     (:,:)   ! 悬沙浓度 [numucat, nsed]
real(r8), allocatable :: a_sedout     (:,:)   ! 悬沙通量 [numucat, nsed]
real(r8), allocatable :: a_bedout     (:,:)   ! 推移质通量 [numucat, nsed]
real(r8), allocatable :: a_sedinp     (:,:)   ! 侵蚀产沙量 [numucat, nsed]
real(r8), allocatable :: a_netflw     (:,:)   ! 净交换通量 [numucat, nsed]
real(r8), allocatable :: a_layer      (:,:)   ! 活动层存量 [numucat, nsed]
real(r8), allocatable :: a_shearvel   (:)     ! 剪切流速 [numucat]
real(r8), allocatable :: a_critshear  (:,:)   ! 临界剪切力 [numucat, nsed]
```

### 6.2 输出文件结构

```
历史文件 (*_unitcat_*.nc):
├── 维度
│   ├── time
│   ├── lon_ucat, lat_ucat
│   └── sed_class (新增，大小=nsed)
│
├── 坐标变量
│   └── sed_diameter(sed_class)  ! 各粒径直径 [m]
│
└── 泥沙变量 (均为 [lon_ucat, lat_ucat, sed_class, time])
    ├── f_sedcon      - 悬沙浓度 [m³/m³]
    ├── f_sedout      - 悬沙通量 [m³/s]
    ├── f_bedout      - 推移质通量 [m³/s]
    ├── f_sedinp      - 侵蚀产沙 [m³/s]
    ├── f_netflw      - 净交换通量 [m³/s]
    ├── f_layer       - 活动层存量 [m³]
    ├── f_shearvel    - 剪切流速 [m/s] (无sed_class维)
    └── f_critshear   - 临界剪切力 [m/s]
```

### 6.3 输出控制 Namelist

```fortran
! 添加到 DEF_hist_vars
logical :: sedcon    = .true.   ! 悬沙浓度
logical :: sedout    = .true.   ! 悬沙通量
logical :: bedout    = .true.   ! 推移质通量
logical :: sedinp    = .true.   ! 侵蚀产沙
logical :: netflw    = .true.   ! 净交换通量
logical :: sedlayer  = .true.   ! 活动层存量
logical :: shearvel  = .false.  ! 剪切流速 (诊断)
logical :: critshear = .false.  ! 临界剪切力 (诊断)
```

---

## 7. 实施步骤

### 阶段一：框架搭建
1. 创建 `MOD_Grid_RiverLakeSediment.F90` 模块框架
2. 添加预处理开关 `GridRiverLakeSediment` 到 `define.h`
3. 添加 namelist 参数到 `MOD_Namelist.F90`
4. 实现静态数据读取（sed_frc, sed_slope）

### 阶段二：核心计算
1. 移植物理计算函数
   - `calc_settling_velocity`
   - `calc_shear_velocity`
   - `calc_critical_shear`
   - `calc_suspend_velocity`
2. 实现平流输运计算 `calc_sediment_advection`
3. 实现悬浮-沉积交换 `calc_sediment_exchange`
4. 实现产沙计算 `calc_sediment_yield`
5. 实现沉积层重分布 `calc_layer_redistribution`

### 阶段三：集成测试
1. 与 `MOD_Grid_RiverLakeFlow` 集成
2. 添加历史输出到 `MOD_Grid_RiverLakeHist`
3. 添加 restart 支持到 `MOD_Grid_RiverLakeTimeVars`
4. 更新 Makefile/CMake 编译配置

### 阶段四：验证
1. 单点测试 - 验证物理计算正确性
2. 区域测试 - 与 CaMa-Flood 泥沙结果对比
3. 质量守恒检验

---

## 8. 参考资料

- CaMa-Flood 泥沙模块源码: `extends/CaMa/src/cmf_ctrl_sed_mod.F90`
- Uchida & Fukuoka (2019) - 悬浮速度公式
- Egiazoroff 方程 - 混合粒径临界剪切力
- Shields 曲线 - 单粒径临界剪切力
