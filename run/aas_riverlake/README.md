# AAS 河湖流实验 —— 运行套件

CoLM 网格化河湖流的全球评估与 bifurcation / levee / dam 影响实验。
完整的论文提纲、实验设计与方法学见
[`docs/superpowers/specs/2026-08-08-riverlake-bif-levee-dam-experiment-design.md`](../../docs/superpowers/specs/2026-08-08-riverlake-bif-levee-dam-experiment-design.md)。

## 配置速览

| 项 | 值 |
|---|---|
| 陆面分辨率 | 0.5° |
| 汇流分辨率 | 0.25°（15 min unit catchment）|
| 强迫 | JRA3Q |
| Spinup | 1996–1999 共享 + 2000 各组独立 |
| 评估期 | 2001–2005 |
| 实验 | dam × levee × bif 全因子 2³ = 8 组 |

## 文件

| 文件 | 说明 |
|---|---|
| `aas_spinup.nml` | 阶段 0：地表数据 + 初始场 + 1996–1999 共享 spinup |
| `aas_E7_ALL.nml` | **8 组实验的模板**（注释最详尽，先读这份）|
| `aas_E{0..6}_*.nml` | 由 `gen_experiments.sh` 从模板生成，勿手工编辑 |
| `hist_vars_routing.nml` | 历史变量清单（E1–E7 用）：只留汇流场 + `rnof` |
| `hist_vars_land_full.nml` | 历史变量清单（E0_CTL 用）：汇流场 + 陆面水循环诊断 |
| `gen_experiments.sh` | 从模板生成 8 份 nml，并自检一致性 |
| `gen_jobs.sh` | 生成 9 份 LSF 作业脚本 |
| `check_inputs.sh` | **前置检查**：确认参数文件含 levee/bif/dam 字段 |
| `distribute_restart.sh` | 把共享 spinup 的 restart 分发给 8 组 |

## 运行流程

```bash
# ---- 0. 前置检查（几秒钟，务必先跑）----
./check_inputs.sh
#   缺字段不会在 namelist 解析时报错，而是在 river_lake_network_init 读参数时
#   才失败——那时 480 核作业已排队、启动、读完地表数据，浪费小时级机时。

# ---- 1. 生成 nml 与作业脚本 ----
./gen_experiments.sh          # 8 份 nml + 一致性自检
NCORE=480 ./gen_jobs.sh       # 9 份 LSF 脚本

# ---- 2. 共享 spinup ----
bsub < submit_SPINUP.lsf

# ---- 3. 分发 restart ----
./distribute_restart.sh        # 先预览
./distribute_restart.sh --run  # 确认后执行

# ---- 4. 单年标定（把 aas_E7_ALL.nml 的 end_year 临时改成 2000）----
bsub < submit_E7_ALL.lsf
#   E7_ALL 是最慢的一组，用它标定得到墙钟时间上界。

# ---- 5. 恢复 end_year=2005，提交全部 8 组 ----
for e in E0_CTL E1_D E2_L E3_B E4_DL E5_DB E6_LB E7_ALL; do
    bsub < submit_${e}.lsf
done
```

## H22 水库调度

原生网格河湖流可用 `DEF_Reservoir_Method = 2` 选择 H22；`1` 仍为原有 Funato 调度，`0` 关闭水库。
先用 `./check_inputs.sh /path/to/reservoir_parameters.nc 2` 检查 H22 字段。

H22 使用与现有水库目录逐行对齐的 `dam_TotalVol_mcm`、`dam_ConVol_mcm`、`dam_Qn`、`dam_Qf`，
另需 `dam_upreal`（上游集水面积，km²）。库容单位为百万 m³，流量单位为 m³/s。
常规库容 `C` 直接取原始 `dam_ConVol_mcm`，防洪专用库容 `F = TotalVol - C`；
低库容阈值为 `0.5 C`，紧急阈值为 `C + 0.95 F`，泄流系数为 `max(1 - F_mcm / dam_upreal / 0.2, 0)`。
每座坝单独计算系数；不得把主河道单元面积代作上游集水面积。
要求有效水库满足 `TotalVol >= C > 0`、`dam_upreal > 0`、`Qf >= Qn >= 0`，且参数有限。

这里移植的是 [CaMa-Flood PR #124 对应版本](https://github.com/global-hydrodynamics/CaMa-Flood_v4/blob/775853bc5e0291e2637b7ae970759a6783fdd397/src/cmf_ctrl_damout_mod.F90) 的 H22 实现。
其紧急阈值使用防洪专用库容的 `0.95`，与 [Hanazaki et al. (2022)](https://doi.org/10.1029/2021MS002944)
论文文字中的 `0.8` 不同；保留该版本原始分段公式和边界。
冷启动已建坝使用低库容阈值 `0.5 C` 作为初始库容；未建坝仍按原有年份门控处理。
重启沿用保存的库容和上游实际泄流缓存，不新增动态状态；切换调度方案会改变后续出流，不重置库容。

## 编译前提

`include/define.h` 需要：

```c
#define GRIDBASED
#define LULC_IGBP_PFT
#define USEMPI
#define GridRiverLakeFlow
#define vanGenuchten_Mualem_SOIL_MODEL

#undef BGC          /* 仓库默认是 #define，必须改 */
#undef CROP         /* 同上，BGC 关闭后须一并关闭 */
#undef TRACER
#undef CaMa_Flood
#undef URBAN_MODEL
#undef LULCC
```

改完 `define.h` **必须 `make clean` 后完整重编**——Makefile 不追踪 `define.h` 的依赖，增量编译会产生 ABI 不一致的目标文件。

## 三条容易踩的坑

1. **`DEF_hist_vars_out_default` 不是有效变量。** 它出现在 `run/Amazon.nml:52` 等旧 nml 里，但 `share/MOD_Namelist.F90` 从未定义它。本套件不使用它。真正的机制是：`DEF_hist_vars` 各字段默认全为 `.true.`，由 `DEF_HIST_vars_namelist` 指向的文件覆盖（`MOD_Namelist.F90:2133-2142`），**漏列的字段会以 `.true.` 生效**。所以两份变量清单都完整列出了全部 482 个字段。

2. **history nml 里的变量名必须全部有效。** 读取失败会直接 `CoLM_Stop`（`MOD_Namelist.F90:2138-2142`），不是警告。

3. **开启 BIF 会强制全局 dt 同步。** `MOD_Grid_RiverLakeFlow.F90:974` 显示 `DEF_USE_BIFURCATION = .true.` 时走 `sync_global_routing_dt`（全部 worker），关闭时只在河系通信子 `p_comm_rivsys` 内 allreduce（`:980`）。因为分汊路径跨河系配对供受水单元，单一全局 dt 才能保证配对体积转移守恒。所以 E3/E5/E6/E7 四组会明显慢于另外四组——这是物理必然，不是 bug，并且是论文 Section 4.5 要量化的结果之一。
