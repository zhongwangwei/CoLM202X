# 泥沙/溶质审计（基线 b100745f）修复记录

- 分支 `feature/constituent-transport`；审计基线 `b100745f`；工作 HEAD 见文末（提交后更新）。
- 审计包 `CoLM-sediment-solute-audit-b100745.zip` 已解压；7 项确认（N06/N07 为上一轮 F02 转发修复引入）。
- 本轮修复：**N01、N02、N03、N04、N06**；**N07 部分**；N05 未实施。
- 未运行真实 MPI / 完整模型 / 真实 NetCDF 集成；以下验证为源码 + 生产例程提取/宿主接口级。

## N02 冠层净截留误当总截留 —— 已修复并数值验证

- 宿主 `qintr_rain/snow` 是净截留（P − direct throughfall − drainage）；材料层需要
  gross interception（P − direct throughfall）。修复（`#ifdef MaterialTransport`，非 MT 路径不变）：
  - `MOD_LeafInterception.F90` 新增 MT-only 模块状态 `mt_leaf_qdrp_rain/snow`，在
    `LEAF_interception_CoLM2014` 主分支置为真实 `tex/deltim`（无降水释放分支置 0）；
    PFT wrapper 按 `pftfrac` 加权聚合；CoLM2024 经 2014 自动继承；
  - `CoLMMAIN.F90` MT 块：
    `gross_intr_rain = max(max(qintr_rain,0)+max(mt_leaf_qdrp_rain,0),0)`（snow 同）。
- 验证：用审计包 `repro_canopy_gross_interface.py` 的**原样宿主例程 + 原样 mt_precip**（仅把
  define.h 打开 `MaterialTransport` 并补 import），CoLM2014 方案，初始冠层水 0.05 mm、旧 CL
  0.0001 kg/m²、1800 s 输入 1 mm 纯雨：

```
host pure-gap throughfall          = 9.009069189841E-01 mm
host post-mix canopy drainage      = 4.908317171802E-02 mm
caller-supplied gross interception = 9.909308101592E-02 mm
host actual gross interception     = 9.909308101592E-02 mm   (修复前为 5.00099E-02)
transported old canopy CL          = 3.292115997843E-05 kg/m2  (修复前 0)
retained old canopy CL             = 6.707884002157E-05 kg/m2
```

  `actual_ground == corrected_ground`（误差 0），且 `tex_rain>0`；宿主水量与材料质量由同一路径闭合。
- 编译：`gfortran -cpp -E`（`#define MaterialTransport`）对 `MOD_LeafInterception.F90`、
  `CoLMMAIN.F90` rc=0。
- 未做：全模型/真实数据验证；durable repo 回归测试（当前用审计包 harness 适配版，正式测试待补）。

## N06 干河段转发 MPI 顺序 —— 已修复

- 位置 `MOD_MT_RiverLake.F90` 转发 push：原按本 rank `trc_dry_fwd` 是否为零 `CYCLE`，与
  `worker_push_data` 的固定点对点拓扑不匹配可挂起。
- 修复：每个 generic tracer **固定顺序、无条件**调用一次 push；本地无物质时发送零向量，
  接收端按拓扑 irecv/isend 匹配；顺序与数据无关。
- 验证：`tests/material_transport/check_audit_p1_fixes.py` 断言“无数据相关 CYCLE、零值发送、
  push 在每 tracer 循环内” PASS；预处理 rc=0。
- 未做：真实多 worker/跨进程河系/局部零库存 MPI 死锁测试（需 MPI 运行环境授权）。

## N07 干河段转发忽略真实流向 —— 部分修复

- 修复：`mt_rlc_diag_accumulate` 把本步 `hflux_fc` 传入 `mt_riverlake_diag_accumulate`；只有
  `hflux_fc(i) >= 0`（实际流向静态下游）才走下游 push。反向/方向未知时，真实质量进入
  `trc_dry_pending`：不送错方向、不删除；下一子步回 `trc_mass` 参与输运，并在
  `mt_rlc_end_routing` 前统一释放回库存（`mt_riverlake_release_pending_drain`，无 restart schema 变化）。
- 未完成：反向流下应把残质送到真实上游接收端（多出边按通量分配）。现有 push 映射只有下游方向，
  需新增反向拓扑 pushdata 或在 coupling 层按 `ucat_ups`+负通量构造；并做真实 MPI 验证。
  当前实现避免了旧错误方向与删除，但不是最终物理投递，故记为**部分修复**。

## N03 洪水 late pool 先溶解/析出再分配 —— 已修复并数值验证

- `MOD_MT_SoilWater.F90` late flood 分支：在 residue 加入后、`late_ratio`/入渗/积水/溢流分配**之前**，
  对有限溶解度溶质调用 `mt_equilibrate_dissolved(itrc, late_water, late_tracer, trc_surface_solid)`；
  干地表固相可随 late 水溶出，饱和洪水多余部分留在同位置固相，浓度不超过上限。
- 验证（审计包 `repro_solute_land.py` 实际 late 混合/入渗赋值 + 完整平衡例程，修后适配断言）：
  - 0.1 kg/m² 固态 CL + 0.1 mm 纯洪水：入渗 CL = **0.0218**（=0.218×0.1），残固相 **0.0782**，
    浓度 **0.218 kg/kg**（上限）；
  - 饱和洪水 0.0218 kg/m²、90% 载体蒸发剩 0.01 mm：入渗 **0.00218**，残固相 **0.01962**，
    浓度 **0.218 kg/kg**；与修复前 0 和 2.18 kg/kg 两个反例对照。
- 未做：河陆边界账/完整耦合运行；durable 测试目前为断言+数值不变量（`check_audit_p1_fixes.py`）。

## N04 冰载体消失后非挥发溶质可回收 —— 已修复并数值验证

- `MOD_MT_Evapo.F90`：每个 nonvolatile solute 在 `mt_evapo` 入口检查各层；`trc_wice>tiny` 且
  `wice_soisno<=trc_water_min_for_ratio`（冰载体已空）时，按同层液水存在与否转入
  `trc_wliq_soisno` 或 `trc_solid_soisno`，清空 ice 库存，并立即与当前液水做溶解度平衡；
  `MOD_MT_SoilWater.F90` 末尾有同一 fallback（防其它路径漏网），ice 库存清空保证不重复转移。
- 验证（审计包完整 Evapo 例程）：1 mm 土壤冰含 0.001 kg/m² CL，完全升华后连续 5 步湿、无冰状态 →
  liquid=**0.001**、ice=**0**、solid=0（修复前 ice 一直 0.001、liquid=0）；
  无液水时转入 solid=0.001（可复湿），ice=0；库存守恒、不重复释放。
- 未做：雪层/混合相分支的完整集成对照；真实模型验证待授权。

## N01 泥沙双向粗时间步破坏恒浓度 —— 已修复

- `MOD_MT_Sediment.F90` `calc_sediment_advection`：普通（非 BIF/levee/joint）路径按
  donor gross outflow 水量做材料子步：`nsub = ceil(max_i(rivout_abs(i)*dt/rivsto_donor(i)))`，
  上限 240，随后以 `dt/nsub` 递归调用，每个子步重新计算 donor 浓度；BIF/levee/joint 路径
  保持原有共享初始库存算法不变。水动力步长、宿主水量、morphology 间隔均未改。
- 验证（审计三河段反例，10 m3/段、C=0.001、dt=10 s、每界面正/反向各 1 m3/s）：
  - 单步旧算法复现 [0, 0.002, 0.001]；
  - 本修复给出 nsub=2（2 m3/s×10 s / 10 m3），恒浓度保持 **0.001/0.001/0.001**，
    总质量 0.03 守恒；
  - 结构断言：子步按 gross outflow 触发、上限 240、仅普通路径、每子步重算 donor 浓度。
- 未做：真实河网/MPI 集成与性能实测（子步上限保证最坏 ~240× 平流计算，需在 7920 记录耗时）；
  床层/推移质在同一子步逻辑内一并推进，未单独改变其参数化。

## N05 —— 未实施

- N01：`MOD_MT_Sediment.F90` 双向粗时间步 donor 限流；需按 gross 出流/donor 水量控制材料子步或
  守恒双向时间积分；`repro_sediment_bidirectional.py` 可直接作为回归（修复后改断言）。
- N03：`MOD_MT_SoilWater.F90` late-pool 洪水加入后先溶解度平衡再分配；两个复现在审计包。
- N04：`MOD_MT_SoilWater.F90`/`MOD_MT_Evapo.F90` 冰完全升华后非挥发溶质转可复湿残留固相。
- N05：`MOD_HistAccRestart.F90`/`MOD_Grid_RiverLakeHist.F90` 材料累计保存/恢复策略与部分窗口语义。
- 这些项目当前保留为未完成，不计入通过。

## 本地检查命令

```sh
python3 tests/material_transport/check_audit_p1_fixes.py          # 0 failure（含 N06/N07 断言）
PATH=/opt/homebrew/bin:$PATH python3 tests/material_transport/check_mt_frac_physics.py   # 15 PASS
PATH=/opt/homebrew/bin:$PATH python3 tests/material_transport/check_rnof_handoff.py      # PASS
# gfortran -cpp -E（define.h 打开 MaterialTransport）对全部改动文件 rc=0
```
