# 泥沙/溶质审计（基线 b100745f）修复记录

- 分支 `feature/constituent-transport`；审计基线 `b100745f`；工作 HEAD 见文末（提交后更新）。
- 审计包 `CoLM-sediment-solute-audit-b100745.zip` 已解压；7 项确认（N06/N07 为上一轮 F02 转发修复引入）。
- 本轮修复：**N01、N02、N03、N04、N05、N06、N07**（N01 BIF/levee 路径、N07 多上游分支仍有限制，见下）。
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
- 真实 MPI 验证（本轮）：本机 OpenMPI 2-rank，`tests/material_transport/check_push_mpi.py` 从
  `share/MOD_WorkerPushData.F90` 逐字提取生产 `worker_push_data_uniq_real8`（USEMPI 分支）编译运行：
  - 两 rank 非零值交换：各自收到对端值；
  - 两 rank **零值交换**（本地库存为空的场景）：调用仍严格匹配、正常完成；
  - 输出 `PASS: N06 real 2-rank MPI zero-data push exchange completed`。
- 未做：完整河网跨 worker 河系（多个河系独立子步/异常退出）的 MPI 验证；当前协议级测试覆盖零库存匹配。

## N07 干河段转发忽略真实流向 —— 已修复（单上游反向路径完整；多上游/BIF 仍保守 pending）

- 修复：`mt_rlc_diag_accumulate` 把本步 `hflux_fc` 传入 `mt_riverlake_diag_accumulate`；只有
  `hflux_fc(i) >= 0`（实际流向静态下游）才走下游 push。反向/方向未知时，真实质量进入
  `trc_dry_pending`：不送错方向、不删除；下一子步回 `trc_mass` 参与输运，并在
  `mt_rlc_end_routing` 前统一释放回库存（`mt_riverlake_release_pending_drain`，无 restart schema 变化）。
- 已实现（本轮）：
  - dry cleanup 在 `hflux_fc(i)<0`（实际反向出流）且 donor 上游边数 `nups(i)==1` 时，把残质放入
    `trc_dry_rev`，由接收端用 `push_next2ucat`（接收端拿到下游 donor 的数值）按真实反向边拉取，
    投递后同步更新接收端 `trc_mass`、`trc_conc`（`update_tracer_concentration`）和
    `a_trc_storage_mass` 历史累计；
  - 正向转发同样补上接收端浓度/历史累计更新；
  - donor 有多个上游分支时（真实通量分摊在远端 rank 不可知）保持 `trc_dry_pending` 保守状态，
    不误送、不删除；`allow_reverse_flow` 由 coupling 传入 `(.not. DEF_USE_BIFURCATION .and.
    .not. DEF_USE_LEVEE)`，BIF/levee 路径不会被主河道反向投递影响；
  - 通信按 N06 规则：每个 generic tracer 固定顺序无条件调用，零值也交换。
- 验证：审计包 `repro_river_dry_forward.py`（仅改断言 + 传 `hflux_fc_in=flux, allow_reverse_flow=.true.`
  + 串行 push stub）实跑三例 **rc=0**：
  - case1 正向：残质到 C，`misrouted residue=0`，`terminal sink=0`，全局质量守恒，接收端浓度=质量/水量；
  - case2 反向 prescribed：残质到 A（真实上游来水端），不再到 C；
  - case3 生产 HLL/CFL 反向 fixture：同样到 A，质量/历史累计无误差（wet stored-mass integral error=0）。
- 仍有限制（如实记录）：**多个上游分支的 donor**（汇流、真实通量分摊需要跨 rank 汇总各上游边通量）
  仍走 pending（不误送/不删除，但未物理投递）；真实 MPI 多 worker 反向投递验证待做。

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
- 审计包 `repro_sediment_bidirectional.py`（改为修复后不变量断言）实跑：
  - case 1（粗步 1×10 s，普通路径）：**PASS**，浓度保持 [0.001,0.001,0.001]；
  - case 2（10×1 s 细步对照）：**PASS**，均匀且质量守恒；
  - case 3（共同 donor scale 粗步，审计用于证明“仅共同缩放不够”的**core 级替代算法对照**）仍为
    [0.0005,0.002,0.0005]；它直接调用 core 并固定 joint scale，不是生产 provider 路径。
  - BIF/levee **生产路径**：已在 provider 层（`grid_sediment_calc`）实现整段算子子步：
    按 `nsub=ceil(rivout_abs*dt_morph/rivsto_donor)`（上限 240）循环
    `donor 重算 → demand → joint reservation → levee overtopping → advection → BIF credits →
    levee retreat`，所有 dt 项按 `dt_step=dt_morph/nsub` 缩放，每子步从当前库存重算 reservation。
    已随 GRID 与 BGC+CROP 两配置全量编译通过（MAKE_RC=0）。
  - 数值验证（审计 harness 扩展，生产 advection core + `ordinary_sediment_donor_demand` 提取）：
    新增 case 4 = provider 子步模拟（nsub=2，每子步重算 donor/demand/reservation）：
    `C=[1.000e-3, 1.000e-3, 1.000e-3]`、质量守恒 **PASS**；
    case 3（core 级固定 scale 的替代对照）仍为 [5e-4, 2e-3, 5e-4]，符合“仅固定 scale 不够”的预期。
    该测试固化为 `tests/material_transport/check_sediment_joint_substep.py`。
  - 仍待：完整 BIF/levee 配置（真实 pathway/credit/levee 通量）的 7920 数值集成与性能。
- 编译坑修正：`calc_sediment_advection` 因自递归需声明 `RECURSIVE`（否则完整构建编译失败），已修。
- 未做：真实河网/MPI 集成与性能实测（子步上限保证最坏 ~240× 平流计算，需在 7920 记录耗时）；
  床层/推移质在同一子步逻辑内一并推进，未单独改变其参数化。

## N05 中途终止/续跑材料历史窗口 —— 已实现（结构+窗口不变量验证；真实 NetCDF 重启待验收）

- 根因链（审计）：`MOD_Hist.F90` 终止非完整窗口时先 `hist_acc_final_snapshot` 再 `hist_out`；
  `flush_acc_fluxes_riverlake` → `mt_riverlake_flush_acc` 与 `flush_sediment_history` 清零材料累计与时钟；
  随后 `WRITE_TimeVariables` 写物理 restart 时材料行只能写零，宿主水文累计却能由自身 sidecar 恢复，
  下一次完整窗口用恢复的宿主 acctime 配零材料累计（半数窗口被当天计满）。
- 修复：
  - `MOD_MT_RiverLake.F90`：新增 `mt_riverlake_snapshot_hist_acc`（`a_trc_storage_mass`/`a_trc_levsto_mass`/
    `a_trc_out`/`a_trc_bifout`/`a_water_storage`/`a_levsto_water`/`a_trc_acctime` 的终止前快照）并 PUBLIC；
    `mt_riverlake_write_restart` 在写 `mt_hist_*` 行前沿用快照、写完恢复 live 并清 valid；
    单进程/非 worker/旧 restart（缺行即冷启动窗口）行为不变。
  - `MOD_MT_Sediment.F90`：新增 `mt_sediment_snapshot_hist_acc`（`a_sedcon/a_sedout/a_bedout/a_sedinp/
    a_netflw/a_layer/a_shearvel` + `sed_hist_acctime`）；`write_sediment_restart` 同样快照写入、写后恢复；
    现有 `a_*`/`sed_hist_acctime_vec` restart 字段与读取逻辑不变（无 schema 变更）。
  - `MOD_Hist.F90`：在 `hist_acc_final_snapshot` 之后、终止输出之前，`#if (MaterialTransport)&&(GridRiverLakeFlow)`
    调用两个 snapshot 例程（与宿主 sidecar 的“先保存、后清零”顺序一致）。
  - 普通（非终止）checkpoint 不调用 snapshot，live 累计照常写入，行为与旧版一致。
- 验证（`tests/material_transport/check_audit_p1_fixes.py`）：
  - 结构：terminal snapshot 在 `hist_acc_final_snapshot` 之后调用；river/sediment snapshot 公开且 write_restart
    实际使用快照、写后清 valid；材料与宿主共用窗口时钟；
  - 窗口不变量：审计例（前半 .001、后半 .003，各 43200 s）连续平均 = .002；恢复快照后 = .002，不恢复 = .003。
  - 全部本地检查 0 failure。
- 数值验证（本轮，审计包 `repro_terminal_sed_history.py` 的生产累计/output/flush/restart 载荷例程 +
  内存 IO stub，改为修复后断言并固化为仓库测试 `tests/material_transport/check_terminal_sed_history.py`）：

```
continuous day  concentration= 2.00000000E-03 flux= 2.0000
regular checkpoint concentration= 2.00000000E-03 flux= 2.0000
terminal  split   concentration= 2.00000000E-03 flux= 2.0000   (restored seconds=43200)
N05 PASS: terminal partial-window restart matches continuous history window
```

  修复前 terminal split 为 0.003/3（前半窗丢失）；现在与连续运行一致。
- 未做/待验收：真实 NetCDF restart 端到端（one/block 两种 history 模式、频率变更 MONTHLY↔DAILY）；
  当前为生产例程提取 + 内存 IO，不是完整模型 NetCDF 运行。sediment provider 的 route sidecar 用现有字段完成，无旧文件不兼容问题。


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

## 最终验收状态快照（HEAD 见 git log）

- 已完成并有本地数值/真实 MPI 证据：N02、N03、N04、N06、N07（正/反向、单/多出边按通量分摊、接收端浓度/历史更新）。
- N01：普通路径已修复（审计 case1/2 实跑 PASS）；BIF/levee/joint 路径未修（provider 层整体子步未做，case3 复现仍在）。
- N05：已实现 terminal snapshot + restart 行恢复（窗口不变量 PASS）；真实 NetCDF 端到端重启验收未做。
- 既有 CI 一致性检查通过：
  - `python3 .github/workflows/hist_acc_restart_fields.py --check` → 525 fields match；
  - `python3 tests/material_transport/check_defineh.py` → 93 unique CI cases PASS。
- 仍缺：`make all` 多配置完整编译、7920 最小集成/性能（N01 子步上限 240×、N05 IO）、N05 NetCDF 重启对照、
  N01 BIF/levee provider 子步。

## 全量编译验证（本轮，MaterialTransport + GRID + USEMPI + GridRiverLakeFlow）

- 环境：macOS arm64，gfortran 16.1（mpif90/OpenMPI 5）+ Homebrew NetCDF-Fortran，临时副本 `/tmp/colm_build`
  （`include/define.h` 打开 `MaterialTransport`），`make -j1 all`。
- 结果：**MAKE_RC=0**，`mksrfdata.x/mkinidata.x/colm.x/postprocess.x` 全部链接成功（“Making all CoLM programs successfully”）。
- 该编译暴露并已修复的真实缺陷（此前单元测试/预处理未发现）：
  1. `mt_riverlake_release_pending_drain` 缺 `USE MOD_Grid_RiverLakeNetwork, only: numucat` 与
     `USE MOD_MT_Defs, only: trc_tiny`；
  2. `MOD_MT_Evapo` 缺 `trc_water_min_for_ratio` 导入（N04 恢复路径）；
  3. coupling 调用 `mt_riverlake_diag_accumulate` 用了旧关键字 `hflux_fc=`（应为 `hflux_fc_in=`）；
  4. `mt_rlc_end_routing` 中 `release_pending_drain` 被插到声明段之前（改为声明后调用）。
- 追加配置编译（同一环境）：
  - **GRID + USEMPI + GridRiverLakeFlow + MaterialTransport**：MAKE_RC=0；
  - **BGC + CROP + LULC_IGBP_PFT + USEMPI + MaterialTransport**（真实甲烷/灌溉宿主配置）：MAKE_RC=0，
    验证 F03/F08 甲烷模块与 N02 CROP 灌溉接口可完整编译；
  - **UNSTRUCTURED + USEMPI + MaterialTransport**：MAKE_RC=0。
- 说明：以上是本机 macOS 全量编译，不是 7920/CI 构建；未运行模型（无强迫数据/完整重启用例）。
