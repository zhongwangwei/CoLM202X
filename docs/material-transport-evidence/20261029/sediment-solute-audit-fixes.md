# 泥沙/溶质审计（基线 b100745f）修复记录

- 分支 `feature/constituent-transport`；审计基线 `b100745f`；工作 HEAD 见文末（提交后更新）。
- 审计包 `CoLM-sediment-solute-audit-b100745.zip` 已解压；7 项确认（N06/N07 为上一轮 F02 转发修复引入）。
- 本轮修复：**N02、N06**；**N07 部分**；N01/N03/N04/N05 未实施。
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

## N01 / N03 / N04 / N05 —— 未实施

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
