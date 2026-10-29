# 甲烷完整窗口所需 15 个边界文件的补传计划（**尚未执行**）

状态：用户 2026-10-29 已批准“可以，补上”。全部测试自 10:03 暂停，**本补传尚未启动**；本文是恢复时可直接照做的
执行单。来源与文件清单来自只读审计（`../methane-closure/GAP_ANALYSIS.md` §1，原始请求台账
`boundary_request_table.md`），当时已确认 BD 源文件存在且尺寸与 1995 同变量一致。

## 目标窗口与缺口

- 目标窗口：**1995-01-01 00:00 → 1996-12-31 21:00 = 35,088 步**（731 天 × 48）。
- 缺 15 个文件：1994 年 7 个变量（`dlwrf,pre,pres,spfh,tmp,ugrd,vgrd`，无 `dswrf`）+ 1997 年全部 8 变量；
  每个文件只读 1 条记录（1994 年末记录 1460 / 1997 年首记录 1），但 reader 按绝对 record 下标寻址，
  因此安全单位是整年文件。**预计总字节 22,706,115,044 B ≈ 21.15 GiB**（7×1994 = 9.87 GiB，8×1997 = 11.28 GiB）。

## 源/目标

| 项 | 路径 |
| --- | --- |
| 源（BD，只读） | `bd:/share/home/dq013/zhwei/colm/data/CoLM_Forcing/crujra_v2.5/` |
| 目标（候选 forcing 根） | `/media/zhwei/data01/compute_runs/data-stage/mt-ch4-global2yr-20261028-2053-1793191930/forcing/` |
| 文件清单（相对路径，15 行） | `transfer-file-list.txt`（本目录） |
| 目标已有 1995/1996 | 保持不变；只新增 1994/1997 文件，不覆盖已有文件 |

## 执行命令（恢复时在 7920l 上执行；建议用 run-7920 小作业或登录节点并保存日志）

```sh
STAGE=/media/zhwei/data01/compute_runs/data-stage/mt-ch4-global2yr-20261028-2053-1793191930
SRC=bd:/share/home/dq013/zhwei/colm/data/CoLM_Forcing/crujra_v2.5/
FILES=$STAGE/transfer-file-list.txt        # 15 行，内容见本目录 transfer-file-list.txt

# 1) 传输前：BD 端 15 文件尺寸（只读）
rsync -rlt --dry-run --files-from="$FILES" "$SRC" "$STAGE/forcing/" > "$STAGE/forcing-补传.dryrun.log" 2>&1

# 2) 传输（沿用原传输方式）
rsync -rlt --safe-links --partial --stats --files-from="$FILES" "$SRC" "$STAGE/forcing/" \
  > "$STAGE/forcing-补传.rsync.log" 2>&1

# 3) 传输后：文件数/字节数/尺寸对照，写入状态 JSON
#    - 期望新增 15 个文件、22,706,115,044 B
#    - 每个文件尺寸与 BD `ls -l` 及 1995 同变量文件对照
#    - 抽查读 NetCDF 第 1 条记录（时间坐标/变量范围）
```

## 验证（必须全过才算输入闭合）

1. 15 个文件均存在且可读；大小与 BD 端一致（记录到 `forcing-补传.target.json` / `…source-after.json`）。
2. 目标原有 1995/1996 文件数量/尺寸不变（与 `forcing.source-before.json` 对照）。
3. 用生产时间索引例程重放 **35,088 步窗口**：LB/UB 请求全部落在已存在文件内，缺失请求数 = **0**。
4. 更新候选包的 `manifest.json`/provenance（新增文件清单、字节数、传输 argv、日志路径、校验结果），
   仍标注科学身份 **UNVERIFIED_CANDIDATE**（传输一致 ≠ 科学正确）。
5. 任何缺失/损坏/尺寸不符：保留失败日志，**不得写成闭合**，不得启动两年计算。

## 依据的原始审计记录（只读，未重新验证 BD）

- `GAP_ANALYSIS.md` §1：15 文件清单、记录号、总字节估算、BD 源存在性（当时只读 `ls -l` 实测）。
- `evidence/boundary_request_table.md`：23,407 条请求台账的逐行汇总（窗口 LB/UB 需求）。
- 候选 staging 的 `forcing.rsync-status.json`：原传输 argv（`rsync -rlt --safe-links --partial --stats`）。
