# 物质输送不支持 CATCHMENT（流域单元）网格：限制与代码清理

状态记录：2026-10-29。用户要求：“物质输送模块不支持流域单元网格。加入限制先。define.h 里面要加好。
如果已经有基于流域网格的物质输送的相关的代码，全部删掉。”本轮只做代码/配置/文档改动与**本机预处理级验证**；
7920 模型测试仍按用户 10:03 的要求暂停，全量构建/CI 验证待恢复测试后执行。

## 1. 限制（三层）

| 层 | 实现 | 行为 |
| --- | --- | --- |
| 受控 `define.h` | `include/define.h`（仓库跟踪的默认头）中的 MaterialTransport 段加入
`#if (defined MaterialTransport && defined CATCHMENT)` → `#error "MaterialTransport does not support the catchment (unit-catchment) grid (CATCHMENT); use GRID or UNSTRUCTURED (vector)."` | 使用该头的构建在预处理期失败 |
| 源码兜底 | `main/MATERIAL/MOD_MT_Defs.F90` 顶部同一条 `#error` | 任何手写/生成的 `define.h` 打开两宏都会 fail-fast |

注：`create_defineh.bash` 不属于本仓库维护范围（用户指示“退回修改”）；本轮**未修改**该脚本，
配置期拒绝仅在 `define.h`/源码层实现。

`GRID`（gridded）与 `UNSTRUCTURED`（vector）仍被允许；`CATCHMENT` 在 **MTOFF** 下不受影响
（只是不能与 MTON 组合）。

## 2. 删除的基于流域网格的 MT 代码

`main/MATERIAL/MOD_MT_Hist.F90` 中所有 CATCHMENT 相关代码已删除（当前 `grep CATCHMENT main/MATERIAL/`
只剩 `MOD_MT_Defs.F90` 的那条限制）：

- 8 处 `#if (defined UNSTRUCTURED || defined CATCHMENT)` → 改为 `#ifdef UNSTRUCTURED`；
- 18 个 `#ifdef CATCHMENT ... #else ... #endif` 分支（HRU/`hydrounit` 维度、`MOD_LandHRU`/`MOD_HRUVector`、
  `hru_patch%substt/subend/subfrc`、`hru_data_address`、`numhru/totalnumhru`、`hydrounit` 输出维度等）
  → 删除 CATCHMENT 分支，保留 UNSTRUCTURED（`MOD_LandElm`/`MOD_ElmVector`/`elm_patch`/`element`）路径；
- 结果：`grep -rn "MOD_LandHRU\|MOD_HRUVector\|hru_patch\|hru_data_address\|numhru\|hydrounit" main/MATERIAL/`
  **无命中**；vector history writers 仅由 `UNSTRUCTURED` 编译。

## 3. 本机验证（预处理/脚本级，不是模型运行）

环境：macOS 本机 `gfortran`（仅预处理 `-E`）、`bash`；用干净临时目录生成三套 `define.h`：

| 测试 | 命令要点 | 结果 |
| --- | --- | --- |
| `CATCHMENT`+MTOFF 仍可用 | 临时 define.h：`#define CATCHMENT`、`#undef MaterialTransport` | 预处理通过（允许） |
| 编译期兜底 | 手工在 CATCHMENT define.h 里打开 `MaterialTransport`（含仓库 `include/define.h` 的同一限制），`gfortran -cpp -E` 预处理 `MOD_MT_Defs.F90` | **rc=1**，报出限制消息 |
| GRID+MT 正常 | `create_defineh.bash GRID ... MTON` + `gfortran -cpp -E` 预处理 `MOD_MT_Defs.F90` / `MOD_MT_Hist.F90` | rc=0；`MOD_MT_Hist` 中 vector writers **0** 处（gridded 路径） |
| UNSTRUCTURED+MT 正常 | `create_defineh.bash UNSTRUCTURED ... MTON` + 预处理 `MOD_MT_Hist.F90` | rc=0；vector writers **4** 处（2d/3d/ratio/delta，UNSTRUCTURED 路径保留） |
| 无流域残留 | `grep -rn "CATCHMENT" main/MATERIAL/`；`grep -rn "MOD_LandHRU\|hru_patch\|hru_data_address\|numhru\|hydrounit" main/MATERIAL/` | 前者只剩 `MOD_MT_Defs.F90` 限制行；后者**无命中** |

## 4. 待恢复测试后的验证（明确未做）

- 7920 上真实 `make all`：`GRID+MT` 与 `UNSTRUCTURED+MT` 应仍编译通过；
  `CATCHMENT+MTON` 必须在配置期/预处理期 fail-fast（不可进入编译/运行）。
- 现有 CI（93 个 GRID 用例）应保持全 PASS；本轮未改 `TestCaseLists`，未加 CI 用例。
- 本轮不改 MT 物理实现，不涉及模型数值；在恢复测试前不声称构建/CI 验证完成。
