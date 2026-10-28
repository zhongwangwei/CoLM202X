# 3 个 cancelled job 的状态与待补日志（run 40725517808）

终态来源：公开 GitHub web job 页面经 7920l SSH 只读中继解析（`web-terminal.json`），
因本机 HTTPS 中断，Actions REST 日志端点暂不可用；以下 job 的完整日志摘录待网络恢复后由
`poll_curl.py` -> `collect_evidence.py` 自动补入本目录（错误行 + 末尾 120 行 + full.log.gz）。

| job | 结论 | 取消时间(UTC) | 运行时长 | 已完成的步骤 | 被取消的步骤 |
| --- | --- | --- | --- | --- | --- |
| Build-CoLM-gnu (ubuntu-latest, mpich, gcc, 12) | cancelled | 2026-10-28T22:23:56Z | 6h00m12s | Set up job / checkout / setup-mpi / setup-fortran / install netcdf+blas / Test mpi / Check the history-accumulator restart field list / Check effective test case macros（均 success） | Build CoLM202X（93 例逐例 make clean && make all，未在 6h 内完成） |
| Build-CoLM-gnu (ubuntu-latest, mpich, gcc, 13) | cancelled | 2026-10-28T22:23:57Z | 6h00m13s | 同上（8/8 success） | Build CoLM202X |
| Build-CoLM-gnu (ubuntu-latest, mpich, gcc, 14) | cancelled | 2026-10-28T22:23:56Z | 6h00m12s | 同上（8/8 success） | Build CoLM202X |

判定：取消原因是 GitHub Actions 单 job 默认 360 分钟超时（3 个 job 起于 16:23:44Z，
结束于 22:23:56–57Z = 6h00m12–13s），不是编译报错；但这 3 个配置的 93 例完整构建**未完成**，
按验收口径必须记为未通过/未完成，不得并入“完整矩阵通过”。

待补：本目录下将出现 `<job_name>.txt`（错误行摘录 + 末尾日志）与 `<job_name>.full.log.gz`（完整日志），
并在 `log-markers.tsv` 增加对应行。
