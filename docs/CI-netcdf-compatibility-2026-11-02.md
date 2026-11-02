# 系统 NetCDF-Fortran 与矩阵编译器的兼容性修复

## 问题与授权

基于个人仓库 master 的 `f73fc85b748a8c7eccf7c9cd0f29b7a8b508a8c4`，本报告对应本地依赖兼容性修复，尚未推送。

恢复后的运行 https://github.com/zhongwangwei/CoLM202X/actions/runs/41439327726 的9个 ubuntu-latest 作业失败于模型 USE netcdf，代表 GNU12/MPICH、GNU13/IntelMPI、GNU14/OpenMPI 日志均为系统 netcdf.mod 来自不同 GNU Fortran 版本。该次 latest 实际为 Ubuntu26.04.1 / ubuntu-26.04 / 20261025.179.1；此前旧运行的 latest 为 Ubuntu24.04，不能沿用旧镜像判断。

Ubuntu26.04 的系统 libnetcdff-dev 4.6.2+ds-1 中 netcdf.mod 解压头为 `GFORTRAN module version '16'`。GNU12–14 读取失败已实测复现；GNU15读取成功。GCC 官方说明了版本16格式及新编译器读取旧格式的方向：https://gcc.gnu.org/pipermail/gcc-cvs/2025-January/415463.html 。不能编辑模块头伪装兼容。

用户明确同意：保持官方 workflow 和矩阵不变，在系统 Fortran 依赖不兼容时重建匹配版本，复用系统 NetCDF-C；先验证，再展示摘要等待推送确认。这是一条新增依赖构建路径，不能描述成“仍只使用系统 Fortran 包”。官方文档：https://docs.unidata.ucar.edu/netcdf-c/4.9.3/building_netcdf_fortran.html 。

## 逐文件变更

| 文件 | 改动与边界 |
|---|---|
| include/Makeoptions.github | GNU配置启用依赖兼容性准备；保留指定编译器、MPI后端与用户覆盖 |
| Makefile | 在已确定CaMa是否启用后接入准备逻辑，匹配include路径优先；不回退模型对象或依赖，不更改模型源码 |
| tools/netcdf.mk | 正常构建前准备并重验依赖；使用实际FF，启用CaMa时额外验证FCMP；clean/dry-run不执行下载或构建；只在原LDFLAGS前添加匹配库搜索目录，保留静态链接等用户选项顺序 |
| tools/prepare_netcdf.py | 先实际编译、链接和运行NetCDF-4标量/数组I/O。兼容则直接使用指定系统库；不兼容则校验官方源码SHA256，执行完整configure/build/make check/install及I/O复验；失败不发布可用配置。缓存依据编译器、系统模块/库、NetCDF-C配置、脚本与MPI后端，复用前重新探测；CaMa开关不造成重复构建 |
| .gitignore | 忽略项目私有 .netcdf 缓存；make clean不删该缓存，因此88配置循环可复用；不写系统目录 |
| tests/check_netcdf_prepare.py | 校验损坏源码拒绝、依赖测试失败禁止安装、缓存重验/损坏拒绝、CaMa不兼容拒绝、CaMa开关共享依赖缓存键 |
| tests/check_netcdf_make.py | 验证准备先于模型、重启次数、缓存失败阻断模型、clean/dry-run无下载、FF选择、CaMa参数传递与-static链接参数顺序 |
| tests/netcdf_mpi_probe.F90 | 可独立编译运行的真实两进程NetCDF-4标量/数组I/O与MPI归约检查，不新增到官方CI |

`.github/workflows` 无本次改动，仍与PR base `ac659a3d09cf378b3e90f7b34880ceb78c22a40d` 完全一致。甲烷MPI修复、模型功能及测试宏配置均保留。

## 使用与缓存行为

正常 `make all` 会自动探测；兼容时无需下载。不兼容时下载并校验 NetCDF-Fortran4.6.4（系统C>=4.9.2）；对较旧的受支持C库使用4.5.4。依赖编译不继承模型的默认实数精度或其他模型编译标志，完整上游 `make check` 不得跳过。NetCDF-C/HDF5继续使用系统库。

缓存默认在 `.netcdf/<工具链指纹>/`，含系统探测日志、完整依赖构建日志和上游test-suite.log。显式损坏或失败缓存会报错并指出目录，不伪装通过。用户可提供兼容的 NETCDF_INC/NETCDF_LIB；原自定义Makeoptions不会被整体覆盖，只有GNU配置默认接入自动准备。NETCDF_AUTO=0可用于明确自行管理依赖的用户构建，不用于官方CI规避检查。

## 验证记录

本地：六项依赖准备测试、make集成测试、编译器覆盖测试、9种材料宏测试、525个历史字段检查、material_naming标识检查与18项单/四rank历史I/O检查、run_coupled辅助函数自检、Python语法与diff-check通过。兼容的macOS/Homebrew库实测首次和缓存使用均不重建。

远程重计算全部通过run-7920私有快照，使用Ubuntu26.04.1、系统NetCDF-C4.9.3，GNU12.5/13.4/14.3由官方发行版deb提取到私有目录，不更改服务器系统。系统Fortran包同样提取验证，不上传外部科研数据。

- `codex-889f952341e643ae8a0444f5c2ac1855`：GNU12系统模块不兼容复现；匹配NetCDF-Fortran4.6.4完整构建、make check、安装和I/O验证，终态succeeded。六个test-suite.log合计61个报告测试，零skip/fail/error。
- `codex-03ade774d6e44203bb0c2379b5a957de`：GNU12+MPICH4.3.2，官方第一条配置GRID/PFT/URBANON/Campbell/CaMaON/BGC/CROP的make all通过，包含模型、预处理、后处理和静态库。后续ldd发现Conda wrapper优先选择了Conda C4.9.2，因此该结果单独不算系统C验证。
- `codex-3e99f819b4074600811218fbf79fdb62`：GNU12两rank检查与GNU15系统模块首次/缓存探测通过；GNU15 manifest明确fallback=false。该早期MPI结果也受上述Conda C搜索路径影响。
- `codex-6394f295603944efaf19d15a08dc6c16`：复用GNU12已编译模型对象，给MPICH后端明确匹配Fortran路径及系统C目录优先级，重新make all及两rank I/O，终态succeeded。模型与I/O程序ldd均确认 `/usr/lib/x86_64-linux-gnu/libnetcdf.so.22`，两rank检查通过。此为增量重链接证据，新编译证据来自03ade任务。
- GNU13：`codex-61ad7db980914db99d369159d6d6350f` 终态succeeded，依赖完整make check、两rank探测及官方第一配置的make all通过；后两者的系统C路径在下一项纠正复验。
- GNU14：`codex-664fbdd2182f4846925de6f1f1f58ea6` 终态succeeded，同上。
- GNU13系统C复验：`codex-9a556b075dc84008a78faf3521f3cadc`；GNU14系统C复验：`codex-9ad9244032724a1eb97b51548bd34458`。两者均终态succeeded，复用对应已编译模型对象重新make all及两rank I/O；模型与探测程序ldd均确认系统libnetcdf.so.22路径，两rank检查通过。

GNU12/13/14每套依赖报告均为61个测试、零skip/fail/error。make集成测试在GNU make3.81与4.x上通过。全部远程结果保留在任务目录；已取回的日志/摘要位于 `/tmp/colm-netcdf-checks/`。本次实际完整模型配置为官方第一条，未把它宣称为每版本88条配置全编译。

常用复验命令：

```sh
python3 tests/check_netcdf_prepare.py
python3 tests/check_netcdf_make.py
python3 tests/check_makeoptions.py
python3 tests/create_material_defineh.py --self-test
python3 tests/history_restart/hist_acc_restart_fields.py --check
python3 tests/material_naming/check_history_io.py
python3 tests/material_host_cleanup/run_coupled.py --self-test
git diff --exit-code ac659a3d09cf378b3e90f7b34880ceb78c22a40d -- .github/workflows
```

实际依赖验证调用 `python3 tools/prepare_netcdf.py --fc '<指定编译器>' --check-fc '<指定CaMa编译器>' --include '<系统模块目录>' --lib '<系统库目录>' --cache-dir '<私有输出目录>'`；MPI程序用对应mpif90和匹配Fortran库编译，`mpiexec -n 2` 运行。远程任务的日志保留完整实际编译命令与前缀，避免把示意命令误作执行证据。

## 尚存边界

尚未推送本次修复，恢复后23项hosted矩阵仍不能宣称通过。远程MPICH/私有提取编译器不等同GitHub全部镜像及MPI安装方式；本次尚未实际验证OpenMPI/IntelMPI的新增依赖构建路径、GNU9–11回退构建或4.5.4分支。原始88配置不分片的六小时风险仍存在，下载/构建也会增加首次耗时。无完整科学运行或科学正确性新增结论；此前旧私有NetCDF自检失败的根因未被这份新环境成功记录解释。

## a65d623a hosted失败后的修复

个人运行41447569138的代表GNU12/OpenMPI及GNU14/MPICH日志证实：自动重建启动后，make check失败，尚未编译模型。旧助手只保留runner文件日志，没有输出具体失败项，不能据此断言所有9项都是同一测试失败。

配对复现发现前次验证遗漏的命令差异：CI传入 `mpif90 -fopenmp`，此前完整依赖验证用普通GNU命令。任务 `codex-768a0a362dd44d8492cf8954736cd46a` 给相同GNU12追加-fopenmp后，61项只有nf03_test4/f90tst_io退出139，其余60项通过。该程序有25,000,000字节实数数组；GNU的-fopenmp隐含-frecursive，使数组从静态存储移到栈。任务 `codex-377f9a60a589431394ebb2fadb223503` 证实原失败二进制在默认8192KiB栈下失败、仅诊断性放宽栈后成功；放宽栈不进入产品修复。

修复仅从串行NetCDF依赖的FC/F77命令移除精确的-fopenmp参数，保留MPI wrapper、所选编译器、模型OpenMP选项和原兼容性探针；依赖源码没有OpenMP指令，不关闭模型功能。完整make check仍必须通过。新增失败日志输出会打印构建日志尾部及test-suite.log，使后续hosted失败可直接定位。

- `codex-edf44f693bd24c7cb0f4ba616a14c856`：相同带-fopenmp输入经过修复后构建依赖，61/61、零skip/fail/error，安装及原带-fopenmp的I/O探测通过，终态succeeded。
- `codex-d4a603e1f8b248a7858b05b10de8ddd1`：从官方镜像下载与CI一致的libnetcdf22 4.9.3-1build2 amd64v3包，在私有目录提取；其缓存默认值与普通包相同，ftst_vars/f03tst_vars均通过，未证实缓存默认值假设。
- `codex-d013f22232884b29bf999f70efb78be1`：私有复制修正后的依赖构建树，删除全部测试日志与.trs后，指定上述amd64v3运行库重新执行完整make check；61/61、零skip/fail/error，终态succeeded。
- 两次服务器镜像下载因网络失败；后由本地下载该约0.5MB公共依赖包并随验证源码快照传入私有任务，不含外部科研数据，不改变模型依赖下载地址。
- 七项助手单测及make集成检查通过，包含依赖编译不携带-fopenmp和失败测试详情输出。官方workflow及模型代码未改变。

本次新改动仍须推送后验证完整hosted矩阵，尚不能把本地61项成功称为9个失败作业已全部修好。
