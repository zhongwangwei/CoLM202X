# 00a：VIC 注释闭合修复

## 范围与理由

用户已明确要求直接修复。唯一模型源码改动是在 `main/HYDRO/MOD_Hydro_VIC_Variables.F90` 的既有 `bubble` 行末补齐 `*/`。

该修复与锁定添加源同文件一致，来自添加源历史 `ebd701e5` 中的这一独立修改；不整体 cherry-pick 该提交，不夹带其他功能。相关测试 `tests/test_hydro_vic_cpp_schema.py` 原样取自锁定添加源。

原来的 `/*` 被 C 预处理器识别后跨行吞掉 `zwtvmoist_zwt` 声明，导致消费者无法编译。补齐终止符恢复已有声明，不增加物理计算，不改接口参数或默认值；不是新增正文解释性注释。

## 验证

- 未修复主线作为负对照，预处理回归测试按预期失败。
- 修复版本运行两个预处理/结构成员编译测试：`2 passed`。
- GNU 13.3.0 / MPICH 5.0.1 / netCDF-Fortran 4.6.3，标准 `-cpp` 下 `make -j1 all` 完成，退出 0；包含 surface、初始化、主程序、后处理和静态库。
- 修复后的标准预处理与未修改主线 `-Wp,-C` 预处理，剔除 Fortran 注释后该模块的可执行声明/语句逐行一致。
- 原版标准构建没有可运行结果，不能宣称与其数值逐位一致。整模型运行及与注释保留原版的同环境对照仍待完成。

测试命令：

```bash
PATH=/opt/mpich-gnu/bin:/usr/bin:/bin \
LD_LIBRARY_PATH=/opt/netcdf-gnu/lib:/opt/hdf5-gnu/lib:/opt/mpich-gnu/lib \
/usr/bin/python3 -m pytest -q tests/test_hydro_vic_cpp_schema.py
```

证据存放在工作区 `artifacts/vic-comment-fix/`，包括完整构建日志、测试输出、负对照及预处理等价结果。默认 Conda Python 未安装 pytest，故使用已存在的系统 Python/pytest，未安装或升级共享环境。

## 状态

编译缺陷修复已验证；整模型兼容回归尚未通过。本批不包含冠层物理迁移。安全回退采用对该独立修复提交的 `git revert`，不重置其他工作。
