# 组分输运宿主钩子清单

> **历史记录。** 本文保留原阶段的设计、决定和验收记录；“待做”“未修”“进行中”以及旧命名不代表当前状态。最新进度以[物质输送 Draft PR 状态](../../material-transport-status.md)为准。旧比较器存在缺文件、缺变量和失败退出码处理缺陷，本文中的“逐位一致”记录须重新核验，不能直接作为本次合并依据。历史命令未作为本次可复现入口重新执行。

源：CoLM202X 工作树；目标：`feature/constituent-transport` 工作树。
由 `docs/superpowers/tools/ct_hook_inventory.py` 生成，移植时在“状态”列打勾：
`ok` 已挂，`drop` 不需要（写理由），`new` 换成新写法（写位置）。

合计 35 个宿主文件，419 个最外层 TRACER 块，块内 5083 行；另有块外 USE 0 处。

## 汇总

| 文件 | 目标分支有此文件 | 块数 | 块内行数 | 块外 USE |
| --- | --- | --- | --- | --- |
| main/CoLMMAIN.F90 | 是 | 31 | 680 | 0 |
| main/HYDRO/MOD_Grid_RiverLakeFlow.F90 | 是 | 54 | 643 | 0 |
| mksrfdata/Aggregation_MethanePH.F90 | **否** | 1 | 562 | 0 |
| main/MOD_SnowLayersCombineDivide.F90 | 是 | 84 | 435 | 0 |
| main/MOD_LeafInterception.F90 | 是 | 17 | 408 | 0 |
| main/MOD_LeafTemperature.F90 | 是 | 11 | 272 | 0 |
| main/LULCC/MOD_Lulcc_Driver.F90 | 是 | 9 | 261 | 0 |
| main/MOD_SoilSnowHydrology.F90 | 是 | 41 | 247 | 0 |
| main/HYDRO/MOD_Hydro_SoilWater.F90 | 是 | 17 | 236 | 0 |
| mksrfdata/Aggregation_LakeSoilC.F90 | **否** | 1 | 197 | 0 |
| main/BGC/MOD_BGC_Soil_BiogeochemNStateUpdate1.F90 | 是 | 7 | 186 | 0 |
| main/BGC/MOD_BGC_CNCStateUpdate1.F90 | 是 | 5 | 116 | 0 |
| main/MOD_Hist.F90 | 是 | 15 | 110 | 0 |
| share/MOD_Namelist.F90 | 是 | 6 | 102 | 0 |
| main/MOD_LeafTemperaturePC.F90 | 是 | 11 | 98 | 0 |
| main/MOD_Thermal.F90 | 是 | 18 | 94 | 0 |
| main/MOD_NewSnow.F90 | 是 | 3 | 88 | 0 |
| mkinidata/MOD_Initialize.F90 | 是 | 6 | 61 | 0 |
| main/BGC/MOD_BGC_Soil_BiogeochemCompetition.F90 | 是 | 2 | 53 | 0 |
| main/CoLM.F90 | 是 | 18 | 51 | 0 |
| main/MOD_Vars_TimeVariables.F90 | 是 | 4 | 37 | 0 |
| main/CoLMDRIVER.F90 | 是 | 12 | 21 | 0 |
| main/BGC/MOD_BGC_CNSummary.F90 | 是 | 2 | 21 | 0 |
| mksrfdata/MOD_RegionClip.F90 | 是 | 2 | 20 | 0 |
| main/MOD_PhaseChange.F90 | 是 | 8 | 18 | 0 |
| main/MOD_Vars_1DAccFluxes.F90 | 是 | 7 | 18 | 0 |
| main/LULCC/MOD_Lulcc_Initialize.F90 | 是 | 4 | 9 | 0 |
| main/MOD_GroundTemperature.F90 | 是 | 4 | 7 | 0 |
| mksrfdata/MKSRFDATA.F90 | 是 | 3 | 7 | 0 |
| main/MOD_Vars_TimeInvariants.F90 | 是 | 5 | 6 | 0 |
| main/HYDRO/MOD_Grid_RiverLakeHist.F90 | 是 | 3 | 6 | 0 |
| main/MOD_GroundFluxes.F90 | 是 | 3 | 5 | 0 |
| include/define.h | 是 | 2 | 4 | 0 |
| main/HYDRO/MOD_Grid_RiverLakeTimeVars.F90 | 是 | 2 | 2 | 0 |
| extends/interception/MOD_Thermal_CanopyPhase_Extended.F90 | **否** | 1 | 2 | 0 |

## main/CoLM.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 105–112 | `TRACER` | 6 | `USE MOD_Tracer_LandPhase, only: land_tracer_init, land_tracer_final` | |
| 113–117 | `TRACER` | 3 | `USE MOD_Tracer_Forcing, only: tracer_forcing_init, read_tracer_forcing, &` | |
| 165–167 | `(defined TRACER) && (defined CaMa_Flood)` | 1 | `integer :: itrc_cama` | |
| 173–176 | `TRACER` | 2 | `logical :: tracer_loaded_restart` | |
| 342–347 | `TRACER` | 4 | `CALL land_tracer_init (numpatch, maxsnl, nl_soil, s_month, lc_year, jdate, &` | |
| 372–375 | `TRACER` | 2 | `CALL tracer_forcing_init (gforc, numpatch)` | |
| 393–398 | `TRACER` | 4 | `IF (p_is_master .and. ntracers > 0) THEN` | |
| 415–425 | `TRACER` | 9 | `IF (DEF_NDEP_FREQUENCY==1)THEN ! Initial annual ndep data readin` | drop：宿主 bug（ndep 年份滞后）的半个修复，不照搬；按计划 P2“宿主 bug 修复”重新做 |
| 433–437 | `TRACER` | 3 | `CALL init_fire_data (s_year)` | drop：同上（火模块人口密度年份滞后） |
| 493–495 | `TRACER` | 1 | `CALL read_tracer_forcing (jdate, dir_forcing)` | |
| 569–573 | `TRACER` | 3 | `CALL CoLMDRIVER (idate,deltim,dolai,doalb,dosst,oroflag,istep)` | |
| 667–670 | `TRACER` | 2 | `CALL tracer_forcing_lulcc_save ()` | |
| 684–687 | `TRACER` | 2 | `CALL tracer_forcing_init (gforc, numpatch)` | |
| 815–818 | `TRACER` | 2 | `CALL tracer_forcing_lulcc_save ()` | |
| 830–833 | `TRACER` | 2 | `CALL tracer_forcing_init (gforc, numpatch)` | |
| 841–843 | `TRACER` | 1 | `CALL tracer_forcing_reset ()` | |
| 860–863 | `TRACER` | 2 | `CALL tracer_forcing_final ()` | |
| 883–886 | `TRACER` | 2 | `CALL tracer_lifecycle_reset ()` | |

## main/CoLMDRIVER.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 3–7 | `TRACER` | 3 | `SUBROUTINE CoLMDRIVER (idate,deltim,dolai,doalb,dosst,oro,istep_in)` | |
| 33–39 | `TRACER` | 5 | `USE MOD_Tracer_LandPhase, only: tracer_resolve_step, tracer_lake_step, &` | |
| 53–55 | `(defined TRACER) && (defined OPENMP)` | 1 | `` | |
| 66–68 | `TRACER` | 1 | `integer,  intent(in), optional :: istep_in` | |
| 73–75 | `TRACER` | 1 | `integer  :: istep_local` | |
| 79–81 | `TRACER` | 1 | `CALL tracer_resolve_step (istep_in, istep_local)` | |
| 105–110 | `(defined URBAN_MODEL) && (defined TRACER)` | 4 | `IF (DEF_URBAN_RUN .and. m.eq.URBAN .and. ntracers > 0) THEN` | |
| 122–124 | `TRACER` | 1 | `waterbody_hist_sample = k == steps_in_one_deltim` | |
| 250–252 | `TRACER` | 1 | `CALL tracer_lake_step (istep_local, i, idate, deltim_phy, k, steps_in_one_deltim)` | |
| 266–268 | `TRACER` | 1 | `CALL tracer_wetland_decomp (i, deltim)` | |
| 270–272 | `TRACER` | 1 | `CALL tracer_soil_step (istep_local, i, idate, deltim)` | |
| 420–422 | `TRACER` | 1 | `CALL tracer_report ()` | |

## main/CoLMMAIN.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 169–271 | `TRACER` | 101 | `USE MOD_LandPFT, only: patch_pft_s, patch_pft_e` | |
| 686–694 | `TRACER` | 7 | `qintr_snow  ,&` | |
| 723–730 | `TRACER` | 6 | `real(r8) :: a, aa, gwat, fevpa_wb, flood_input_wb` | |
| 751–792 | `TRACER` | 40 | `real(r8) :: wextra, t_rain, t_snow` | |
| 854–857 | `TRACER` | 2 | `ldew_rain_old_trc = ldew_rain` | |
| 913–915 | `TRACER` | 1 | `waterstorage_trc_beg = 0._r8` | |
| 918–920 | `TRACER` | 1 | `waterstorage_trc_beg = max(waterstorage(ipatch), 0._r8)` | |
| 927–929 | `TRACER` | 1 | `waterstorage_trc_ground = max(waterstorage_trc_beg - max(qflx_irrig_sprinkler, 0._r8) * de` | |
| 933–950 | `TRACER` | 16 | `ldew_smelt_trc    = 0._r8` | |
| 958–967 | `TRACER` | 8 | `pg_snow,qintr,qintr_rain,qintr_snow,gross_intr_rain,gross_intr_snow,&` | |
| 974–983 | `TRACER` | 8 | `pg_snow,qintr,qintr_rain,qintr_snow,gross_intr_rain,gross_intr_snow,&` | |
| 991–1000 | `TRACER` | 8 | `pg_snow,qintr,qintr_rain,qintr_snow,gross_intr_rain,gross_intr_snow,&` | |
| 1005–1024 | `TRACER` | 18 | `IF (DEF_USE_IRRIGATION) THEN` | |
| 1032–1040 | `TRACER` | 7 | `IF (.not. allocated(wice_snow_bef_trc)) THEN` | |
| 1046–1058 | `TRACER` | 11 | `scv_bef_trc = scv` | |
| 1066–1104 | `TRACER` | 37 | `allocate(soil_thaw_mass_th(lb:nl_soil), soil_frzc_mass_th(lb:nl_soil))` | |
| 1106–1110 | `defined(TRACER) && defined(GridRiverLakeFlow)` | 3 | `flood_evap_temp_trc = t_soisno(lb)` | |
| 1170–1176 | `TRACER` | 5 | `,canopy_smelt_mass_th = canopy_smelt_mass_th, &` | |
| 1179–1247 | `TRACER` | 67 | `IF (allocated(flood_tracer_evap_patch)) THEN` | |
| 1275–1288 | `TRACER` | 12 | `qflx_irrig_drip=qflx_irrig_drip      ,qflx_irrig_flood=qflx_irrig_flood    ,&` | |
| 1320–1444 | `TRACER` | 123 | `qflx_irrig_drip   ,qflx_irrig_flood  ,qflx_irrig_paddy  ,&` | |
| 1459–1485 | `TRACER` | 25 | `IF (ntracers > 0) THEN` | |
| 1487–1505 | `TRACER` | 17 | `IF (ntracers > 0) THEN` | |
| 1511–1580 | `TRACER` | 68 | `IF (ntracers > 0) THEN` | |
| 1633–1670 | `TRACER` | 36 | `flood_input_wb = 0._r8` | |
| 1677–1717 | `TRACER` | 39 | `CALL tracer_apply_reactive_processes(ipatch, snl, nl_soil, deltim)` | |
| 1741–1743 | `TRACER` | 1 | `glacier_overflow_mass_trc = 0._r8` | |
| 1780–1782 | `TRACER` | 1 | `glacier_overflow_mass_trc = glacier_overflow_mass_trc + wextra*deltim` | |
| 1792–1794 | `TRACER` | 1 | `glacier_overflow_mass_trc = glacier_overflow_mass_trc + wextra*deltim` | |
| 1920–1926 | `TRACER` | 5 | `CALL tracer_glacier_patch(ipatch, maxsnl, nl_soil, deltim, &` | |
| 2137–2143 | `TRACER` | 5 | `CALL tracer_waterbody_patch(ipatch, maxsnl, nl_soil, snl, deltim, &` | |

## main/MOD_GroundFluxes.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 28–32 | `TRACER` | 3 | `z0m, z0hg, zol, rib, ustar, qstar, tstar, fm, fh, fq, raw_out)` | |
| 116–118 | `TRACER` | 1 | `real(r8), intent(out), optional :: raw_out` | |
| 234–236 | `TRACER` | 1 | `IF (present(raw_out)) raw_out = raw` | |

## main/MOD_GroundTemperature.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 38–40 | `TRACER` | 1 | `,qphs_thaw_lay,qphs_frzc_lay &` | |
| 166–169 | `TRACER` | 2 | `real(r8), intent(out), optional :: qphs_thaw_lay(lb:nl_soil)` | |
| 429–432 | `TRACER` | 2 | `,qphs_thaw_lay = qphs_thaw_lay, &` | |
| 455–458 | `TRACER` | 2 | `,qphs_thaw_lay = qphs_thaw_lay, &` | |

## main/MOD_Hist.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 21–24 | `TRACER` | 2 | `USE MOD_Tracer_Vars, only: flush_Tracer_Acc` | |
| 37–43 | `TRACER` | 5 | `USE MOD_Grid_RiverLakeHist, only: hist_grid_riverlake_init, hist_grid_riverlake_out, &` | |
| 56–58 | `TRACER` | 1 | `character(len=256) :: file_last_tracer = 'null'` | |
| 69–81 | `TRACER` | 11 | `logical :: flush_reactive` | |
| 191–193 | `TRACER` | 1 | `logical :: forcmask_empty(0)` | |
| 194–203 | `TRACER` | 8 | `real(r8) :: history_window_seconds` | |
| 254–256 | `TRACER` | 1 | `CALL flush_Tracer_Acc ()` | |
| 345–385 | `TRACER` | 39 | `history_window_seconds=0._r8` | |
| 5414–5422 | `TRACER` | 7 | `IF (allocated(forcmask_pch)) THEN` | |
| 5431–5433 | `TRACER` | 1 | `CALL flush_Tracer_Acc ()` | |
| 5442–5444 | `TRACER` | 1 | `file_last_tracer = file_hist_tracer` | |
| 5481–5492 | `TRACER` | 10 | `IF (p_is_worker) THEN` | |
| 5498–5506 | `TRACER` | 7 | `WHERE (acc_vec/=spval .and. acc_num>0)` | |
| 5682–5691 | `TRACER` | 8 | `WHERE (acc_vec(i1,i2,:)/=spval .and. acc_num>0)` | |
| 5702–5711 | `TRACER` | 8 | `WHERE (acc_vec(i1,i2,:)/=spval .and. acc_num>0)` | |

## main/MOD_LeafInterception.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 34–36 | `TRACER` | 1 | `USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite` | |
| 79–81 | `TRACER` | 1 | `PRIVATE :: xsc_rain, xsc_snow` | |
| 180–188 | `TRACER` | 7 | `ldew,ldew_rain,ldew_snow,z0m,hu,pg_rain,pg_snow,qintr,qintr_rain,qintr_snow,&` | |
| 270–313 | `TRACER` | 42 | `real(r8), intent(out), optional :: gross_intr_rain` | |
| 324–329 | `TRACER` | 4 | `satcap_snow = 48._r8*satcap` | |
| 430–439 | `TRACER` | 8 | `FT = max(0._r8, (tleaf - tfrz) / 1.87e5_r8)` | |
| 505–507 | `TRACER` | 1 | `xsc_rain = max(0._r8, ldew)` | |
| 511–513 | `TRACER` | 1 | `xsc_snow = max(0._r8, ldew)` | |
| 531–540 | `TRACER` | 8 | `IF (present(gross_intr_rain))       gross_intr_rain       = max(0._r8, qintr_rain)` | |
| 547–571 | `TRACER` | 23 | `qintr,qintr_rain,qintr_snow,&` | |
| 576–597 | `TRACER` | 20 | `real(r8), intent(out) :: gross_intr_rain, gross_intr_snow` | |
| 601–608 | `TRACER` | 6 | `gross_intr_rain,gross_intr_snow,xsc_rain_out,xsc_snow_out,&` | |
| 615–745 | `TRACER` | 129 | `pg_snow,qintr,qintr_rain,qintr_snow, &` | |
| 754–781 | `TRACER` | 26 | `ldew,ldew_rain,ldew_snow,z0m,hu,pg_rain,pg_snow,qintr,qintr_rain,qintr_snow,&` | |
| 787–854 | `TRACER` | 66 | `integer,  intent(in)    :: ipatch` | |
| 864–923 | `TRACER` | 58 | `IF (present(canopy_phase_heat_p_out)) canopy_phase_heat_p_out(:) = 0._r8` | |
| 933–941 | `TRACER` | 7 | `IF (present(gross_intr_rain)) gross_intr_rain = gross_rain_tmp` | |

## main/MOD_LeafTemperature.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 61–67 | `TRACER` | 5 | `hk         ,hksati     ,rootflux                                        &` | |
| 271–275 | `TRACER` | 3 | `rstfacsun,  &! factor of soil water stress to transpiration on sunlit leaf` | |
| 279–288 | `TRACER` | 8 | `real(r8), intent(inout) :: &` | |
| 454–457 | `TRACER` | 2 | `IF (present(canopy_smelt_mass_out)) canopy_smelt_mass_out = 0._r8` | |
| 1055–1063 | `TRACER` | 7 | `gssun = 0._r8` | |
| 1114–1123 | `TRACER` | 8 | `IF (etr < 0._r8) THEN` | |
| 1205–1208 | `TRACER` | 2 | `IF (DEF_Interception_scheme .eq. 1) THEN` | |
| 1241–1476 | `TRACER` | 234 | `ELSEIF (DEF_Interception_scheme .eq. 2) THEN !CLM4.5` | |
| 1521–1523 | `TRACER` | 1 | `IF (present(canopy_smelt_mass_out)) canopy_smelt_mass_out = qmelt*deltim` | |
| 1533–1535 | `TRACER` | 1 | `IF (present(canopy_frzc_mass_out)) canopy_frzc_mass_out = qfrz*deltim` | |
| 1547–1549 | `TRACER` | 1 | `IF (present(raw_trc_out)) raw_trc_out = max(raw, 0._r8)` | |

## main/MOD_LeafTemperaturePC.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 72–74 | `TRACER` | 1 | `,canopy_smelt_mass_p_out, canopy_frzc_mass_p_out, raw_trc_out &` | |
| 268–272 | `TRACER` | 3 | `real(r8), dimension(ps:pe), intent(out), optional :: canopy_smelt_mass_p_out` | |
| 536–539 | `TRACER` | 2 | `IF (present(canopy_smelt_mass_p_out)) canopy_smelt_mass_p_out(:) = 0._r8` | |
| 566–568 | `TRACER` | 1 | `IF (present(raw_trc_out)) raw_trc_out = 0._r8` | |
| 1805–1815 | `TRACER` | 9 | `gssun = 0._r8` | |
| 1893–1902 | `TRACER` | 8 | `IF (etr(i) < 0._r8) THEN` | |
| 1924–1928 | `TRACER` | 3 | `IF (DEF_Interception_scheme .eq. 1 .or. DEF_Interception_scheme .eq. 8) THEN !colm2014` | |
| 1961–2030 | `TRACER` | 68 | `IF (.not. DEF_VEG_SNOW) THEN` | |
| 2050–2052 | `TRACER` | 1 | `IF (present(canopy_smelt_mass_p_out)) canopy_smelt_mass_p_out(i) = qmelt(i)*deltim` | |
| 2062–2064 | `TRACER` | 1 | `IF (present(canopy_frzc_mass_p_out)) canopy_frzc_mass_p_out(i) = qfrz(i)*deltim` | |
| 2163–2165 | `TRACER` | 1 | `IF (present(raw_trc_out)) raw_trc_out = max(raw, 0._r8)` | |

## main/MOD_NewSnow.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 12–14 | `TRACER` | 1 | `PUBLIC :: relocate_soil_frost_ice` | |
| 80–84 | `TRACER` | 3 | `IF(patchtype==2 .and. t_grnd>tfrz .and. snl==0)THEN` | |
| 133–218 | `TRACER` | 84 | `SUBROUTINE relocate_soil_frost_ice(maxsnl, porsl1, snl, zi, z, dz, t, wliq, wice, &` | |

## main/MOD_PhaseChange.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 34–36 | `TRACER` | 1 | `,qphs_thaw_lay, qphs_frzc_lay &` | |
| 108–111 | `TRACER` | 2 | `real(r8), intent(out), optional :: qphs_thaw_lay(lb:nl_soil)` | |
| 128–131 | `TRACER` | 2 | `IF (present(qphs_thaw_lay)) qphs_thaw_lay(:) = 0._r8` | |
| 321–326 | `TRACER` | 4 | `IF (present(qphs_thaw_lay)) &` | |
| 355–357 | `TRACER` | 1 | `,qphs_thaw_lay, qphs_frzc_lay &` | |
| 430–433 | `TRACER` | 2 | `real(r8), intent(out), optional :: qphs_thaw_lay(lb:nl_soil)` | |
| 451–454 | `TRACER` | 2 | `IF (present(qphs_thaw_lay)) qphs_thaw_lay(:) = 0._r8` | |
| 650–655 | `TRACER` | 4 | `IF (present(qphs_thaw_lay)) &` | |

## main/MOD_SnowLayersCombineDivide.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 254–256 | `TRACER` | 1 | `,trc_wliq, trc_wice, trc_solid, trc_scv &` | |
| 287–292 | `TRACER` | 4 | `real(r8), intent(inout), optional :: trc_wliq(:, lb:)` | |
| 308–312 | `TRACER` | 3 | `real(r8), allocatable :: zwtrc_ice(:), zwtrc_liq(:), zwtrc_solid(:)` | |
| 317–323 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 332–341 | `TRACER` | 8 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 350–358 | `TRACER` | 7 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 360–366 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 370–378 | `TRACER` | 7 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 389–392 | `TRACER` | 2 | `IF (present(trc_wliq) .and. present(trc_wice)) &` | |
| 404–413 | `TRACER` | 8 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 428–441 | `TRACER` | 12 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 479–488 | `TRACER` | 8 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 498–506 | `TRACER` | 7 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 511–519 | `TRACER` | 7 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 544–547 | `TRACER` | 2 | `IF (present(trc_wliq) .and. present(trc_wice)) &` | |
| 554–556 | `TRACER` | 1 | `,trc_wliq, trc_wice, trc_solid &` | |
| 579–583 | `TRACER` | 3 | `real(r8), intent(inout), optional :: trc_wliq(:, lb:)` | |
| 599–603 | `TRACER` | 3 | `real(r8), allocatable :: strc_wice(:,:), strc_wliq(:,:), strc_solid(:,:)` | |
| 607–615 | `TRACER` | 7 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 623–631 | `TRACER` | 7 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 641–647 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 653–659 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 670–676 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 681–687 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 692–698 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 709–715 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 721–727 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 738–744 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 749–755 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 760–766 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 777–783 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 789–795 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 806–812 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 817–823 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 828–834 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 845–851 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 857–863 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 874–880 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 885–891 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 896–902 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 916–924 | `TRACER` | 7 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 933–938 | `TRACER` | 4 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1015–1017 | `TRACER` | 1 | `,trc_wliq, trc_wice, trc_solid, trc_scv &` | |
| 1061–1066 | `TRACER` | 4 | `real(r8), intent(inout), optional :: trc_wliq(:, lb:)` | |
| 1082–1085 | `TRACER` | 2 | `real(r8), allocatable :: zwtrc_ice(:), zwtrc_liq(:), zwtrc_solid(:)` | |
| 1090–1096 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1105–1114 | `TRACER` | 8 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1137–1145 | `TRACER` | 7 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1158–1164 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1168–1176 | `TRACER` | 7 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1199–1202 | `TRACER` | 2 | `IF (present(trc_wliq) .and. present(trc_wice)) &` | |
| 1214–1223 | `TRACER` | 8 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1249–1262 | `TRACER` | 12 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1300–1309 | `TRACER` | 8 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1331–1339 | `TRACER` | 7 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1355–1363 | `TRACER` | 7 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1388–1391 | `TRACER` | 2 | `IF (present(trc_wliq) .and. present(trc_wice)) &` | |
| 1404–1406 | `TRACER` | 1 | `,trc_wliq, trc_wice, trc_solid &` | |
| 1444–1448 | `TRACER` | 3 | `real(r8), intent(inout), optional :: trc_wliq(:, lb:)` | |
| 1469–1473 | `TRACER` | 3 | `real(r8), allocatable :: strc_wice(:,:), strc_wliq(:,:), strc_solid(:,:)` | |
| 1477–1485 | `TRACER` | 7 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1505–1513 | `TRACER` | 7 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1527–1533 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1541–1547 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1564–1570 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1578–1584 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1594–1600 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1614–1620 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1628–1634 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1650–1656 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1664–1670 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1680–1686 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1700–1706 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1714–1720 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1737–1743 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1751–1757 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1767–1773 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1787–1793 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1801–1807 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1824–1830 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1838–1844 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1854–1860 | `TRACER` | 5 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1887–1895 | `TRACER` | 7 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |
| 1905–1910 | `TRACER` | 4 | `IF (present(trc_wliq) .and. present(trc_wice)) THEN` | |

## main/MOD_SoilSnowHydrology.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 7–9 | `TRACER` | 1 | `USE MOD_SPMD_Task, only: CoLM_stop` | |
| 67–73 | `TRACER` | 5 | `qflx_irrig_drip  ,qflx_irrig_flood ,qflx_irrig_paddy, &` | |
| 199–201 | `TRACER` | 1 | `logical, intent(in), optional :: defer_surface_ice_overflow` | |
| 225–227 | `TRACER` | 1 | `real(r8) :: dew_input, dew_capacity, dew_retained, dew_excess, frost_excess, ice_before_fr` | |
| 307–313 | `TRACER` | 5 | `IF (patchtype <= 1) THEN` | |
| 319–322 | `TRACER` | 2 | `gwat = gwat + wdsrf/deltim` | |
| 509–512 | `TRACER` | 2 | `dew_input = 0._r8` | |
| 516–551 | `TRACER` | 34 | `dew_input = max(qsdew*deltim, 0._r8)` | |
| 555–576 | `TRACER` | 20 | `IF ((.not.DEF_SPLIT_SOILSNOW) .or. (patchtype==1 .and. DEF_URBAN_RUN)) THEN` | |
| 655–664 | `TRACER` | 8 | `qflx_irrig_drip  ,qflx_irrig_flood ,qflx_irrig_paddy            ,&` | |
| 682–686 | `TRACER` | 3 | `USE MOD_Vars_TimeInvariants, only: wetwatmax, smpmin` | |
| 792–805 | `TRACER` | 12 | `real(r8), intent(out) :: etroot(1:nl_soil)` | |
| 830–832 | `TRACER` | 1 | `logical, intent(in), optional :: defer_surface_ice_overflow` | |
| 850–856 | `TRACER` | 5 | `real(r8), parameter :: richards_water_tolerance = 1.e-3_r8` | |
| 874–878 | `TRACER` | 3 | `real(r8) :: rsubst_fld` | |
| 884–886 | `TRACER` | 1 | `real(r8) :: fevpg_runoff` | |
| 906–924 | `TRACER` | 17 | `IF (present(permeable_soil_out)) permeable_soil_out = .false.` | |
| 957–959 | `TRACER` | 1 | `,snow_qout_layer(lb:0)&` | |
| 968–970 | `TRACER` | 1 | `,qout_snow_layer=snow_qout_layer(lb:0)&` | |
| 984–986 | `TRACER` | 1 | `,snow_qout_layer(lb:0)&` | |
| 995–997 | `TRACER` | 1 | `,qout_snow_layer=snow_qout_layer(lb:0)&` | |
| 1045–1047 | `TRACER` | 1 | `IF (present(permeable_soil_out)) permeable_soil_out = is_permeable` | |
| 1139–1153 | `TRACER` | 13 | `IF (((patchtype == 2 .AND. DEF_USE_Dynamic_Wetland) .OR. &` | |
| 1268–1305 | `TRACER` | 36 | `qgtop_deficit = -qgtop * deltim` | |
| 1312–1314 | `TRACER` | 1 | `qgtop_out = qgtop` | |
| 1322–1328 | `TRACER` | 5 | `etroot(1:nl_soil),       etroot_actual(1:nl_soil), etroot_aquifer, etroot_surface,       &` | |
| 1350–1356 | `TRACER` | 5 | `rsur_before_late = rsur` | |
| 1360–1364 | `TRACER` | 3 | `dew_input = max(qsdew*deltim, 0._r8)` | |
| 1368–1399 | `TRACER` | 30 | `dew_input = max(qsdew_soil*deltim, 0._r8)` | |
| 1405–1407 | `TRACER` | 1 | `wblc_ice_sink(j) = wblc` | |
| 1412–1414 | `TRACER` | 1 | `wblc_ice_sink(j) = wice_soisno(j)` | |
| 1479–1482 | `TRACER` | 2 | `IF (frost_excess+dew_input-dew_retained > 0._r8) &` | |
| 1555–1570 | `TRACER` | 14 | `DO j = 1, nl_soil` | |
| 1607–1609 | `TRACER` | 1 | `,qout_snow_layer &` | |
| 1651–1654 | `TRACER` | 2 | `real(r8), intent(out), optional :: &` | |
| 1702–1704 | `TRACER` | 1 | `IF (present(qout_snow_layer)) qout_snow_layer(:) = 0._r8` | |
| 1721–1723 | `TRACER` | 1 | `IF (present(qout_snow_layer)) qout_snow_layer(j) = qout` | |
| 1743–1745 | `TRACER` | 1 | `,qout_snow_layer &` | |
| 1792–1795 | `TRACER` | 2 | `real(r8), intent(out), optional :: &` | |
| 1930–1932 | `TRACER` | 1 | `IF (present(qout_snow_layer)) qout_snow_layer(:) = 0._r8` | |
| 1974–1976 | `TRACER` | 1 | `IF (present(qout_snow_layer)) qout_snow_layer(j) = qout` | |

## main/MOD_Thermal.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 81–84 | `TRACER` | 2 | `,canopy_smelt_mass_th, canopy_frzc_mass_th                  ,&` | |
| 145–147 | `TRACER` | 1 | `USE MOD_Namelist, only: DEF_VEG_SNOW` | |
| 286–298 | `TRACER` | 11 | `real(r8), intent(out), optional :: canopy_smelt_mass_th` | |
| 404–411 | `TRACER` | 6 | `real(r8), intent(inout) :: flddepth` | |
| 529–532 | `TRACER` | 2 | `real(r8) :: canopy_smelt_mass_local, canopy_frzc_mass_local, raw_trc_local, raw_trc_pc` | |
| 549–565 | `TRACER` | 15 | `fevpg_fld = 0._r8` | |
| 573–583 | `TRACER` | 9 | `canopy_smelt_mass_local = 0._r8` | |
| 725–730 | `TRACER` | 4 | `z0m_g,z0h_g,zol_g,rib_g,ustar_g,qstar_g,tstar_g,fm_g,fh_g,fq_g, &` | |
| 801–805 | `TRACER` | 3 | `,canopy_smelt_mass_out=canopy_smelt_mass_local, &` | |
| 873–889 | `TRACER` | 15 | `allocate ( canopy_smelt_mass_p_local(ps:pe) )` | |
| 1030–1039 | `TRACER` | 8 | `smp             ,hk(1:)          ,hksati(1:)      ,rootflux_p(1:,i)                &` | |
| 1050–1054 | `TRACER` | 3 | `qstar_p(i),tstar_p(i),fm_p(i),fh_p(i),fq_p(i),raw_out=raw_trc_p(i))` | |
| 1175–1179 | `TRACER` | 3 | `,canopy_smelt_mass_p_out=canopy_smelt_mass_p_local(ps:pe),&` | |
| 1206–1208 | `TRACER` | 1 | `raw_trc_p    (ps:pe) = raw_trc_pc` | |
| 1264–1268 | `TRACER` | 3 | `canopy_smelt_mass_local = sum( canopy_smelt_mass_p_local(ps:pe)*pftfrac(ps:pe) )` | |
| 1327–1331 | `TRACER` | 3 | `deallocate ( canopy_smelt_mass_p_local )` | |
| 1336–1340 | `TRACER` | 3 | `IF (present(canopy_smelt_mass_th)) canopy_smelt_mass_th = canopy_smelt_mass_local` | |
| 1397–1400 | `TRACER` | 2 | `,qphs_thaw_lay=qphs_thaw_lay_th, &` | |

## main/MOD_Vars_1DAccFluxes.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 14–16 | `TRACER` | 1 | `USE MOD_Tracer_LandPhase, only: tracer_flush_acc_fluxes, tracer_accumulate_fluxes` | |
| 1795–1799 | `TRACER` | 3 | `SUBROUTINE FLUSH_acc_fluxes (flush_reactive)` | |
| 1806–1812 | `TRACER` | 5 | `logical, intent(in), optional :: flush_reactive` | |
| 2400–2402 | `TRACER` | 1 | `IF (flush_reactive_active) CALL tracer_flush_acc_fluxes ()` | |
| 2427–2432 | `TRACER` | 4 | `USE MOD_Namelist, only: DEF_USE_CBL_HEIGHT, DEF_USE_OZONESTRESS, DEF_USE_PLANTHYDRAULICS, ` | |
| 2573–2577 | `TRACER` | 3 | `IF (.not. DEF_USE_VariablySaturatedFlow) CALL acc1d (qcharge, a_qcharge)` | |
| 3430–3432 | `TRACER` | 1 | `CALL tracer_accumulate_fluxes ()` | |

## main/MOD_Vars_TimeInvariants.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 234–236 | `(defined TRACER) && (defined BGC)` | 1 | `real(r8), allocatable :: lake_soilc_srf(:,:)` | |
| 380–383 | `(defined TRACER) && (defined BGC)` | 2 | `allocate (lake_soilc_srf(nl_soil,numpatch))` | |
| 544–546 | `(defined TRACER) && (defined BGC)` | 1 | `CALL ncio_read_vector (file_restart, 'lake_soilc_srf', nl_soil, landpatch, lake_soilc_srf,` | |
| 767–769 | `(defined TRACER) && (defined BGC)` | 1 | `CALL ncio_write_vector (file_restart, 'lake_soilc_srf', 'soil', nl_soil, 'patch', landpatc` | |
| 947–949 | `(defined TRACER) && (defined BGC)` | 1 | `deallocate (lake_soilc_srf )` | |

## main/MOD_Vars_TimeVariables.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 1112–1118 | `TRACER` | 5 | `USE MOD_Tracer_Defs, only: ntracers` | |
| 1131–1133 | `TRACER` | 1 | `real(r8) :: empty_patch(0), empty_soilsnow(maxsnl+1:nl_soil,0)` | |
| 1262–1293 | `TRACER` | 30 | `IF (allocated(ldew_rain)) THEN` | |
| 1340–1342 | `TRACER` | 1 | `CALL write_tracer_restart(file_restart)` | |

## main/LULCC/MOD_Lulcc_Driver.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 7–16 | `TRACER` | 8 | `USE MOD_Mesh, only: mesh` | |
| 100–103 | `TRACER` | 2 | `USE MOD_Vars_TimeInvariants, only: patchclass` | |
| 109–117 | `TRACER` | 7 | `USE MOD_Tracer_Defs, only: ntracers` | |
| 130–133 | `TRACER` | 2 | `real(r8), allocatable :: old_patch_area(:), new_patch_area(:), inventory_trace(:,:)` | |
| 149–164 | `TRACER` | 14 | `CALL save_land_tracer_lulcc_state ()` | |
| 176–185 | `TRACER` | 8 | `IF (inventory_active) THEN` | |
| 211–236 | `TRACER` | 24 | `IF (inventory_active) THEN` | |
| 245–271 | `TRACER` | 25 | `IF (p_is_worker .and. allocated(patchclass)) THEN` | |
| 283–455 | `TRACER` | 171 | `SUBROUTINE lulcc_patch_areas(patches, areas)` | |

## main/LULCC/MOD_Lulcc_Initialize.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 53–56 | `TRACER` | 2 | `USE MOD_Tracer_Defs, only: ntracers` | |
| 70–72 | `TRACER` | 1 | `logical :: old_has_patches` | |
| 79–82 | `TRACER` | 2 | `old_has_patches = .false.` | |
| 112–117 | `TRACER` | 4 | `IF (p_is_worker .and. ntracers > 0 .and. DEF_LULCC_SCHEME == 2) THEN` | |

## main/BGC/MOD_BGC_CNCStateUpdate1.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 30–34 | `TRACER` | 3 | `decomp_cpools_vr, I_met_c_vr_acc, I_cel_c_vr_acc, I_lig_c_vr_acc, &` | |
| 147–151 | `TRACER` | 3 | `PUBLIC :: CStateUpdate1, CDecompStateUpdate` | |
| 166–170 | `TRACER` | 3 | `integer j` | |
| 191–237 | `TRACER` | 45 | `CALL CDecompStateUpdate(i, deltim, nl_soil, ndecomp_transitions)` | |
| 508–571 | `TRACER` | 62 | `SUBROUTINE CDecompStateUpdate(i, deltim, nl_soil, ndecomp_transitions, apply_direct)` | |

## main/BGC/MOD_BGC_CNSummary.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 190–192 | `TRACER` | 1 | `PUBLIC CNDriverSummarizeNonvegetatedSoilStates` | |
| 232–253 | `TRACER` | 20 | `SUBROUTINE CNDriverSummarizeNonvegetatedSoilStates(i,nl_soil,dz_soi,ndecomp_pools)` | |

## main/BGC/MOD_BGC_Soil_BiogeochemCompetition.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 34–38 | `TRACER` | 3 | `PUBLIC SoilBiogeochemCompetition, SoilBiogeochemCompetitionNoPlant` | |
| 472–523 | `TRACER` | 50 | `SUBROUTINE SoilBiogeochemCompetitionNoPlant(i, deltim, nl_soil, dz_soi)` | |

## main/BGC/MOD_BGC_Soil_BiogeochemNStateUpdate1.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 31–35 | `TRACER` | 3 | `decomp_npools_vr           , sminn_vr                   , smin_nh4_vr                , smi` | |
| 55–59 | `TRACER` | 3 | `PUBLIC SoilBiogeochemNStateUpdate1, SoilBiogeochemNDecompStateUpdate` | |
| 71–75 | `TRACER` | 3 | `integer j` | |
| 122–177 | `TRACER` | 54 | `CALL SoilBiogeochemNDecompStateUpdate(i, deltim, nl_soil, ndecomp_transitions)` | |
| 184–207 | `TRACER` | 22 | `DO k = 1, ndecomp_transitions` | |
| 226–236 | `TRACER` | 9 | `smin_nh4_vr(j,i) = smin_nh4_vr(j,i) + gross_nmin_vr(j,i)*deltim` | |
| 265–358 | `TRACER` | 92 | `SUBROUTINE SoilBiogeochemNDecompStateUpdate(i, deltim, nl_soil, &` | |

## main/HYDRO/MOD_Grid_RiverLakeFlow.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 25–37 | `TRACER` | 11 | `USE MOD_Tracer_Lifecycle, only: tracer_lifecycle_route_has_active, tracer_lifecycle_route_` | |
| 73–79 | `TRACER` | 5 | `real(r8), allocatable :: flood_tracer_credit_patch(:,:)` | |
| 90–92 | `TRACER` | 1 | `USE MOD_Tracer_Defs,         only: ntracers, tracer_uses_land_water_transport` | |
| 101–106 | `TRACER` | 4 | `logical :: trc_restart_found, has_flood_tracer` | |
| 133–138 | `TRACER` | 4 | `CALL tracer_lifecycle_route_init()` | |
| 172–212 | `TRACER` | 39 | `trc_restart_found = .false.` | |
| 219–232 | `TRACER` | 12 | `has_flood_tracer = .false.` | |
| 257–259 | `TRACER` | 1 | `IF (allocated(flood_tracer_credit_patch)) flood_tracer_credit_patch = 0._r8` | |
| 334–336 | `TRACER` | 1 | `USE MOD_Vars_1DForcing, only: forc_prc, forc_prl` | |
| 337–341 | `TRACER` | 3 | `USE MOD_Tracer_Defs,    only: ntracers, tracer_uses_land_water_transport, &` | |
| 358–362 | `TRACER` | 3 | `real(r8), allocatable :: prcp_gd(:)` | |
| 363–374 | `TRACER` | 10 | `integer  :: itrc, itrc_dep` | |
| 399–402 | `TRACER` | 2 | `integer,  allocatable :: ucat2resv_safe(:)` | |
| 416–419 | `TRACER` | 2 | `real(r8), allocatable :: trc_mass_bef(:), trc_mass_aft(:)` | |
| 445–465 | `TRACER` | 19 | `IF (ntracers > 0) THEN` | |
| 469–471 | `TRACER` | 1 | `IF (ntracers > 0) CALL tracer_input_from_runoff(rnof_uc*1.e-3*deltime, numucat, trc_rnof_u` | |
| 476–479 | `TRACER` | 2 | `IF (allocated(trc_rnof_gd)) deallocate(trc_rnof_gd)` | |
| 481–528 | `TRACER` | 46 | `IF (tracer_lifecycle_route_has_active()) THEN` | |
| 565–573 | `(defined CoLMDEBUG) && (defined TRACER)` | 7 | `allocate (trc_mass_bef(ntracers), trc_mass_aft(ntracers), trc_mass_inp(ntracers), &` | |
| 593–606 | `TRACER` | 12 | `allocate (total_floodarea(numucat), total_flooddepth(numucat))` | |
| 627–654 | `TRACER` | 26 | `IF (allocated(ucat2resv)) THEN` | |
| 703–708 | `TRACER` | 4 | `IF (ntracers > 0) CALL tracer_lifecycle_route_sediment_levee_repartition(i, &` | |
| 1215–1227 | `TRACER` | 11 | `IF (ntracers > 0 .and. DEF_USE_BIFURCATION) &` | |
| 1243–1250 | `TRACER` | 6 | `particle_water_storage_start(i) = max(volwater, 0._r8)` | |
| 1272–1287 | `TRACER` | 14 | `IF (volwater > 1.e-6_r8) THEN` | |
| 1296–1300 | `TRACER` | 3 | `IF (ntracers > 0) CALL tracer_lifecycle_route_sediment_levee_repartition(i, &` | |
| 1330–1333 | `TRACER` | 2 | `IF ((.not. is_built_resv(i)) .and. (wdsrf_ucat(i) >= RIVERMIN)) &` | |
| 1337–1343 | `TRACER` | 5 | `IF (numucat > 0) THEN` | |
| 1351–1356 | `TRACER` | 4 | `DO itrc = 1, ntracers` | |
| 1376–1383 | `TRACER` | 6 | `total_floodarea(i) = floodarea` | |
| 1438–1474 | `TRACER` | 35 | `IF (tracer_lifecycle_route_has_active()) THEN` | |
| 1499–1511 | `TRACER` | 11 | `IF (numucat > 0) THEN` | |
| 1529–1537 | `TRACER` | 7 | `IF (ntracers > 0) THEN` | |
| 1547–1559 | `TRACER` | 11 | `DO itrc = 1, ntracers` | |
| 1561–1563 | `TRACER` | 1 | `deallocate (trc_mass_bef, trc_mass_aft, trc_mass_inp, trc_mass_dis, trc_mass_reactive)` | |
| 1577–1605 | `TRACER` | 27 | `IF (p_is_worker) THEN` | |
| 1615–1618 | `TRACER` | 2 | `IF (allocated(levee_floodarea)) CALL publish_levee_fldfrc_to_patches (levee_floodarea)` | |
| 1641–1648 | `TRACER` | 6 | `IF (allocated(ucat2resv_safe)) deallocate(ucat2resv_safe)` | |
| 1655–1734 | `TRACER` | 78 | `SUBROUTINE publish_levee_fldfrc_to_patches (levee_floodarea_in)` | |
| 1815–1817 | `TRACER` | 1 | `USE MOD_Tracer_Defs, only: ntracers, tracer_uses_land_water_transport` | |
| 1824–1827 | `TRACER` | 2 | `integer :: itrc` | |
| 1835–1843 | `TRACER` | 7 | `IF (allocated(flood_visible_tracer_uc)) THEN` | |
| 1860–1862 | `TRACER` | 1 | `IF (allocated(trc_solid)) CALL equilibrate_river_tracer_cell(i, visible, protected)` | |
| 1873–1879 | `TRACER` | 5 | `IF (allocated(flood_visible_tracer_uc)) THEN` | |
| 1916–1956 | `TRACER` | 39 | `IF (allocated(flood_visible_tracer_uc)) THEN` | |
| 1963–1966 | `TRACER` | 2 | `USE MOD_Tracer_Defs, only: ntracers, tracer_uses_land_water_transport, &` | |
| 1971–1978 | `TRACER` | 6 | `real(r8), allocatable :: gain_patch(:), coefficient_uc(:), gain_uc(:)` | |
| 1988–1994 | `TRACER` | 5 | `IF (allocated(flood_visible_tracer_uc)) THEN` | |
| 2042–2160 | `TRACER` | 117 | `IF (allocated(flood_visible_tracer_uc)) THEN` | |
| 2170–2176 | `TRACER` | 5 | `IF (allocated(trc_solid)) THEN` | |
| 2181–2186 | `TRACER` | 4 | `IF (ntracers > 0) CALL tracer_lifecycle_route_sediment_levee_repartition(j, &` | |
| 2207–2209 | `TRACER` | 1 | `CALL tracer_lifecycle_route_final()` | |
| 2214–2216 | `TRACER` | 1 | `CALL river_lake_tracer_final()` | |
| 2223–2227 | `TRACER` | 3 | `IF (allocated(flood_tracer_credit_patch)) deallocate(flood_tracer_credit_patch, &` | |

## main/HYDRO/MOD_Grid_RiverLakeHist.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 15–18 | `TRACER` | 2 | `USE MOD_Tracer_Lifecycle, only: tracer_lifecycle_route_write_history, tracer_lifecycle_rou` | |
| 577–580 | `TRACER` | 2 | `CALL tracer_lifecycle_route_write_history (file_hist_ucat, itime_in_file_ucat)` | |
| 680–683 | `TRACER` | 2 | `CALL tracer_flush_acc()` | |

## main/HYDRO/MOD_Grid_RiverLakeTimeVars.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 14–16 | `TRACER` | 1 | `USE MOD_Tracer_Lifecycle, only: tracer_lifecycle_route_write_restart` | |
| 683–685 | `TRACER` | 1 | `CALL tracer_lifecycle_route_write_restart(file_restart)` | |

## main/HYDRO/MOD_Hydro_SoilWater.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 21–27 | `TRACER` | 5 | `USE MOD_SPMD_Task, only: CoLM_stop` | |
| 66–68 | `TRACER` | 1 | `PRIVATE :: project_richards_liquid_water` | |
| 177–180 | `TRACER` | 2 | `,etroot_out, etroot_actual_out, etroot_aquifer_out, etroot_surface_out, &` | |
| 230–236 | `TRACER` | 5 | `real(r8), intent(out) :: etroot_out(1:nlev)` | |
| 245–247 | `TRACER` | 1 | `real(r8) :: attempted, ss_vliq_pre, residual_mm` | |
| 261–265 | `TRACER` | 3 | `real(r8) :: exchange_dp_before, pond_exchange` | |
| 320–328 | `TRACER` | 7 | `etroot_out(1:nlev) = etroot(1:nlev)` | |
| 335–363 | `TRACER` | 27 | `attempted   = etroot(ilev)*dt + deficit` | |
| 375–442 | `TRACER` | 66 | `wexchange = rsubst * dt + deficit` | |
| 446–465 | `TRACER` | 18 | `IF (wexchange > 0._r8) pond_exchange = max(exchange_dp_before-ss_dp, 0._r8)` | |
| 562–566 | `TRACER` | 3 | `qinfl = qgtop - (ss_dp - dp_m1 + pond_exchange)/dt` | |
| 616–629 | `TRACER` | 12 | `IF (DEF_USE_PLANTHYDRAULICS .and. &` | |
| 881–886 | `TRACER` | 4 | `real(r8) :: wsum_m1, wsum, werr, mass_budget` | |
| 893–897 | `TRACER` | 3 | `mass_budget = 256._r8 * epsilon(1._r8) * &` | |
| 983–996 | `TRACER` | 12 | `projection_ok = .true.` | |
| 1000–1004 | `TRACER` | 3 | `.or. wet2dry .or. (.not. projection_ok)) THEN` | |
| 1260–1325 | `TRACER` | 64 | `SUBROUTINE project_richards_liquid_water ( &` | |

## share/MOD_Namelist.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 476–479 | `(defined TRACER) && (defined BGC)` | 2 | `character(len=256) :: DEF_file_GIEMS = 'null'` | |
| 1409–1442 | `TRACER` | 32 | `DEF_TRACER_USE_FRACTIONATION, &` | |
| 1533–1536 | `TRACER` | 2 | `IF (len_trim(DEF_TRACER_PARAM_FILES) == len(DEF_TRACER_PARAM_FILES)) &` | |
| 1614–1646 | `TRACER` | 31 | `SELECT CASE (trim(adjustl(DEF_TRACER_KINETIC_SCHEME)))` | |
| 1681–1684 | `TRACER` | 2 | `IF (.not. DEF_USE_VariablySaturatedFlow) &` | |
| 2336–2370 | `TRACER` | 33 | `CALL mpi_bcast (DEF_TRACER_USE_FRACTIONATION            ,1   ,mpi_logical   ,p_address_mas` | |

## mksrfdata/Aggregation_LakeSoilC.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 3–201 | `(defined TRACER) && (defined BGC)` | 197 | `SUBROUTINE Aggregation_LakeSoilC ( &` | |

## mksrfdata/Aggregation_MethanePH.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 3–566 | `(defined TRACER) && (defined BGC)` | 562 | `SUBROUTINE Aggregation_MethanePH (dir_rawdata, dir_model_landdata, lc_year)` | |

## mksrfdata/MKSRFDATA.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 80–82 | `(defined TRACER) && (defined BGC)` | 1 | `USE MOD_Tracer_Reactive_Methane_Preprocessing, only: methane_preprocessing_requirements` | |
| 113–115 | `(defined TRACER) && (defined BGC)` | 1 | `logical   :: requires_lake_soilc, requires_spatial_ph` | |
| 520–526 | `(defined TRACER) && (defined BGC)` | 5 | `CALL methane_preprocessing_requirements (requires_lake_soilc, requires_spatial_ph)` | |

## mksrfdata/MOD_RegionClip.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 40–42 | `(defined TRACER) && (defined BGC)` | 1 | `logical :: fexists` | |
| 344–364 | `(defined TRACER) && (defined BGC)` | 19 | `write(cyear,'(i4.4)') DEF_LC_YEAR` | |

## mkinidata/MOD_Initialize.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 40–42 | `defined(TRACER) && defined(BGC)` | 1 | `USE MOD_BGC_CNSummary, only: CNDriverSummarizeNonvegetatedSoilStates` | |
| 118–123 | `(defined TRACER) && (defined BGC)` | 4 | `real(r8) :: wetland_tot_c` | |
| 404–419 | `(defined TRACER) && (defined BGC)` | 14 | `write(cyear,'(i4.4)') lc_year` | |
| 1081–1085 | `TRACER` | 3 | `IF (patchtype(i) == 0 .or. patchtype(i) == 2)THEN` | |
| 1109–1148 | `TRACER` | 38 | `IF (patchtype(i) == 2) THEN` | |
| 1510–1512 | `defined(TRACER) && defined(BGC)` | 1 | `IF (patchtype(i) == 2) CALL CNDriverSummarizeNonvegetatedSoilStates(i, nl_soil, dz_soi, nd` | |

## extends/interception/MOD_Thermal_CanopyPhase_Extended.F90

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 1331–1334 | `TRACER` | 2 | `,qphs_thaw_lay = qphs_thaw_lay_th, &` | |

## include/define.h

| 行 | 条件 | 行数 | 第一行代码 | 状态 |
| --- | --- | --- | --- | --- |
| 104–106 | `(defined TRACER) && (defined Campbell_SOIL_MODEL)` | 1 | `` | |
| 107–111 | `(defined TRACER) && (defined BGC)` | 3 | `` | |

## 块外改动（不在 `#ifdef TRACER` 里，脚本抓不到，手工补）

| 位置 | 内容 | 去向 | 状态 |
| --- | --- | --- | --- |
| share/MOD_Namelist.F90:447–475 | 29 个 `DEF_TRACER_*` 声明不在 `#ifdef` 里（读/广播在） | P2 换成新 namelist，声明也放进 `#ifdef ConstituentTransport` | |
| main/MOD_Vars_1DAccFluxes.F90:18–19 | `interface history_acc_write` | P1 移进 `MOD_HistAccRestart` | |
| main/MOD_Vars_1DAccFluxes.F90:3525 | `#include <land_history_restart.inc>` | P1 删掉，换成独立模块 | |
| main/CoLM.F90:380 | `read_history_acc_restart`（初始化） | P1 | |
| main/CoLM.F90:451 | `restore_river_history_acc_restart`（`GridRiverLakeFlow` 内） | P1 | |
| main/CoLM.F90:735, 748, 749 | `write_` → `WRITE_TimeVariables` → `mark_` → `complete_`（提交协议） | P1：尽量压成 restart 前后各一个调用 | |
| main/MOD_Hist.F90:262, 298 | 写 history 时同步写侧车 | P1 | |
| include/define.h:103 | `#undef TRACER` | P1 换宏 | |

注：计划原来写“上游文件只加两个调用”，实际侧车有 7 个调用点（CoLM.F90 5 个、MOD_Hist 2 个）。P1 设计时要么压缩调用点，要么在登记表里逐个写理由。

另外，fork 宿主里不带 `#ifdef` 的签名或逻辑改动（例如为了传组分参数改过的接口），要在 P2 用“fork 宿主 vs 上游宿主”的完整 diff 再过一遍，不能只靠本清单。
