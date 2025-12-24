# GridRiverLakeSediment 模块实施计划

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** 将 CaMa-Flood 泥沙模块移植到 GridRiverLakeFlow 框架，实现完整的泥沙输运模拟。

**Architecture:** 创建独立的 `MOD_Grid_RiverLakeSediment.F90` 模块，包含泥沙状态变量、物理计算和诊断输出。通过预处理开关 `GridRiverLakeSediment` 控制编译，通过运行时开关 `DEF_USE_SEDIMENT` 控制启用。泥沙计算在水流计算完成后执行，使用独立的时间子循环。

**Tech Stack:** Fortran 90, NetCDF, MPI (可选)

**Reference Design:** `docs/plans/2025-12-24-grid-riverlake-sediment-design.md`

**Reference Implementation:** `extends/CaMa/src/cmf_ctrl_sed_mod.F90`

---

## Phase 1: 框架搭建

### Task 1.1: 添加预处理开关

**Files:**
- Modify: `include/define.h`

**Step 1: 查看现有预处理开关格式**

```bash
grep -n "GridRiverLake" include/define.h
```

**Step 2: 添加泥沙模块预处理开关**

在 `#define GridRiverLakeFlow` 下方添加：

```fortran
#define GridRiverLakeSediment
```

**Step 3: 提交**

```bash
git add include/define.h
git commit -m "feat: add GridRiverLakeSediment preprocessor switch"
```

---

### Task 1.2: 添加 Namelist 参数

**Files:**
- Modify: `main/MOD_Namelist.F90`

**Step 1: 查找 namelist 变量声明位置**

```bash
grep -n "DEF_Reservoir_Method" main/MOD_Namelist.F90 | head -5
```

**Step 2: 添加泥沙相关 namelist 变量声明**

在水文相关变量声明区域添加：

```fortran
   ! Sediment module parameters
   logical,  public :: DEF_USE_SEDIMENT     = .false.   ! Enable sediment module
   real(r8), public :: DEF_SED_LAMBDA       = 0.4_r8    ! Porosity [-]
   real(r8), public :: DEF_SED_LYRDPH       = 0.00005_r8 ! Active layer depth [m]
   real(r8), public :: DEF_SED_DENSITY      = 2.65_r8   ! Sediment density [g/cm3]
   real(r8), public :: DEF_SED_WATER_DENSITY = 1.0_r8   ! Water density [g/cm3]
   real(r8), public :: DEF_SED_VISKIN       = 1.0e-6_r8 ! Kinematic viscosity [m2/s]
   real(r8), public :: DEF_SED_VONKAR       = 0.4_r8    ! Von Karman coefficient [-]
   integer,  public :: DEF_SED_TOTLYRNUM    = 5         ! Number of deposition layers
   real(r8), public :: DEF_SED_DT_MAX       = 3600._r8  ! Max sediment timestep [s]
   character(len=256), public :: DEF_SED_DIAMETER = "0.0002,0.002,0.02" ! Grain diameters [m]
   ! Sediment yield parameters (Sunada & Hasegawa 1993)
   real(r8), public :: DEF_SED_PYLD         = 0.01_r8   ! Yield coefficient
   real(r8), public :: DEF_SED_PYLDC        = 2.0_r8    ! Slope exponent
   real(r8), public :: DEF_SED_PYLDPC       = 2.0_r8    ! Precipitation exponent
   real(r8), public :: DEF_SED_DSYLUNIT     = 1.0e-6_r8 ! Unit conversion factor
```

**Step 3: 查找 namelist 读取位置**

```bash
grep -n "DEF_Reservoir_Method" main/MOD_Namelist.F90 | grep -i "read\|nml"
```

**Step 4: 添加 namelist 读取**

在相应的 namelist 组中添加读取逻辑（根据项目具体格式调整）。

**Step 5: 提交**

```bash
git add main/MOD_Namelist.F90
git commit -m "feat: add sediment module namelist parameters"
```

---

### Task 1.3: 添加历史输出控制变量

**Files:**
- Modify: `main/MOD_Namelist.F90` (DEF_hist_vars 部分)

**Step 1: 查找历史输出变量定义位置**

```bash
grep -n "discharge\|riv_height" main/MOD_Namelist.F90 | head -10
```

**Step 2: 添加泥沙输出控制变量**

在 `DEF_hist_vars` 类型定义中添加：

```fortran
   ! Sediment output variables
   logical :: sedcon    = .true.    ! Suspended sediment concentration
   logical :: sedout    = .true.    ! Suspended sediment flux
   logical :: bedout    = .true.    ! Bedload flux
   logical :: sedinp    = .true.    ! Erosion/sediment input
   logical :: netflw    = .true.    ! Net exchange flux
   logical :: sedlayer  = .true.    ! Active layer storage
   logical :: shearvel  = .false.   ! Shear velocity (diagnostic)
   logical :: critshear = .false.   ! Critical shear velocity (diagnostic)
```

**Step 3: 提交**

```bash
git add main/MOD_Namelist.F90
git commit -m "feat: add sediment history output control variables"
```

---

### Task 1.4: 创建泥沙模块骨架

**Files:**
- Create: `main/HYDRO/MOD_Grid_RiverLakeSediment.F90`

**Step 1: 创建模块基本结构**

```fortran
#include <define.h>

#ifdef GridRiverLakeSediment
MODULE MOD_Grid_RiverLakeSediment
!-------------------------------------------------------------------------------------
! DESCRIPTION:
!
!   Sediment transport module for GridRiverLakeFlow.
!   Ported from CaMa-Flood sediment module (cmf_ctrl_sed_mod.F90)
!
! Created by: [Your Name], Dec 2025
!-------------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   IMPLICIT NONE

   !-------------------------------------------------------------------------------------
   ! Module Parameters
   !-------------------------------------------------------------------------------------
   integer,  save :: nsed           ! Number of sediment size classes
   integer,  save :: totlyrnum      ! Number of deposition layers
   integer,  save :: nlfp_sed       ! Number of floodplain layers for slope

   real(r8), save :: lambda         ! Porosity [-]
   real(r8), save :: lyrdph         ! Active layer depth [m]
   real(r8), save :: psedD          ! Sediment density [g/cm3]
   real(r8), save :: pwatD          ! Water density [g/cm3]
   real(r8), save :: visKin         ! Kinematic viscosity [m2/s]
   real(r8), save :: vonKar         ! Von Karman coefficient [-]

   ! Sediment yield parameters
   real(r8), save :: pyld           ! Yield coefficient
   real(r8), save :: pyldc          ! Slope exponent
   real(r8), save :: pyldpc         ! Precipitation exponent
   real(r8), save :: dsylunit       ! Unit conversion factor

   !-------------------------------------------------------------------------------------
   ! Static Data (read from DEF_UnitCatchment_file)
   !-------------------------------------------------------------------------------------
   real(r8), allocatable :: sed_frc   (:,:)    ! Sediment fraction [numucat, nsed]
   real(r8), allocatable :: sed_slope (:,:)    ! Floodplain slope [numucat, nlfp]
   real(r8), allocatable :: sDiam     (:)      ! Grain diameter [nsed]
   real(r8), allocatable :: setvel    (:)      ! Settling velocity [nsed]

   !-------------------------------------------------------------------------------------
   ! State Variables
   !-------------------------------------------------------------------------------------
   real(r8), allocatable :: sedcon  (:,:)      ! Suspended sediment concentration [numucat, nsed]
   real(r8), allocatable :: layer   (:,:)      ! Active layer storage [numucat, nsed]
   real(r8), allocatable :: seddep  (:,:,:)    ! Deposition layer storage [numucat, totlyrnum, nsed]

   !-------------------------------------------------------------------------------------
   ! Diagnostic Variables
   !-------------------------------------------------------------------------------------
   real(r8), allocatable :: sedout  (:,:)      ! Suspended sediment outflow [numucat, nsed]
   real(r8), allocatable :: bedout  (:,:)      ! Bedload outflow [numucat, nsed]
   real(r8), allocatable :: sedinp  (:,:)      ! Erosion input [numucat, nsed]
   real(r8), allocatable :: netflw  (:,:)      ! Net exchange flux [numucat, nsed]
   real(r8), allocatable :: shearvel(:)        ! Shear velocity [numucat]
   real(r8), allocatable :: critshearvel(:,:)  ! Critical shear velocity [numucat, nsed]
   real(r8), allocatable :: susvel  (:,:)      ! Suspension velocity [numucat, nsed]

   !-------------------------------------------------------------------------------------
   ! Accumulated Variables for Sediment Time-stepping
   !-------------------------------------------------------------------------------------
   real(r8), save :: sed_acc_time              ! Accumulated time for averaging
   real(r8), allocatable :: sed_acc_veloc(:)   ! Accumulated velocity [numucat]
   real(r8), allocatable :: sed_acc_wdsrf(:)   ! Accumulated water depth [numucat]
   real(r8), allocatable :: sed_precip(:)      ! Precipitation for sediment yield [numucat]

   !-------------------------------------------------------------------------------------
   ! Accumulated Variables for History Output
   !-------------------------------------------------------------------------------------
   real(r8), allocatable :: a_sedcon  (:,:)    ! Accumulated sedcon
   real(r8), allocatable :: a_sedout  (:,:)    ! Accumulated sedout
   real(r8), allocatable :: a_bedout  (:,:)    ! Accumulated bedout
   real(r8), allocatable :: a_sedinp  (:,:)    ! Accumulated sedinp
   real(r8), allocatable :: a_netflw  (:,:)    ! Accumulated netflw
   real(r8), allocatable :: a_layer   (:,:)    ! Accumulated layer
   real(r8), allocatable :: a_shearvel(:)      ! Accumulated shearvel

   !-------------------------------------------------------------------------------------
   ! Public Subroutines
   !-------------------------------------------------------------------------------------
   PUBLIC :: grid_sediment_init
   PUBLIC :: grid_sediment_calc
   PUBLIC :: grid_sediment_final
   PUBLIC :: sediment_diag_accumulate
   PUBLIC :: sediment_forcing_put

CONTAINS

   !-------------------------------------------------------------------------------------
   SUBROUTINE grid_sediment_init()
   ! Initialize sediment module
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE
      ! Placeholder - to be implemented
      WRITE(*,*) 'MOD_Grid_RiverLakeSediment: grid_sediment_init called'
   END SUBROUTINE grid_sediment_init

   !-------------------------------------------------------------------------------------
   SUBROUTINE grid_sediment_calc(deltime)
   ! Main sediment calculation routine
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE
   real(r8), intent(in) :: deltime
      ! Placeholder - to be implemented
      WRITE(*,*) 'MOD_Grid_RiverLakeSediment: grid_sediment_calc called, dt=', deltime
   END SUBROUTINE grid_sediment_calc

   !-------------------------------------------------------------------------------------
   SUBROUTINE sediment_diag_accumulate(dt, veloc, wdsrf)
   ! Accumulate water flow variables for sediment calculation
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE
   real(r8), intent(in) :: dt
   real(r8), intent(in) :: veloc(:)
   real(r8), intent(in) :: wdsrf(:)
      ! Placeholder - to be implemented
   END SUBROUTINE sediment_diag_accumulate

   !-------------------------------------------------------------------------------------
   SUBROUTINE sediment_forcing_put(precip)
   ! Store precipitation for sediment yield calculation
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE
   real(r8), intent(in) :: precip(:)
      ! Placeholder - to be implemented
   END SUBROUTINE sediment_forcing_put

   !-------------------------------------------------------------------------------------
   SUBROUTINE grid_sediment_final()
   ! Cleanup sediment module
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE
      IF (allocated(sed_frc     )) deallocate(sed_frc     )
      IF (allocated(sed_slope   )) deallocate(sed_slope   )
      IF (allocated(sDiam       )) deallocate(sDiam       )
      IF (allocated(setvel      )) deallocate(setvel      )
      IF (allocated(sedcon      )) deallocate(sedcon      )
      IF (allocated(layer       )) deallocate(layer       )
      IF (allocated(seddep      )) deallocate(seddep      )
      IF (allocated(sedout      )) deallocate(sedout      )
      IF (allocated(bedout      )) deallocate(bedout      )
      IF (allocated(sedinp      )) deallocate(sedinp      )
      IF (allocated(netflw      )) deallocate(netflw      )
      IF (allocated(shearvel    )) deallocate(shearvel    )
      IF (allocated(critshearvel)) deallocate(critshearvel)
      IF (allocated(susvel      )) deallocate(susvel      )
      IF (allocated(sed_acc_veloc)) deallocate(sed_acc_veloc)
      IF (allocated(sed_acc_wdsrf)) deallocate(sed_acc_wdsrf)
      IF (allocated(sed_precip  )) deallocate(sed_precip  )
      IF (allocated(a_sedcon    )) deallocate(a_sedcon    )
      IF (allocated(a_sedout    )) deallocate(a_sedout    )
      IF (allocated(a_bedout    )) deallocate(a_bedout    )
      IF (allocated(a_sedinp    )) deallocate(a_sedinp    )
      IF (allocated(a_netflw    )) deallocate(a_netflw    )
      IF (allocated(a_layer     )) deallocate(a_layer     )
      IF (allocated(a_shearvel  )) deallocate(a_shearvel  )
   END SUBROUTINE grid_sediment_final

END MODULE MOD_Grid_RiverLakeSediment
#endif
```

**Step 2: 提交**

```bash
git add main/HYDRO/MOD_Grid_RiverLakeSediment.F90
git commit -m "feat: create MOD_Grid_RiverLakeSediment module skeleton"
```

---

### Task 1.5: 添加模块到 Makefile

**Files:**
- Modify: `Makefile`

**Step 1: 查找 HYDRO 模块编译规则**

```bash
grep -n "MOD_Grid_RiverLake" Makefile | head -10
```

**Step 2: 添加泥沙模块编译规则**

在相应位置添加：

```makefile
$(OBJS_DIR)/MOD_Grid_RiverLakeSediment.o: main/HYDRO/MOD_Grid_RiverLakeSediment.F90
	$(FC) -c $(FOPTS) $(INCLUDE_DIR) -o $@ $<
```

并将 `$(OBJS_DIR)/MOD_Grid_RiverLakeSediment.o` 添加到依赖列表中。

**Step 3: 提交**

```bash
git add Makefile
git commit -m "build: add MOD_Grid_RiverLakeSediment to Makefile"
```

---

## Phase 2: 静态数据读取

### Task 2.1: 实现维度读取

**Files:**
- Modify: `main/HYDRO/MOD_Grid_RiverLakeSediment.F90`

**Step 1: 添加维度读取子程序**

在 `grid_sediment_init` 中添加：

```fortran
   SUBROUTINE grid_sediment_init()
   USE MOD_NetCDFSerial
   USE MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, &
      ucat_data_address, topo_rivwth, topo_rivlen
   IMPLICIT NONE

   character(len=256) :: parafile
   integer :: i

      IF (.not. DEF_USE_SEDIMENT) RETURN

      ! Set parameters from namelist
      lambda    = DEF_SED_LAMBDA
      lyrdph    = DEF_SED_LYRDPH
      psedD     = DEF_SED_DENSITY
      pwatD     = DEF_SED_WATER_DENSITY
      visKin    = DEF_SED_VISKIN
      vonKar    = DEF_SED_VONKAR
      totlyrnum = DEF_SED_TOTLYRNUM
      pyld      = DEF_SED_PYLD
      pyldc     = DEF_SED_PYLDC
      pyldpc    = DEF_SED_PYLDPC
      dsylunit  = DEF_SED_DSYLUNIT

      parafile = DEF_UnitCatchment_file

      ! Read dimensions from NetCDF file
      IF (p_is_master) THEN
         CALL ncio_inquire_length(parafile, 'sed_n', nsed)
         CALL ncio_inquire_length(parafile, 'slope_layers', nlfp_sed)
      ENDIF

#ifdef USEMPI
      CALL mpi_bcast(nsed, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast(nlfp_sed, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif

      WRITE(*,*) 'Sediment module: nsed=', nsed, ' nlfp_sed=', nlfp_sed

      ! Parse grain diameters from string
      CALL parse_grain_diameters()

      ! Calculate settling velocities
      CALL calc_settling_velocities()

      ! Read static data
      CALL read_sediment_static_data(parafile)

      ! Allocate and initialize state variables
      CALL allocate_sediment_vars()

      ! Initialize state from sed_frc
      CALL initialize_sediment_state()

   END SUBROUTINE grid_sediment_init
```

**Step 2: 提交**

```bash
git add main/HYDRO/MOD_Grid_RiverLakeSediment.F90
git commit -m "feat: implement sediment dimension reading"
```

---

### Task 2.2: 实现粒径解析

**Files:**
- Modify: `main/HYDRO/MOD_Grid_RiverLakeSediment.F90`

**Step 1: 添加粒径解析子程序**

```fortran
   SUBROUTINE parse_grain_diameters()
   ! Parse grain diameters from comma-separated string
   IMPLICIT NONE
   character(len=256) :: str
   integer :: i, j, k, n

      allocate(sDiam(nsed))
      str = trim(adjustl(DEF_SED_DIAMETER))

      n = 0
      j = 1
      DO i = 1, len_trim(str)
         IF (str(i:i) == ',' .or. i == len_trim(str)) THEN
            n = n + 1
            IF (i == len_trim(str)) THEN
               k = i
            ELSE
               k = i - 1
            ENDIF
            IF (n <= nsed) THEN
               read(str(j:k), *) sDiam(n)
            ENDIF
            j = i + 1
         ENDIF
      ENDDO

      IF (n /= nsed) THEN
         WRITE(*,*) 'WARNING: Number of diameters does not match nsed'
         WRITE(*,*) '  Parsed:', n, ' Expected:', nsed
      ENDIF

      WRITE(*,*) 'Grain diameters (m):', sDiam

   END SUBROUTINE parse_grain_diameters
```

**Step 2: 提交**

```bash
git add main/HYDRO/MOD_Grid_RiverLakeSediment.F90
git commit -m "feat: implement grain diameter parsing"
```

---

### Task 2.3: 实现沉降速度计算

**Files:**
- Modify: `main/HYDRO/MOD_Grid_RiverLakeSediment.F90`

**Step 1: 添加沉降速度计算子程序**

从 `cmf_ctrl_sed_mod.F90` 移植：

```fortran
   SUBROUTINE calc_settling_velocities()
   ! Calculate settling velocity using Stokes-Rubey formula
   USE MOD_Const_Physical, only: grav
   IMPLICIT NONE
   real(r8) :: sTmp
   integer :: i

      allocate(setvel(nsed))

      DO i = 1, nsed
         sTmp = 6.0_r8 * visKin / sDiam(i)
         setvel(i) = sqrt(2.0_r8/3.0_r8 * (psedD-pwatD)/pwatD * grav * sDiam(i) &
                    + sTmp*sTmp) - sTmp
      ENDDO

      WRITE(*,*) 'Settling velocities (m/s):', setvel

   END SUBROUTINE calc_settling_velocities
```

**Step 2: 提交**

```bash
git add main/HYDRO/MOD_Grid_RiverLakeSediment.F90
git commit -m "feat: implement settling velocity calculation"
```

---

### Task 2.4: 实现静态数据读取

**Files:**
- Modify: `main/HYDRO/MOD_Grid_RiverLakeSediment.F90`

**Step 1: 添加静态数据读取子程序**

使用与 `MOD_Grid_RiverLakeNetwork` 相同的模式：

```fortran
   SUBROUTINE read_sediment_static_data(parafile)
   USE MOD_Grid_RiverLakeNetwork, only: numucat, readin_riverlake_parameter
   IMPLICIT NONE
   character(len=*), intent(in) :: parafile

      IF (p_is_worker) THEN
         IF (numucat > 0) THEN
            allocate(sed_frc  (nsed, numucat))
            allocate(sed_slope(nlfp_sed, numucat))
         ENDIF
      ENDIF

      ! Read sediment fraction
      CALL readin_riverlake_parameter(parafile, 'sed_frc', rdata2d=sed_frc)

      ! Read sediment slope
      CALL readin_riverlake_parameter(parafile, 'sed_slope', rdata2d=sed_slope)

      ! Normalize sed_frc (ensure sum = 1)
      CALL normalize_sed_frc()

   END SUBROUTINE read_sediment_static_data

   SUBROUTINE normalize_sed_frc()
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE
   integer :: i
   real(r8) :: frc_sum

      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      DO i = 1, numucat
         frc_sum = sum(sed_frc(:,i))
         IF (frc_sum > 0._r8) THEN
            sed_frc(:,i) = sed_frc(:,i) / frc_sum
         ELSE
            sed_frc(:,i) = 1._r8 / real(nsed, r8)
         ENDIF
      ENDDO

   END SUBROUTINE normalize_sed_frc
```

**Step 2: 提交**

```bash
git add main/HYDRO/MOD_Grid_RiverLakeSediment.F90
git commit -m "feat: implement sediment static data reading"
```

---

### Task 2.5: 实现变量分配与初始化

**Files:**
- Modify: `main/HYDRO/MOD_Grid_RiverLakeSediment.F90`

**Step 1: 添加变量分配子程序**

```fortran
   SUBROUTINE allocate_sediment_vars()
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE

      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      ! State variables
      allocate(sedcon(nsed, numucat))
      allocate(layer (nsed, numucat))
      allocate(seddep(nsed, totlyrnum, numucat))

      ! Diagnostic variables
      allocate(sedout      (nsed, numucat))
      allocate(bedout      (nsed, numucat))
      allocate(sedinp      (nsed, numucat))
      allocate(netflw      (nsed, numucat))
      allocate(shearvel    (numucat))
      allocate(critshearvel(nsed, numucat))
      allocate(susvel      (nsed, numucat))

      ! Accumulation variables
      allocate(sed_acc_veloc(numucat))
      allocate(sed_acc_wdsrf(numucat))
      allocate(sed_precip   (numucat))

      ! History output variables
      allocate(a_sedcon  (nsed, numucat))
      allocate(a_sedout  (nsed, numucat))
      allocate(a_bedout  (nsed, numucat))
      allocate(a_sedinp  (nsed, numucat))
      allocate(a_netflw  (nsed, numucat))
      allocate(a_layer   (nsed, numucat))
      allocate(a_shearvel(numucat))

      ! Initialize to zero
      sedcon       = 0._r8
      layer        = 0._r8
      seddep       = 0._r8
      sedout       = 0._r8
      bedout       = 0._r8
      sedinp       = 0._r8
      netflw       = 0._r8
      shearvel     = 0._r8
      critshearvel = 0._r8
      susvel       = 0._r8
      sed_acc_veloc = 0._r8
      sed_acc_wdsrf = 0._r8
      sed_acc_time  = 0._r8
      sed_precip    = 0._r8
      a_sedcon     = 0._r8
      a_sedout     = 0._r8
      a_bedout     = 0._r8
      a_sedinp     = 0._r8
      a_netflw     = 0._r8
      a_layer      = 0._r8
      a_shearvel   = 0._r8

   END SUBROUTINE allocate_sediment_vars

   SUBROUTINE initialize_sediment_state()
   USE MOD_Grid_RiverLakeNetwork, only: numucat, topo_rivwth, topo_rivlen
   IMPLICIT NONE
   integer :: i, ilyr

      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      ! Initialize active layer based on sed_frc
      DO i = 1, numucat
         layer(:,i) = lyrdph * topo_rivwth(i) * topo_rivlen(i) * sed_frc(:,i)

         ! Initialize deposition layers
         DO ilyr = 1, totlyrnum - 1
            seddep(:,ilyr,i) = layer(:,i)
         ENDDO
         ! Bottom layer gets extra depth
         seddep(:,totlyrnum,i) = max(10._r8 - lyrdph*totlyrnum, 0._r8) &
            * topo_rivwth(i) * topo_rivlen(i) * sed_frc(:,i)
      ENDDO

   END SUBROUTINE initialize_sediment_state
```

**Step 2: 提交**

```bash
git add main/HYDRO/MOD_Grid_RiverLakeSediment.F90
git commit -m "feat: implement sediment variable allocation and initialization"
```

---

## Phase 3: 核心物理计算

### Task 3.1: 实现剪切流速计算

**Files:**
- Modify: `main/HYDRO/MOD_Grid_RiverLakeSediment.F90`

**Step 1: 添加剪切流速计算函数**

```fortran
   FUNCTION calc_shear_velocity(rivvel, rivdph, rivman) RESULT(svel)
   ! Calculate shear velocity using Manning's equation
   USE MOD_Const_Physical, only: grav
   IMPLICIT NONE
   real(r8), intent(in) :: rivvel   ! River velocity [m/s]
   real(r8), intent(in) :: rivdph   ! River depth [m]
   real(r8), intent(in) :: rivman   ! Manning coefficient
   real(r8) :: svel

      IF (rivdph > 0._r8) THEN
         svel = sqrt(grav * rivman**2 * rivvel**2 * rivdph**(-1._r8/3._r8))
      ELSE
         svel = 0._r8
      ENDIF

   END FUNCTION calc_shear_velocity
```

**Step 2: 提交**

```bash
git add main/HYDRO/MOD_Grid_RiverLakeSediment.F90
git commit -m "feat: implement shear velocity calculation"
```

---

### Task 3.2: 实现临界剪切力计算

**Files:**
- Modify: `main/HYDRO/MOD_Grid_RiverLakeSediment.F90`

**Step 1: 添加临界剪切力计算函数**

```fortran
   FUNCTION calc_critical_shear_velocity(diam) RESULT(csvel)
   ! Calculate critical shear velocity using Shields curve
   IMPLICIT NONE
   real(r8), intent(in) :: diam    ! Grain diameter [m]
   real(r8) :: csvel               ! [(cm/s)^2]
   real(r8) :: cA, cB

      cB = 1._r8
      IF (diam >= 0.00303_r8) THEN
         cA = 80.9_r8
      ELSEIF (diam >= 0.00118_r8) THEN
         cA = 134.6_r8
         cB = 31._r8 / 32._r8
      ELSEIF (diam >= 0.000565_r8) THEN
         cA = 55._r8
      ELSEIF (diam >= 0.000065_r8) THEN
         cA = 8.41_r8
         cB = 11._r8 / 32._r8
      ELSE
         cA = 226._r8
      ENDIF

      csvel = cA * (diam * 100._r8) ** cB

   END FUNCTION calc_critical_shear_velocity

   SUBROUTINE calc_critical_shear_egiazoroff(i, svel, csvel_out)
   ! Calculate critical shear velocity using Egiazoroff equation for mixed-size sediment
   IMPLICIT NONE
   integer,  intent(in)  :: i           ! Unit catchment index
   real(r8), intent(in)  :: svel        ! Shear velocity [m/s]
   real(r8), intent(out) :: csvel_out(nsed)

   real(r8) :: dMean, csVel0
   integer  :: ised
   real(r8) :: layer_sum

      layer_sum = sum(layer(:,i))

      IF (layer_sum <= 0._r8) THEN
         csvel_out(:) = 1.e20_r8
         RETURN
      ENDIF

      ! Calculate mean diameter
      dMean = 0._r8
      DO ised = 1, nsed
         dMean = dMean + sDiam(ised) * layer(ised,i) / layer_sum
      ENDDO

      csVel0 = calc_critical_shear_velocity(dMean)

      DO ised = 1, nsed
         IF (sDiam(ised) / dMean >= 0.4_r8) THEN
            csvel_out(ised) = sqrt(csVel0 * sDiam(ised) / dMean) * &
               (log10(19._r8) / log10(19._r8 * sDiam(ised) / dMean)) * 0.01_r8
         ELSE
            csvel_out(ised) = sqrt(0.85_r8 * csVel0) * 0.01_r8
         ENDIF
      ENDDO

   END SUBROUTINE calc_critical_shear_egiazoroff
```

**Step 2: 提交**

```bash
git add main/HYDRO/MOD_Grid_RiverLakeSediment.F90
git commit -m "feat: implement critical shear velocity calculation"
```

---

### Task 3.3: 实现悬浮速度计算

**Files:**
- Modify: `main/HYDRO/MOD_Grid_RiverLakeSediment.F90`

**Step 1: 添加悬浮速度计算函数**

```fortran
   SUBROUTINE calc_suspend_velocity(csvel, svel, susvel_out)
   ! Calculate suspension velocity using Uchida & Fukuoka (2019) Eq.44
   IMPLICIT NONE
   real(r8), intent(in)  :: csvel(nsed)  ! Critical shear velocity [m/s]
   real(r8), intent(in)  :: svel         ! Shear velocity [m/s]
   real(r8), intent(out) :: susvel_out(nsed)

   real(r8) :: alpha, a, cB, sTmp
   integer  :: ised

      alpha = vonKar / 6._r8
      a = 0.08_r8
      cB = 1._r8 - lambda

      susvel_out(:) = 0._r8

      DO ised = 1, nsed
         IF (csvel(ised) > svel) CYCLE
         IF (svel <= 0._r8) CYCLE

         sTmp = setvel(ised) / alpha / svel
         susvel_out(ised) = max(setvel(ised) * cB / (1._r8 + sTmp) * &
            (1._r8 - a*sTmp) / (1._r8 + (1._r8-a)*sTmp), 0._r8)
      ENDDO

   END SUBROUTINE calc_suspend_velocity
```

**Step 2: 提交**

```bash
git add main/HYDRO/MOD_Grid_RiverLakeSediment.F90
git commit -m "feat: implement suspension velocity calculation"
```

---

### Task 3.4: 实现平流输运计算

**Files:**
- Modify: `main/HYDRO/MOD_Grid_RiverLakeSediment.F90`

**Step 1: 添加平流输运计算子程序**

这是最复杂的部分，从 `cmf_ctrl_sed_mod.F90` 的 `calc_advection` 移植：

```fortran
   SUBROUTINE calc_sediment_advection(dt, rivout, rivsto)
   ! Calculate suspended sediment and bedload advection
   USE MOD_Const_Physical, only: grav
   USE MOD_Grid_RiverLakeNetwork, only: numucat, ucat_next, topo_rivwth
   USE MOD_WorkerPushData
   IMPLICIT NONE

   real(r8), intent(in) :: dt           ! Time step [s]
   real(r8), intent(in) :: rivout(:)    ! River outflow [m3/s]
   real(r8), intent(in) :: rivsto(:)    ! River storage [m3]

   real(r8), allocatable :: sedsto(:,:)       ! Sediment storage [numucat, nsed]
   real(r8), allocatable :: bOut(:,:), sOut(:,:)
   real(r8), allocatable :: brate(:,:), srate(:,:)
   real(r8), allocatable :: sedcon_next(:,:)

   integer  :: i, ised, i0, i1
   real(r8) :: plusVel, minusVel
   real(r8) :: layer_sum

      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      allocate(sedsto(nsed, numucat))
      allocate(bOut  (nsed, numucat))
      allocate(sOut  (nsed, numucat))
      allocate(brate (nsed, numucat))
      allocate(srate (nsed, numucat))
      allocate(sedcon_next(nsed, numucat))

      ! Calculate sediment storage from concentration
      DO i = 1, numucat
         sedsto(:,i) = sedcon(:,i) * max(rivsto(i), 0._r8)
      ENDDO

      ! Get downstream sediment concentration
      DO ised = 1, nsed
         CALL worker_push_data(push_next2ucat, sedcon(ised,:), sedcon_next(ised,:), &
            fillvalue = 0._r8)
      ENDDO

      bOut(:,:) = 0._r8
      sOut(:,:) = 0._r8

      ! Calculate outflows
      DO i = 1, numucat
         IF (rivout(i) >= 0._r8) THEN
            i0 = i
            i1 = ucat_next(i)
         ELSE
            i0 = ucat_next(i)
            i1 = i
         ENDIF

         IF (rivout(i) == 0._r8) THEN
            sedout(:,i) = 0._r8
            bedout(:,i) = 0._r8
            CYCLE
         ENDIF

         ! Suspended sediment outflow
         IF (i0 < 0) THEN
            sedout(:,i) = sedcon(:,i1) * rivout(i)
         ELSE
            sedout(:,i) = sedcon(:,i0) * rivout(i)
            sOut(:,i0) = sOut(:,i0) + abs(sedout(:,i)) * dt
         ENDIF

         ! Bedload outflow
         layer_sum = sum(layer(:,i))
         IF (all(critshearvel(:,i) >= shearvel(i)) .or. layer_sum == 0._r8 .or. i0 < 0) THEN
            bedout(:,i) = 0._r8
         ELSE
            DO ised = 1, nsed
               IF (critshearvel(ised,i) >= shearvel(i) .or. layer(ised,i) == 0._r8) THEN
                  bedout(ised,i) = 0._r8
                  CYCLE
               ENDIF
               plusVel = shearvel(i) + critshearvel(ised,i)
               minusVel = shearvel(i) - critshearvel(ised,i)
               bedout(ised,i) = 17._r8 * topo_rivwth(i) * plusVel * minusVel * minusVel &
                  / ((psedD-pwatD)/pwatD) / grav * layer(ised,i) / layer_sum
               IF (i0 > 0) bOut(ised,i0) = bOut(ised,i0) + bedout(ised,i) * dt
            ENDDO
         ENDIF
      ENDDO

      ! Adjust outflow if larger than available storage
      brate(:,:) = 1._r8
      srate(:,:) = 1._r8

      DO i = 1, numucat
         DO ised = 1, nsed
            IF (sOut(ised,i) > 1.e-8_r8) THEN
               srate(ised,i) = min(sedsto(ised,i) / sOut(ised,i), 1._r8)
            ENDIF
            IF (bOut(ised,i) > 1.e-8_r8) THEN
               brate(ised,i) = min(layer(ised,i) / bOut(ised,i), 1._r8)
            ENDIF
         ENDDO
      ENDDO

      ! Apply adjusted outflows and update storage
      DO i = 1, numucat
         IF (rivout(i) >= 0._r8) THEN
            i0 = i
            i1 = ucat_next(i)
         ELSE
            i0 = ucat_next(i)
            i1 = i
         ENDIF

         IF (i0 > 0) THEN
            sedout(:,i) = sedout(:,i) * srate(:,i0)
            sedsto(:,i0) = max(sedsto(:,i0) - abs(sedout(:,i)) * dt, 0._r8)
            bedout(:,i) = bedout(:,i) * brate(:,i0)
            layer(:,i0) = max(layer(:,i0) - abs(bedout(:,i)) * dt, 0._r8)
         ENDIF

         IF (i1 > 0) THEN
            sedsto(:,i1) = max(sedsto(:,i1) + abs(sedout(:,i)) * dt, 0._r8)
            layer(:,i1) = max(layer(:,i1) + abs(bedout(:,i)) * dt, 0._r8)
         ENDIF
      ENDDO

      ! Update concentration from storage
      DO i = 1, numucat
         IF (rivsto(i) > 0._r8) THEN
            sedcon(:,i) = sedsto(:,i) / rivsto(i)
         ENDIF
      ENDDO

      deallocate(sedsto, bOut, sOut, brate, srate, sedcon_next)

   END SUBROUTINE calc_sediment_advection
```

**Step 2: 提交**

```bash
git add main/HYDRO/MOD_Grid_RiverLakeSediment.F90
git commit -m "feat: implement sediment advection calculation"
```

---

### Task 3.5: 实现悬浮-沉积交换计算

**Files:**
- Modify: `main/HYDRO/MOD_Grid_RiverLakeSediment.F90`

**Step 1: 添加交换计算子程序**

```fortran
   SUBROUTINE calc_sediment_exchange(dt, rivsto, rivwth, rivlen)
   ! Calculate suspension-deposition exchange
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE

   real(r8), intent(in) :: dt
   real(r8), intent(in) :: rivsto(:)
   real(r8), intent(in) :: rivwth(:)
   real(r8), intent(in) :: rivlen(:)

   real(r8) :: Es(nsed), D(nsed), Zd(nsed)
   real(r8) :: sedsto(nsed), dTmp(nsed)
   real(r8) :: layer_sum, area, dTmp1
   integer  :: i, ised
   real(r8), parameter :: IGNORE_DPH = 0.05_r8

      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      DO i = 1, numucat
         IF (rivsto(i) < rivwth(i) * rivlen(i) * IGNORE_DPH) THEN
            netflw(:,i) = 0._r8
            CYCLE
         ENDIF

         layer_sum = sum(layer(:,i))
         area = rivwth(i) * rivlen(i)

         ! Calculate entrainment (suspension)
         IF (layer_sum == 0._r8 .or. all(susvel(:,i) == 0._r8)) THEN
            Es(:) = 0._r8
         ELSE
            Es(:) = susvel(:,i) * (1._r8 - lambda) * area * layer(:,i) / layer_sum
            Es(:) = max(Es(:), 0._r8)
         ENDIF

         ! Calculate deposition
         IF (shearvel(i) == 0._r8 .or. all(setvel(:) == 0._r8)) THEN
            D(:) = 0._r8
         ELSE
            DO ised = 1, nsed
               Zd(ised) = 6._r8 * setvel(ised) / vonKar / shearvel(i)
               D(ised) = setvel(ised) * area * sedcon(ised,i) * &
                  Zd(ised) / (1._r8 - exp(-Zd(ised)))
            ENDDO
            D(:) = max(D(:), 0._r8)
         ENDIF

         netflw(:,i) = Es(:) - D(:)

         ! Apply exchange with mass conservation
         sedsto(:) = sedcon(:,i) * rivsto(i)

         DO ised = 1, nsed
            IF (netflw(ised,i) == 0._r8) THEN
               CYCLE
            ELSEIF (netflw(ised,i) > 0._r8) THEN
               ! Suspension: transfer from layer to water
               dTmp1 = netflw(ised,i) * dt / (1._r8 - lambda)
               IF (dTmp1 < layer(ised,i)) THEN
                  layer(ised,i) = layer(ised,i) - dTmp1
               ELSE
                  netflw(ised,i) = layer(ised,i) * (1._r8 - lambda) / dt
                  layer(ised,i) = 0._r8
               ENDIF
               sedsto(ised) = sedsto(ised) + netflw(ised,i) * dt
            ELSE
               ! Deposition: transfer from water to layer
               IF (abs(netflw(ised,i)) * dt < sedsto(ised)) THEN
                  sedsto(ised) = max(sedsto(ised) - abs(netflw(ised,i)) * dt, 0._r8)
               ELSE
                  netflw(ised,i) = -sedsto(ised) / dt
                  sedsto(ised) = 0._r8
               ENDIF
               layer(ised,i) = layer(ised,i) + abs(netflw(ised,i)) * dt / (1._r8 - lambda)
            ENDIF
         ENDDO

         ! Add erosion input
         sedsto(:) = sedsto(:) + sedinp(:,i) * dt

         ! Limit concentration to 1%
         IF (sum(sedsto(:)) > rivsto(i) * 0.01_r8) THEN
            dTmp(:) = (sum(sedsto(:)) - rivsto(i) * 0.01_r8) * sedsto(:) / sum(sedsto(:))
            netflw(:,i) = netflw(:,i) - dTmp(:) / dt
            sedsto(:) = sedsto(:) - dTmp(:)
            layer(:,i) = layer(:,i) + dTmp(:) / (1._r8 - lambda)
         ENDIF

         ! Update concentration
         IF (rivsto(i) > 0._r8) THEN
            sedcon(:,i) = sedsto(:) / rivsto(i)
         ENDIF
      ENDDO

   END SUBROUTINE calc_sediment_exchange
```

**Step 2: 提交**

```bash
git add main/HYDRO/MOD_Grid_RiverLakeSediment.F90
git commit -m "feat: implement suspension-deposition exchange"
```

---

### Task 3.6: 实现沉积层重分布

**Files:**
- Modify: `main/HYDRO/MOD_Grid_RiverLakeSediment.F90`

**Step 1: 添加沉积层重分布子程序**

```fortran
   SUBROUTINE calc_layer_redistribution(rivwth, rivlen)
   ! Redistribute sediment into vertical bed layers
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE

   real(r8), intent(in) :: rivwth(:)
   real(r8), intent(in) :: rivlen(:)

   real(r8) :: lyrvol, diff
   real(r8) :: layerP(nsed), seddepP(totlyrnum+1, nsed), tmp(nsed)
   integer  :: i, ilyr, jlyr, slyr

      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      DO i = 1, numucat
         lyrvol = lyrdph * rivwth(i) * rivlen(i)

         ! Ensure non-negative values
         layer(:,i) = max(layer(:,i), 0._r8)
         seddep(:,:,i) = max(seddep(:,:,i), 0._r8)

         ! If total bed storage less than layer volume
         IF (sum(layer(:,i)) + sum(seddep(:,:,i)) <= lyrvol) THEN
            layer(:,i) = layer(:,i) + sum(seddep(:,:,i), dim=2)
            seddep(:,:,i) = 0._r8
            CYCLE
         ENDIF

         ! Distribute into top exchange layer
         layerP(:) = layer(:,i)
         IF (sum(layerP(:)) >= lyrvol) THEN
            layer(:,i) = layerP(:) * min(lyrvol / sum(layerP(:)), 1._r8)
            layerP(:) = max(layerP(:) - layer(:,i), 0._r8)
            slyr = 0
         ELSEIF (sum(seddep(:,:,i)) > 0._r8) THEN
            layerP(:) = 0._r8
            DO ilyr = 1, totlyrnum
               diff = lyrvol - sum(layer(:,i))
               IF (diff <= 0._r8) EXIT
               IF (sum(seddep(:,ilyr,i)) <= diff) THEN
                  layer(:,i) = layer(:,i) + seddep(:,ilyr,i)
                  seddep(:,ilyr,i) = 0._r8
                  slyr = ilyr + 1
               ELSE
                  tmp(:) = diff * seddep(:,ilyr,i) / sum(seddep(:,ilyr,i))
                  layer(:,i) = layer(:,i) + tmp(:)
                  seddep(:,ilyr,i) = max(seddep(:,ilyr,i) - tmp(:), 0._r8)
                  slyr = ilyr
                  EXIT
               ENDIF
            ENDDO
         ELSE
            seddep(:,:,i) = 0._r8
            CYCLE
         ENDIF

         IF (sum(seddep(:,:,i)) == 0._r8) CYCLE

         ! Distribute remaining bedload into vertical deposition layers
         seddepP(1,:) = layerP(:)
         seddepP(2:,:) = seddep(:,:,i)
         seddep(:,:,i) = 0._r8

         DO ilyr = 1, totlyrnum - 1
            IF (sum(seddep(:,ilyr,i)) == lyrvol) CYCLE
            DO jlyr = slyr + 1, totlyrnum + 1
               diff = lyrvol - sum(seddep(:,ilyr,i))
               IF (diff <= 0._r8) EXIT
               IF (sum(seddepP(jlyr,:)) <= diff) THEN
                  seddep(:,ilyr,i) = seddep(:,ilyr,i) + seddepP(jlyr,:)
                  seddepP(jlyr,:) = 0._r8
               ELSE
                  tmp(:) = diff * seddepP(jlyr,:) / sum(seddepP(jlyr,:))
                  seddep(:,ilyr,i) = seddep(:,ilyr,i) + tmp(:)
                  seddepP(jlyr,:) = max(seddepP(jlyr,:) - tmp(:), 0._r8)
                  EXIT
               ENDIF
            ENDDO
         ENDDO

         IF (sum(seddepP) > 0._r8) THEN
            seddep(:,totlyrnum,i) = sum(seddepP, dim=1)
         ENDIF
      ENDDO

   END SUBROUTINE calc_layer_redistribution
```

**Step 2: 提交**

```bash
git add main/HYDRO/MOD_Grid_RiverLakeSediment.F90
git commit -m "feat: implement layer redistribution"
```

---

### Task 3.7: 实现侵蚀产沙计算

**Files:**
- Modify: `main/HYDRO/MOD_Grid_RiverLakeSediment.F90`

**Step 1: 添加产沙计算子程序**

```fortran
   SUBROUTINE calc_sediment_yield(fldfrc, grarea)
   ! Calculate sediment yield from precipitation
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE

   real(r8), intent(in) :: fldfrc(:)   ! Flooded fraction
   real(r8), intent(in) :: grarea(:)   ! Grid area [m2]

   real(r8) :: precip_mm
   integer  :: i, ilyr

      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      sedinp(:,:) = 0._r8

      DO i = 1, numucat
         ! Convert precip from mm/s or kg/m2/s to mm/day for threshold
         precip_mm = sed_precip(i) * 86400._r8

         IF (precip_mm <= 10._r8) CYCLE

         ! Calculate erosion for each floodplain layer
         DO ilyr = 1, nlfp_sed
            IF (fldfrc(i) * nlfp_sed > real(ilyr, r8)) CYCLE  ! No erosion if submerged

            sedinp(:,i) = sedinp(:,i) + &
               pyld * (sed_precip(i) * 3600._r8)**pyldpc * sed_slope(ilyr,i)**pyldc / 3600._r8 &
               * grarea(i) * min(real(ilyr, r8)/real(nlfp_sed, r8) - fldfrc(i), 1._r8/real(nlfp_sed, r8)) &
               * dsylunit * sed_frc(:,i)
         ENDDO
      ENDDO

   END SUBROUTINE calc_sediment_yield
```

**Step 2: 提交**

```bash
git add main/HYDRO/MOD_Grid_RiverLakeSediment.F90
git commit -m "feat: implement sediment yield calculation"
```

---

## Phase 4: 模块集成

### Task 4.1: 实现主计算子程序

**Files:**
- Modify: `main/HYDRO/MOD_Grid_RiverLakeSediment.F90`

**Step 1: 完善 grid_sediment_calc 子程序**

```fortran
   SUBROUTINE grid_sediment_calc(deltime)
   USE MOD_Grid_RiverLakeNetwork, only: numucat, topo_rivwth, topo_rivlen, &
      topo_rivman, topo_area
   USE MOD_Grid_RiverLakeTimeVars, only: wdsrf_ucat, veloc_riv
   USE MOD_Grid_RiverLakeHist, only: a_floodarea
   IMPLICIT NONE

   real(r8), intent(in) :: deltime

   real(r8) :: sed_time_remaining, dt_sed
   real(r8) :: avg_veloc, avg_wdsrf
   real(r8), allocatable :: rivsto(:), rivout(:), fldfrc(:)
   integer  :: i, ised

      IF (.not. DEF_USE_SEDIMENT) RETURN
      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      allocate(rivsto(numucat))
      allocate(rivout(numucat))
      allocate(fldfrc(numucat))

      sed_time_remaining = deltime

      DO WHILE (sed_time_remaining > 0._r8)

         dt_sed = min(sed_time_remaining, DEF_SED_DT_MAX)

         ! Calculate average water flow variables
         IF (sed_acc_time > 0._r8) THEN
            DO i = 1, numucat
               avg_veloc = sed_acc_veloc(i) / sed_acc_time
               avg_wdsrf = sed_acc_wdsrf(i) / sed_acc_time

               ! Calculate shear velocity
               shearvel(i) = calc_shear_velocity(avg_veloc, avg_wdsrf, topo_rivman(i))

               ! Calculate critical shear velocity (Egiazoroff)
               CALL calc_critical_shear_egiazoroff(i, shearvel(i), critshearvel(:,i))

               ! Calculate suspension velocity
               CALL calc_suspend_velocity(critshearvel(:,i), shearvel(i), susvel(:,i))

               ! Estimate river storage and outflow
               rivsto(i) = avg_wdsrf * topo_rivwth(i) * topo_rivlen(i)
               rivout(i) = avg_veloc * avg_wdsrf * topo_rivwth(i)

               ! Estimate flooded fraction
               IF (a_floodarea(i) > 0._r8 .and. topo_area(i) > 0._r8) THEN
                  fldfrc(i) = a_floodarea(i) / topo_area(i) / sed_acc_time
               ELSE
                  fldfrc(i) = 0._r8
               ENDIF
            ENDDO
         ENDIF

         ! Calculate sediment yield from precipitation
         CALL calc_sediment_yield(fldfrc, topo_area)

         ! Calculate advection
         CALL calc_sediment_advection(dt_sed, rivout, rivsto)

         ! Calculate suspension-deposition exchange
         CALL calc_sediment_exchange(dt_sed, rivsto, topo_rivwth, topo_rivlen)

         ! Redistribute bed layers
         CALL calc_layer_redistribution(topo_rivwth, topo_rivlen)

         ! Accumulate for output
         CALL accumulate_sediment_output(dt_sed)

         sed_time_remaining = sed_time_remaining - dt_sed
      ENDDO

      ! Reset accumulation variables
      sed_acc_time = 0._r8
      sed_acc_veloc(:) = 0._r8
      sed_acc_wdsrf(:) = 0._r8
      sed_precip(:) = 0._r8

      deallocate(rivsto, rivout, fldfrc)

   END SUBROUTINE grid_sediment_calc
```

**Step 2: 添加累积输出子程序**

```fortran
   SUBROUTINE accumulate_sediment_output(dt)
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE
   real(r8), intent(in) :: dt
   integer :: i

      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      DO i = 1, numucat
         a_sedcon(:,i) = a_sedcon(:,i) + sedcon(:,i) * dt
         a_sedout(:,i) = a_sedout(:,i) + sedout(:,i) * dt
         a_bedout(:,i) = a_bedout(:,i) + bedout(:,i) * dt
         a_sedinp(:,i) = a_sedinp(:,i) + sedinp(:,i) * dt
         a_netflw(:,i) = a_netflw(:,i) + netflw(:,i) * dt
         a_layer(:,i)  = a_layer(:,i)  + layer(:,i)  * dt
         a_shearvel(i) = a_shearvel(i) + shearvel(i) * dt
      ENDDO

   END SUBROUTINE accumulate_sediment_output
```

**Step 3: 提交**

```bash
git add main/HYDRO/MOD_Grid_RiverLakeSediment.F90
git commit -m "feat: implement main sediment calculation routine"
```

---

### Task 4.2: 实现诊断累积子程序

**Files:**
- Modify: `main/HYDRO/MOD_Grid_RiverLakeSediment.F90`

**Step 1: 完善 sediment_diag_accumulate 子程序**

```fortran
   SUBROUTINE sediment_diag_accumulate(dt, veloc, wdsrf)
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE
   real(r8), intent(in) :: dt
   real(r8), intent(in) :: veloc(:)
   real(r8), intent(in) :: wdsrf(:)
   integer :: i

      IF (.not. DEF_USE_SEDIMENT) RETURN
      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      sed_acc_time = sed_acc_time + dt

      DO i = 1, numucat
         sed_acc_veloc(i) = sed_acc_veloc(i) + veloc(i) * dt
         sed_acc_wdsrf(i) = sed_acc_wdsrf(i) + wdsrf(i) * dt
      ENDDO

   END SUBROUTINE sediment_diag_accumulate
```

**Step 2: 完善 sediment_forcing_put 子程序**

```fortran
   SUBROUTINE sediment_forcing_put(precip)
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE
   real(r8), intent(in) :: precip(:)

      IF (.not. DEF_USE_SEDIMENT) RETURN
      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      sed_precip(:) = precip(:)

   END SUBROUTINE sediment_forcing_put
```

**Step 3: 提交**

```bash
git add main/HYDRO/MOD_Grid_RiverLakeSediment.F90
git commit -m "feat: implement sediment diagnostic accumulation"
```

---

### Task 4.3: 在 Flow 模块中添加调用

**Files:**
- Modify: `main/HYDRO/MOD_Grid_RiverLakeFlow.F90`

**Step 1: 添加模块引用**

在文件开头 USE 语句区域添加：

```fortran
#ifdef GridRiverLakeSediment
   USE MOD_Grid_RiverLakeSediment, only: grid_sediment_calc, sediment_diag_accumulate
#endif
```

**Step 2: 在水流子循环中添加诊断累积**

在 `DO WHILE (any(dt_res > 0))` 循环内，时间步完成后添加：

```fortran
            ! ... existing code to update dt_res ...
            dt_res = dt_res - dt_all

#ifdef GridRiverLakeSediment
            ! Accumulate water variables for sediment calculation
            IF (DEF_USE_SEDIMENT) THEN
               DO i = 1, numucat
                  IF (ucatfilter(i)) THEN
                     CALL sediment_diag_accumulate(dt_all(irivsys(i)), &
                        veloc_riv(i:i), wdsrf_ucat(i:i))
                  ENDIF
               ENDDO
            ENDIF
#endif
```

**Step 3: 在水流计算完成后调用泥沙计算**

在 `acctime_rnof = 0.` 之前添加：

```fortran
#ifdef GridRiverLakeSediment
      IF (DEF_USE_SEDIMENT) THEN
         CALL grid_sediment_calc(acctime_rnof_max)
      ENDIF
#endif
```

**Step 4: 提交**

```bash
git add main/HYDRO/MOD_Grid_RiverLakeFlow.F90
git commit -m "feat: integrate sediment module into flow calculation"
```

---

## Phase 5: 历史输出与 Restart

### Task 5.1: 添加泥沙历史输出

**Files:**
- Modify: `main/HYDRO/MOD_Grid_RiverLakeHist.F90`

**Step 1: 添加模块引用和变量声明**

**Step 2: 在 hist_grid_riverlake_init 中分配泥沙输出变量**

**Step 3: 在 hist_grid_riverlake_out 中写出泥沙变量**

**Step 4: 在 flush_acc_fluxes_riverlake 中重置泥沙累积变量**

**Step 5: 在 hist_grid_riverlake_final 中释放泥沙变量**

**Step 6: 提交**

```bash
git add main/HYDRO/MOD_Grid_RiverLakeHist.F90
git commit -m "feat: add sediment history output"
```

---

### Task 5.2: 添加泥沙 Restart 支持

**Files:**
- Modify: `main/HYDRO/MOD_Grid_RiverLakeTimeVars.F90`

**Step 1: 在 READ_GridRiverLakeTimeVars 中读取泥沙状态**

**Step 2: 在 WRITE_GridRiverLakeTimeVars 中写出泥沙状态**

**Step 3: 提交**

```bash
git add main/HYDRO/MOD_Grid_RiverLakeTimeVars.F90
git commit -m "feat: add sediment restart support"
```

---

## 验证与测试

### Task 6.1: 编译测试

```bash
make clean
make -j4
```

确保无编译错误。

### Task 6.2: 单点测试

使用小区域配置文件运行，验证：
1. 泥沙模块正确初始化
2. 输出变量有合理值
3. 质量守恒

### Task 6.3: 与 CaMa-Flood 对比

使用相同输入数据，对比：
1. 悬沙浓度分布
2. 推移质通量
3. 沉积层变化

---

## 文件清单

**新建文件：**
- `main/HYDRO/MOD_Grid_RiverLakeSediment.F90`

**修改文件：**
- `include/define.h`
- `main/MOD_Namelist.F90`
- `main/HYDRO/MOD_Grid_RiverLakeFlow.F90`
- `main/HYDRO/MOD_Grid_RiverLakeHist.F90`
- `main/HYDRO/MOD_Grid_RiverLakeTimeVars.F90`
- `Makefile`
