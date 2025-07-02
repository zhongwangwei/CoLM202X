#include <define.h>

MODULE MOD_Tracer_Forcing
   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_Grid
   USE MOD_SpatialMapping
   USE MOD_TimeManager
   USE MOD_SPMD_Task
   USE MOD_MonthlyinSituCO2MaunaLoa
   USE MOD_Vars_Global, only: pi
   USE MOD_OrbCoszen
   USE MOD_UserDefFun
   USE MOD_Tracer_Namelist_Defs, only: DEF_Tracers, nl_tracer_forcing_type, DEF_Tracer_Forcings_NL, MAX_TRACER_FORCING_VARS
   USE MOD_Tracer_Vars_1DForcing, only: forc_prc_O18, forc_prl_O18, forc_q_O18, &
                                         Tracer_allocate_1D_Forcing_O18, Tracer_allocate_1D_Forcing_H2

   IMPLICIT NONE
   SAVE
   type :: tracer_dataset_type
      character(len=256) :: tracer_name              ! 数据集名称
      character(len=256) :: dataset_name             ! 数据集名称
      character(len=16)  :: tracer_type              ! 数据集类型
      character(len=256) :: tracer_dir               ! 数据源目录
      integer            :: startyr                  ! 开始年份
      integer            :: startmo                  ! 开始月份
      integer            :: endyr                    ! 结束年份
      integer            :: endmo                    ! 结束月份
      logical            :: leapyear                          
      logical            :: data2d                             
      logical            :: hightdim                           
      logical            :: dim2d                             
      character(len=256) :: latname   
      character(len=256) :: lonname   
      character(len=256) :: groupby   

      integer            :: NVAR                             
      logical            :: has_missing_value                
      character(len=256) :: missing_value_name               
      logical            :: regional                          
      real(r8)           :: regbnd(4)                          

      character(len=64), allocatable :: vnames(:)      ! 变量名称列表
      character(len=64), allocatable :: tintalgo(:)    ! 时间插值算法
      character(len=64), allocatable :: timelog(:)     ! 时间记录方式
      character(len=256), allocatable :: fprefix(:)    ! 文件前缀列表
      integer,           allocatable :: dtime(:)       ! 时间步长
      integer,           allocatable :: offset(:)      ! 时间偏移
      type(grid_type)    :: tracegrid              ! 网格结构
      type(spatial_mapping_type) :: mg2p      ! 空间映射
      ! 数据存储
      type(block_data_real8_2d), allocatable :: forcn(:)      ! 当前强迫数据
      type(block_data_real8_2d), allocatable :: forcn_LB(:)   ! 下边界数据
      type(block_data_real8_2d), allocatable :: forcn_UB(:)   ! 上边界数据
      
      ! 时间戳
      type(timestamp), allocatable :: tstamp_LB(:)
      type(timestamp), allocatable :: tstamp_UB(:)
      
      ! 掩码
      logical, allocatable :: forcmask_pch(:)

      type(block_data_real8_2d) :: avgcos   ! time-average of cos(zenith)
      type(block_data_real8_2d) :: metdata  ! forcing data
   end type tracer_dataset_type


   type :: multi_tracer_forcing_type
      integer :: ndataset                                    ! 数据集数量
      type(tracer_dataset_type), allocatable :: datasets(:) ! 数据集数组
   
      ! 变量到数据集的映射
      integer, allocatable :: var_to_dataset(:)  ! 每个变量对应的数据集ID
      integer, allocatable :: var_to_index(:)    ! 每个变量在数据集中的索引
   end type multi_tracer_forcing_type


   type(multi_tracer_forcing_type), PUBLIC :: mtf
   type(grid_type)    :: gforc              ! 网格结构

   logical, allocatable :: forcmask_pch (:)  ! 全局forcing掩码，供历史输出使用

   ! local variables
   integer  :: deltim_int                ! model time step length
   ! real(r8) :: deltim_real             ! model time step length
   !real(r8), allocatable :: dummy_data(:)    ! dummy data for MPI communication in range check

   ! Variables for array allocation

   !  for SinglePoint
   PUBLIC :: tracer_forcing_init
   PUBLIC :: read_tracer_forcing
   PUBLIC :: tracer_forcing_final
   PUBLIC :: tracer_forcing_reset
   PUBLIC :: initialize_multi_tracer_forcing
   PUBLIC :: cleanup_multi_tracer_forcing
   PUBLIC :: tracer_setstampLB
   PUBLIC :: tracer_setstampUB
   PUBLIC :: tracer_metreadLBUB
   PUBLIC :: tracer_metpreprocess
   PUBLIC :: get_tracer_forcing_data
   PUBLIC :: forcmask_pch

CONTAINS

   SUBROUTINE tracer_forcing_init (deltatime, ststamp, lc_year, etstamp, lulcc_call)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Block
   USE MOD_DataType
   USE MOD_Mesh
   USE MOD_LandElm
   USE MOD_LandPatch
#ifdef CROP
   USE MOD_LandCrop
#endif
   USE MOD_NetCDFSerial
   USE MOD_NetCDFVector
   USE MOD_NetCDFBlock
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_1DForcing
   IMPLICIT NONE

   real(r8),         intent(in) :: deltatime  ! model time step
   type(timestamp),  intent(in) :: ststamp
   integer,          intent(in) :: lc_year    ! which year of land cover data used
   type(timestamp),  intent(in), optional :: etstamp
   logical,          intent(in), optional :: lulcc_call ! whether it is a lulcc CALL

   ! Local variables
   integer            :: idate(3)
   type(timestamp)    :: tstamp
   character(len=256) :: filename, lndname, cyear
   integer            :: ivar, year, month, day, time_i
   real(r8)           :: missing_value
   integer            :: ielm, istt, iend
   integer            :: i, j, k, di

   integer :: iblkme, xblk, yblk, xloc, yloc

   integer :: num_tracers_with_forcing, max_nvar
   integer :: dataset_idx, var_count, total_vars

   ! Check allocation status of tracer arrays

   num_tracers_with_forcing = SIZE(DEF_Tracer_Forcings_NL, 1)
   
   ! 计算最大变量数
   max_nvar = 0
   DO k = 1, num_tracers_with_forcing
      IF (DEF_Tracer_Forcings_NL(k)%NVAR > max_nvar) THEN
         max_nvar = DEF_Tracer_Forcings_NL(k)%NVAR
      ENDIF
   ENDDO

   IF (p_is_master) THEN
      WRITE(*,*) "=== 初始化 Tracer Forcing ==="
      WRITE(*,*) "检测到的tracer数量:", num_tracers_with_forcing
      WRITE(*,*) "最大变量数:", max_nvar
   ENDIF

   ! 保持原有的1D forcing分配不变
   DO k = 1, num_tracers_with_forcing
         SELECT CASE (trim(DEF_Tracer_Forcings_NL(k)%tracer_name))
         CASE ('O18')
            CALL Tracer_allocate_1D_Forcing_O18
         CASE ('H2')
            CALL Tracer_allocate_1D_Forcing_H2
         END SELECT
   ENDDO

   ! get value of deltim (moved outside loop)
   deltim_int  = int(deltatime)
   idate = (/ststamp%year, ststamp%day, ststamp%sec/)

   CALL adj2begin (idate)
   
   ! 初始化multi_tracer_forcing结构
   CALL initialize_multi_tracer_forcing(idate)

   ! 为每个数据集建立空间映射
   DO di = 1, mtf%ndataset
      IF (trim(DEF_Forcing_Interp_Method) == 'arealweight') THEN
         IF (present(lulcc_call)) CALL mtf%datasets(di)%mg2p%forc_free_mem
         CALL mtf%datasets(di)%mg2p%build_arealweighted (mtf%datasets(di)%tracegrid, landpatch)
      ELSEIF (trim(DEF_Forcing_Interp_Method) == 'bilinear') THEN
         IF (present(lulcc_call)) CALL mtf%datasets(di)%mg2p%forc_free_mem
         CALL mtf%datasets(di)%mg2p%build_bilinear (mtf%datasets(di)%tracegrid, landpatch)
      ENDIF

      ! 处理缺失值
      IF (mtf%datasets(di)%has_missing_value) THEN
         tstamp = idate
         CALL tracer_setstampLB(tstamp, di, 1, year, month, day, time_i)
         filename = trim(mtf%datasets(di)%tracer_dir)//&
                   trim(tracerfilename(year, month, day, di, 1))

         IF (p_is_master) THEN
            CALL ncio_get_attr (filename, mtf%datasets(di)%vnames(1), &
                               trim(mtf%datasets(di)%missing_value_name), missing_value)
         ENDIF
#ifdef USEMPI
         CALL mpi_bcast (missing_value, 1, MPI_REAL8, p_address_master, p_comm_glb, p_err)
#endif

         CALL ncio_read_block_time (filename, mtf%datasets(di)%vnames(1), &
                                   mtf%datasets(di)%tracegrid, time_i, mtf%datasets(di)%metdata)

         CALL mtf%datasets(di)%mg2p%set_missing_value (mtf%datasets(di)%metdata, &
                                                      missing_value, mtf%datasets(di)%forcmask_pch)
      ENDIF
   ENDDO

   ! 初始化全局forcmask_pch变量，供历史输出使用
   IF (p_is_worker) THEN
      IF (numpatch > 0) THEN
         allocate(forcmask_pch(numpatch))
         forcmask_pch(:) = .true.  ! 默认所有patch都有效
         
         ! 如果有数据集有missing value处理，则使用第一个数据集的掩码
         ! （或者可以实现所有数据集掩码的逻辑AND）
         DO di = 1, mtf%ndataset
            IF (mtf%datasets(di)%has_missing_value .and. allocated(mtf%datasets(di)%forcmask_pch)) THEN
               forcmask_pch(:) = forcmask_pch(:) .and. mtf%datasets(di)%forcmask_pch(:)
               EXIT  ! 使用第一个有缺失值处理的数据集的掩码
            ENDIF
         ENDDO
      ENDIF
   ENDIF

   IF (p_is_master) WRITE(*,*) "    --- Tracer Forcing初始化完成 ---"
   END SUBROUTINE tracer_forcing_init

   ! ---- forcing finalize ----
   SUBROUTINE tracer_forcing_final ()

   USE MOD_LandPatch, only: numpatch
   IMPLICIT NONE

   integer :: i, j, di

   ! 清理每个数据集的空间映射
   DO di = 1, mtf%ndataset
      IF (allocated(mtf%datasets)) THEN
         ! 清理空间映射内存
         CALL mtf%datasets(di)%mg2p%forc_free_mem
      ENDIF
   ENDDO

   ! 清理全局forcmask_pch
   IF (allocated(forcmask_pch)) deallocate(forcmask_pch)

   ! 清理multi_tracer_forcing结构
   CALL cleanup_multi_tracer_forcing()

   END SUBROUTINE tracer_forcing_final

   ! ------------
   SUBROUTINE tracer_forcing_reset ()

   IMPLICIT NONE

   integer :: di, ivar

   ! 重置每个数据集的时间戳
   IF (allocated(mtf%datasets)) THEN
      DO di = 1, mtf%ndataset
         IF (allocated(mtf%datasets(di)%tstamp_LB)) THEN
            DO ivar = 1, mtf%datasets(di)%NVAR
               mtf%datasets(di)%tstamp_LB(ivar) = timestamp(-1, -1, -1)
               mtf%datasets(di)%tstamp_UB(ivar) = timestamp(-1, -1, -1)
            ENDDO
         ENDIF
      ENDDO
   ENDIF

   END SUBROUTINE tracer_forcing_reset


!-----------------------------------------------------------------------
   SUBROUTINE read_tracer_forcing (idate)
   USE MOD_OrbCosazi
   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_Const_Physical, only: rgas, grav
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_TimeVariables, only: alb
   USE MOD_Vars_1DForcing
   USE MOD_Block
   USE MOD_SPMD_Task
   USE MOD_DataType
   USE MOD_Mesh
   USE MOD_LandPatch
   USE MOD_RangeCheck
   USE MOD_NetCDFVector

   IMPLICIT NONE

   integer, intent(in) :: idate(3)

   ! local variables:
   integer  :: di, ivar, istt, iend, id(3)
   integer  :: iblkme, ib, jb, i, j, ilon, ilat, np, ipart, ne
   real(r8) :: calday                             ! Julian cal day (1.xx to 365.xx)
   real(r8) :: sunang, cloud, difrat, vnrat
   real(r8) :: a, hsolar, ratio_rvrf
   integer  :: ii
   character(10) :: cyear = "2005"
   character(256):: lndname

   type(timestamp) :: mtstamp
   integer  :: dtLB, dtUB
   real(r8) :: cosz, coszen(numpatch), cosa, cosazi(numpatch), balb
   integer  :: year, month, mday
   logical  :: has_u,has_v
   real solar, frl, prcp, tm, us, vs, pres, qm
   real(r8) :: pco2m
   integer target_server, ierr

   IF (p_is_io) THEN
      !------------------------------------------------------------
      ! READ in THE TRACER FORCING DATA
      ! 对每个数据集进行处理
      DO di = 1, mtf%ndataset
         
         ! read lower and upper boundary forcing data
         CALL tracer_metreadLBUB(idate, di)
         
         ! set model time stamp
         id(:) = idate(:)
         mtstamp = id
         
         ! loop for variables in this dataset
         DO ivar = 1, mtf%datasets(di)%NVAR
            IF (trim(mtf%datasets(di)%vnames(ivar)) == 'NULL') CYCLE     ! no data, CYCLE
            IF (trim(mtf%datasets(di)%tintalgo(ivar)) == 'NULL') CYCLE

            ! to make sure the forcing data calculated is in the range of time
            ! interval [LB, UB]
            IF ( (mtstamp < mtf%datasets(di)%tstamp_LB(ivar)) .or. &
                 (mtf%datasets(di)%tstamp_UB(ivar) < mtstamp) ) THEN
               write(6, *) "the tracer data required is out of range! STOP!"; CALL CoLM_stop()
            ENDIF

            ! calculate distance to lower/upper boundary
            dtLB = mtstamp - mtf%datasets(di)%tstamp_LB(ivar)
            dtUB = mtf%datasets(di)%tstamp_UB(ivar) - mtstamp

            ! linear method, for T, Pres, Q, W, LW
            IF (mtf%datasets(di)%tintalgo(ivar) == 'linear') THEN
               IF ( (dtLB+dtUB) > 0 ) THEN
                  CALL block_data_linear_interp ( &
                     mtf%datasets(di)%forcn_LB(ivar), real(dtUB,r8)/real(dtLB+dtUB,r8), &
                     mtf%datasets(di)%forcn_UB(ivar), real(dtLB,r8)/real(dtLB+dtUB,r8), &
                     mtf%datasets(di)%forcn(ivar))
               ELSE
                  CALL block_data_copy (mtf%datasets(di)%forcn_LB(ivar), mtf%datasets(di)%forcn(ivar))
               ENDIF
            ENDIF

            ! for precipitation, two algorithms available
            ! nearest method, for precipitation
            IF (mtf%datasets(di)%tintalgo(ivar) == 'nearest') THEN
               IF (dtLB <= dtUB) THEN
                  CALL block_data_copy (mtf%datasets(di)%forcn_LB(ivar), mtf%datasets(di)%forcn(ivar))
               ELSE
                  CALL block_data_copy (mtf%datasets(di)%forcn_UB(ivar), mtf%datasets(di)%forcn(ivar))
               ENDIF
            ENDIF

            ! set all the same value, for precipitation
            IF (mtf%datasets(di)%tintalgo(ivar) == 'uniform') THEN
               IF (trim(mtf%datasets(di)%timelog(ivar)) == 'forward') THEN
                  CALL block_data_copy (mtf%datasets(di)%forcn_LB(ivar), mtf%datasets(di)%forcn(ivar))
               ELSE
                  CALL block_data_copy (mtf%datasets(di)%forcn_UB(ivar), mtf%datasets(di)%forcn(ivar))
               ENDIF
            ENDIF

            ! coszen method, for SW
            IF (mtf%datasets(di)%tintalgo(ivar) == 'coszen') THEN
               DO iblkme = 1, gblock%nblkme
                  ib = gblock%xblkme(iblkme)
                  jb = gblock%yblkme(iblkme)

                  DO j = 1, mtf%datasets(di)%tracegrid%ycnt(jb)
                     DO i = 1, mtf%datasets(di)%tracegrid%xcnt(ib)

                        ilat = mtf%datasets(di)%tracegrid%ydsp(jb) + j
                        ilon = mtf%datasets(di)%tracegrid%xdsp(ib) + i
                        IF (ilon > mtf%datasets(di)%tracegrid%nlon) &
                           ilon = ilon - mtf%datasets(di)%tracegrid%nlon

                        calday = calendarday(mtstamp)
                        cosz = orb_coszen(calday, mtf%datasets(di)%tracegrid%rlon(ilon), &
                                         mtf%datasets(di)%tracegrid%rlat(ilat))
                        cosz = max(0.001, cosz)
                        
                        IF (trim(mtf%datasets(di)%timelog(ivar)) == 'forward') THEN
                           mtf%datasets(di)%forcn(ivar)%blk(ib,jb)%val(i,j) = &
                              cosz / mtf%datasets(di)%avgcos%blk(ib,jb)%val(i,j) * &
                              mtf%datasets(di)%forcn_LB(ivar)%blk(ib,jb)%val(i,j)
                        ELSE
                           mtf%datasets(di)%forcn(ivar)%blk(ib,jb)%val(i,j) = &
                              cosz / mtf%datasets(di)%avgcos%blk(ib,jb)%val(i,j) * &
                              mtf%datasets(di)%forcn_UB(ivar)%blk(ib,jb)%val(i,j)
                        ENDIF

                     ENDDO
                  ENDDO
               ENDDO
            ENDIF

         ENDDO

         ! tracer specific preprocessing (if needed)
         CALL tracer_metpreprocess (mtf%datasets(di)%tracegrid, mtf%datasets(di)%forcn, di)

      ENDDO
   ENDIF

   ! 将数据从网格映射到patch - 参考MOD_Forcing.F90的简洁风格
   ! Mapping the 2d tracer fields [lon_points]x[lat_points]
   !     -> the 1d vector of subgrid points [numpatch]
   
   DO di = 1, mtf%ndataset
      DO ivar = 1, mtf%datasets(di)%NVAR
         IF (trim(mtf%datasets(di)%vnames(ivar)) == 'NULL') CYCLE     
         ! 直接映射，类似MOD_Forcing.F90的方式
         SELECT CASE (trim(mtf%datasets(di)%tracer_name))
         CASE ('O18')
            SELECT CASE (trim(mtf%datasets(di)%vnames(ivar)))
            CASE ('prc_O18', 'prc_o18', 'PRC_O18', 'prate1sfc')
               CALL mtf%datasets(di)%mg2p%grid2pset (mtf%datasets(di)%forcn(ivar), forc_prc_O18)
               CALL mtf%datasets(di)%mg2p%grid2pset (mtf%datasets(di)%forcn(ivar), forc_prl_O18)
            CASE ('q_O18', 'q_o18', 'Q_O18', 'spfh12m')
               CALL mtf%datasets(di)%mg2p%grid2pset (mtf%datasets(di)%forcn(ivar), forc_q_O18)
            END SELECT
         CASE ('H2')
            ! H2的映射可以在这里添加
             SELECT CASE (trim(mtf%datasets(di)%vnames(ivar)))
            CASE ('prc_H2', 'prc_h2', 'PRC_H2', 'prate2sfc')
               CALL mtf%datasets(di)%mg2p%grid2pset (mtf%datasets(di)%forcn(ivar), forc_prc_H2)
            CASE ('q_H2', 'q_h2', 'Q_H2', 'spfh2m')
               CALL mtf%datasets(di)%mg2p%grid2pset (mtf%datasets(di)%forcn(ivar), forc_q_H2)
            END SELECT
         END SELECT
      ENDDO
   ENDDO
   
   ! 后处理 - 类似MOD_Forcing.F90在worker进程上进行数据处理
   IF (p_is_worker) THEN
      ! O18 precipitation splitting (如果需要从总降水分割为对流和大尺度降水)
      IF (allocated(forc_prc_O18) .AND. allocated(forc_prl_O18)) THEN
         ! 可选：如果只有一个总降水变量，可以在这里进行分割
         ! 例如，如果prc_O18包含总降水，可以这样分割：
         ! forc_prl_O18(:) = forc_prc_O18(:)  ! 先复制
          forc_prc_O18 = forc_prc_O18 * (1.0_r8/3.0_r8)  ! convective (1/3)
          forc_prl_O18 = forc_prl_O18 * (2.0_r8/3.0_r8)  ! large scale (2/3)
      ENDIF
   ENDIF

#ifdef RangeCheck
#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif
   IF (p_is_master) write(*,'(/, A25)') 'Checking tracer forcing ...'
   ! Check O18 tracer forcing data
      CALL check_vector_data ('Tracer O18 prc [mm/s]  ', forc_prc_O18)
      CALL check_vector_data ('Tracer O18 prl [mm/s]  ', forc_prl_O18)
      CALL check_vector_data ('Tracer O18 q   [kg/kg] ', forc_q_O18)
   
   ! Check O18 specific humidity  
  !    IF (p_is_worker .AND. allocated(forc_q_O18)) THEN
  !       CALL check_vector_data ('Tracer O18 q1   [kg/kg] ', forc_q_O18)
  !     ELSEIF (p_is_master) THEN
  !        ! Master process participates in MPI communication with dummy data
  !        IF (.NOT. allocated(dummy_data)) THEN
  !           allocate(dummy_data(1))
  !           dummy_data(1) = 0.0_r8
  !        ENDIF
  !        CALL check_vector_data ('dummy', dummy_data)
  !    ENDIF
      
   ! Check H2 tracer forcing data (when available)
   ! IF (allocated(forc_prc_H2)) THEN
   !    CALL check_vector_data ('Tracer H2 prc  [mm/s]  ', forc_prc_H2)
   ! ENDIF
   ! IF (allocated(forc_prl_H2)) THEN
   !    CALL check_vector_data ('Tracer H2 prl  [mm/s]  ', forc_prl_H2)
   ! ENDIF
   ! IF (allocated(forc_q_H2)) THEN
   !    CALL check_vector_data ('Tracer H2 q    [kg/kg] ', forc_q_H2)
   ! ENDIF

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif
#endif

   END SUBROUTINE read_tracer_forcing





   
   


   FUNCTION tracerfilename(year, month, day, tracer_idx,var_i) RESULT(metfilename)

      USE MOD_Namelist
      USE MOD_SPMD_Task
      IMPLICIT NONE
   
      integer, intent(in) :: year
      integer, intent(in) :: month
      integer, intent(in) :: day
      integer, intent(in) :: tracer_idx
      integer, intent(in) :: var_i
      character(len=256)  :: metfilename
      character(len=256)  :: yearstr
      character(len=256)  :: monthstr
      character(len=256)  :: dataset_name_lower
   
         ! IF (p_is_master) THEN
         !    WRITE(*,*) "      === DEBUG: tracerfilename called ==="
         !    WRITE(*,*) "        Input: year=", year, " month=", month, " day=", day, " var_i=", var_i
         ! ENDIF
   
         write(yearstr, '(I4.4)') year
         write(monthstr, '(I2.2)') month
         
         dataset_name_lower = adjustl(trim(mtf%datasets(tracer_idx)%dataset_name))
         CALL to_lower(dataset_name_lower)
   
         ! IF (p_is_master) THEN
         !    WRITE(*,*) "        Formatted strings: yearstr='", TRIM(yearstr), "' monthstr='", TRIM(monthstr), "'"
         !    WRITE(*,*) "        Dataset name (original): '", TRIM(mtf%datasets(tracer_idx)%dataset_name), "'"
         !    WRITE(*,*) "        Dataset name (lower): '", dataset_name_lower, "'"
         ! ENDIF
   
         select CASE (dataset_name_lower)
         
         CASE ('isogsm')
            !DESCRIPTION
            !===========
               !--- Isotopes-incorporated Global Spectral Model (IsoGSM)
      
            !data source:
            !-------------------
               !---https://isotope.iis.u-tokyo.ac.jp/about-our-lab?lang=en
      
            !References:
            !-------------------
               !---Bong, H., Cauquoin, A., Okazaki, A., Chang, E.-C., Werner, M., Wei, Z., et al. (2024). 
               !   Process-based intercomparison of water isotope-enabled models and reanalysis nudging effects. 
               !   Journal of Geophysical Research: Atmospheres, 129, e2023JD038719. 
               !   https://doi.org/10.1029/2023JD038719
      
            !REVISION HISTORY
            !----------------
               !---2025.03.23   Zhongwang Wei @ SYSU: add the isotope forcing data
      
               ! Check if var_i is within bounds for fprefix array


               ! Construct filename based on the specific prefix for the variable, year, and .nc suffix
               metfilename = '/'//trim(mtf%datasets(tracer_idx)%fprefix(var_i))//'_'//trim(yearstr)//'.nc'

   
         
         CASE ('POINT')
            metfilename = '/'//trim(mtf%datasets(tracer_idx)%fprefix(1))
            ! IF (p_is_master) THEN
            !    WRITE(*,*) "        CASE POINT selected"
            !    WRITE(*,*) "        fprefix_tracer_forcing(1) = '", TRIM(mtf%datasets(tracer_idx)%fprefix(1)), "'"
            !    WRITE(*,*) "        Generated filename: '", TRIM(metfilename), "'"
            ! ENDIF
            
         CASE DEFAULT
            ! IF (p_is_master) THEN
            !    WRITE(*,*) "        WARNING: Unknown dataset name '", TRIM(mtf%datasets(tracer_idx)%dataset_name), "'"
            !    WRITE(*,*) "        Using default POINT format"
            ! ENDIF
            metfilename = '/'//trim(mtf%datasets(tracer_idx)%fprefix(1))
         END select
         
         ! IF (p_is_master) THEN
         !    WRITE(*,*) "      === DEBUG: tracerfilename returning '", TRIM(metfilename), "' ==="
         ! ENDIF
         
         ! IF (DEF_USE_CBL_HEIGHT) THEN
         !    select CASE (var_i)
         !    CASE (9)
         !       metfilename = '/'//trim(fprefix_tracer_forcing(9))//'_'//trim(yearstr)//'_'//trim(monthstr)//&
         !          '_boundary_layer_height.nc4'
         !    END select
         ! ENDIF
      END FUNCTION tracerfilename

   ! 初始化multi_tracer_forcing结构
   
   
   SUBROUTINE initialize_multi_tracer_forcing(idate)
   
   USE MOD_SPMD_Task
   IMPLICIT NONE
   
   integer, intent(in) :: idate(3)
   
   integer :: i, j, k, dataset_idx, var_idx
   integer :: num_tracers_with_forcing, total_vars
   character(len=256) :: current_dataset
   logical :: dataset_exists
   
   ! 获取tracer数量
   num_tracers_with_forcing = SIZE(DEF_Tracer_Forcings_NL, 1)
   
   IF (p_is_master) THEN
      WRITE(*,*) "开始初始化multi_tracer_forcing结构..."
   ENDIF
   
   ! 首先统计有多少个不同的数据集
   mtf%ndataset = num_tracers_with_forcing
   
   allocate(mtf%datasets(mtf%ndataset))


   DO i = 1, mtf%ndataset
      CALL initialize_single_dataset(i,idate)
   ENDDO

   
   ! 计算总变量数并建立映射
   total_vars = 0
   DO i = 1, num_tracers_with_forcing
      IF (trim(DEF_Tracer_Forcings_NL(i)%tracer_name) /= '') THEN
         total_vars = total_vars + DEF_Tracer_Forcings_NL(i)%NVAR
      ENDIF
   ENDDO
   
   ! 分配映射数组
   IF (total_vars > 0) THEN
      allocate(mtf%var_to_dataset(total_vars))
      allocate(mtf%var_to_index(total_vars))
      
      ! 建立变量到数据集的映射
      var_idx = 0
      DO i = 1, num_tracers_with_forcing
         IF (trim(DEF_Tracer_Forcings_NL(i)%tracer_name) /= '') THEN
            current_dataset = trim(mtf%datasets(i)%dataset_name)
            
            ! 找到对应的数据集索引
            dataset_idx = 0
            DO j = 1, mtf%ndataset
               IF (trim(mtf%datasets(j)%dataset_name) == trim(current_dataset)) THEN
                  dataset_idx = j
                  EXIT
               ENDIF
            ENDDO
            
            ! 为这个tracer的所有变量建立映射
            DO k = 1, DEF_Tracer_Forcings_NL(i)%NVAR
               var_idx = var_idx + 1
               mtf%var_to_dataset(var_idx) = dataset_idx
               mtf%var_to_index(var_idx) = k
            ENDDO
         ENDIF
      ENDDO
   ENDIF
   
   IF (p_is_master) THEN
      WRITE(*,*) "multi_tracer_forcing初始化完成!"
      WRITE(*,*) "总变量数:", total_vars
   ENDIF
   
   END SUBROUTINE initialize_multi_tracer_forcing

   ! 初始化单个数据集
   SUBROUTINE initialize_single_dataset(di,idate)
   
   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial
   USE MOD_LandPatch, only: numpatch
   IMPLICIT NONE
   
   integer, intent(in) :: di    ! 在DEF_Tracer_Forcings_NL中的索引
   integer, intent(in) :: idate(3)

   integer :: i
   type(timestamp) :: tstamp
   character(len=256) :: filename
   integer :: year, month, day, time_i
   real(r8), allocatable :: latxy (:,:)    ! latitude values in 2d
   real(r8), allocatable :: lonxy (:,:)    ! longitude values in 2d
   real(r8), allocatable :: lon_in(:)
   real(r8), allocatable :: lat_in(:)
   
   ! 复制基本信息从namelist到dataset结构
   mtf%datasets(di)%tracer_name = DEF_Tracer_Forcings_NL(di)%tracer_name
   mtf%datasets(di)%dataset_name = DEF_Tracer_Forcings_NL(di)%dataset_name
   mtf%datasets(di)%tracer_type = DEF_Tracer_Forcings_NL(di)%tracer_type
   mtf%datasets(di)%tracer_dir = DEF_Tracer_Forcings_NL(di)%tracer_dir
   mtf%datasets(di)%startyr = DEF_Tracer_Forcings_NL(di)%startyr
   mtf%datasets(di)%startmo = DEF_Tracer_Forcings_NL(di)%startmo
   mtf%datasets(di)%endyr = DEF_Tracer_Forcings_NL(di)%endyr
   mtf%datasets(di)%endmo = DEF_Tracer_Forcings_NL(di)%endmo
   mtf%datasets(di)%leapyear = DEF_Tracer_Forcings_NL(di)%leapyear
   mtf%datasets(di)%data2d = DEF_Tracer_Forcings_NL(di)%data2d
   mtf%datasets(di)%hightdim = DEF_Tracer_Forcings_NL(di)%hightdim
   mtf%datasets(di)%dim2d = DEF_Tracer_Forcings_NL(di)%dim2d
   mtf%datasets(di)%latname = DEF_Tracer_Forcings_NL(di)%latname
   mtf%datasets(di)%lonname = DEF_Tracer_Forcings_NL(di)%lonname
   mtf%datasets(di)%groupby = DEF_Tracer_Forcings_NL(di)%groupby
   mtf%datasets(di)%NVAR = DEF_Tracer_Forcings_NL(di)%NVAR
   mtf%datasets(di)%has_missing_value = DEF_Tracer_Forcings_NL(di)%has_missing_value
   mtf%datasets(di)%missing_value_name = DEF_Tracer_Forcings_NL(di)%missing_value_name
   mtf%datasets(di)%regional = DEF_Tracer_Forcings_NL(di)%regional
   mtf%datasets(di)%regbnd = DEF_Tracer_Forcings_NL(di)%regbnd
   
   ! 分配并复制可变大小数组
   IF (DEF_Tracer_Forcings_NL(di)%NVAR > 0) THEN
      allocate(mtf%datasets(di)%vnames(DEF_Tracer_Forcings_NL(di)%NVAR))
      allocate(mtf%datasets(di)%tintalgo(DEF_Tracer_Forcings_NL(di)%NVAR))
      allocate(mtf%datasets(di)%timelog(DEF_Tracer_Forcings_NL(di)%NVAR))
      allocate(mtf%datasets(di)%fprefix(DEF_Tracer_Forcings_NL(di)%NVAR))
      allocate(mtf%datasets(di)%dtime(DEF_Tracer_Forcings_NL(di)%NVAR))
      allocate(mtf%datasets(di)%offset(DEF_Tracer_Forcings_NL(di)%NVAR))
      
      DO i = 1, DEF_Tracer_Forcings_NL(di)%NVAR
         mtf%datasets(di)%vnames(i) = DEF_Tracer_Forcings_NL(di)%vname(i)
         mtf%datasets(di)%tintalgo(i) = DEF_Tracer_Forcings_NL(di)%tintalgo(i)
         mtf%datasets(di)%timelog(i) = DEF_Tracer_Forcings_NL(di)%timelog(i)
         mtf%datasets(di)%fprefix(i) = DEF_Tracer_Forcings_NL(di)%fprefix(i)
         mtf%datasets(di)%dtime(i) = DEF_Tracer_Forcings_NL(di)%dtime(i)
         mtf%datasets(di)%offset(i) = DEF_Tracer_Forcings_NL(di)%offset(i)
      ENDDO
      
      ! 分配时间戳数组 - 所有进程都需要
      allocate(mtf%datasets(di)%tstamp_LB(DEF_Tracer_Forcings_NL(di)%NVAR))
      allocate(mtf%datasets(di)%tstamp_UB(DEF_Tracer_Forcings_NL(di)%NVAR))
      
      ! 初始化时间戳为无效值
      DO i = 1, DEF_Tracer_Forcings_NL(di)%NVAR
         mtf%datasets(di)%tstamp_LB(i) = timestamp(-1, -1, -1)
         mtf%datasets(di)%tstamp_UB(i) = timestamp(-1, -1, -1)
      ENDDO
   ENDIF

   ! 读取网格信息 - 参考主forcing模块的metread_latlon
   ! 这需要在所有进程上进行，因为空间映射需要网格信息
   IF (trim(mtf%datasets(di)%dataset_name) == 'POINT' .or. &
       trim(mtf%datasets(di)%dataset_name) == 'CPL7' ) THEN
      CALL mtf%datasets(di)%tracegrid%define_by_ndims (360, 180)
   ELSE
      ! 使用传入的idate参数获取网格信息
      tstamp = idate
      
      ! 获取年月日信息来构建文件名
      year = tstamp%year
      CALL julian2monthday(year, tstamp%day, month, day)
      time_i = 1
      
      filename = trim(DEF_Tracer_Forcings_NL(di)%tracer_dir)//&
                 trim(tracerfilename(year, month, day, di, 1))

      IF (mtf%datasets(di)%dim2d) THEN
         CALL ncio_read_bcast_serial (filename, mtf%datasets(di)%latname, latxy)
         CALL ncio_read_bcast_serial (filename, mtf%datasets(di)%lonname, lonxy)

         allocate (lat_in (size(latxy,2)))
         allocate (lon_in (size(lonxy,1)))
         lat_in = latxy(1,:)
         lon_in = lonxy(:,1)

         deallocate (latxy)
         deallocate (lonxy)
      ELSE
         CALL ncio_read_bcast_serial (filename, mtf%datasets(di)%latname, lat_in)
         CALL ncio_read_bcast_serial (filename, mtf%datasets(di)%lonname, lon_in)
      ENDIF

      IF (.not. mtf%datasets(di)%regional) THEN
         CALL mtf%datasets(di)%tracegrid%define_by_center (lat_in, lon_in)
      ELSE
         CALL mtf%datasets(di)%tracegrid%define_by_center (lat_in, lon_in, &
            south = mtf%datasets(di)%regbnd(1), &
            north = mtf%datasets(di)%regbnd(2), &
            west  = mtf%datasets(di)%regbnd(3), &
            east  = mtf%datasets(di)%regbnd(4))
      ENDIF

      deallocate (lat_in)
      deallocate (lon_in)
   ENDIF

   CALL mtf%datasets(di)%tracegrid%set_rlon ()
   CALL mtf%datasets(di)%tracegrid%set_rlat ()

   ! 在所有进程上分配数据存储数组，以支持grid2pset分布式操作
   IF (DEF_Tracer_Forcings_NL(di)%NVAR > 0) THEN
      allocate(mtf%datasets(di)%forcn(DEF_Tracer_Forcings_NL(di)%NVAR))
      
      DO i = 1, DEF_Tracer_Forcings_NL(di)%NVAR
         CALL allocate_block_data (mtf%datasets(di)%tracegrid, &
                                  mtf%datasets(di)%forcn(i))
      ENDDO
      
      ! 只在IO进程上分配边界数据和辅助数组，因为只有IO进程需要读取数据
      IF (p_is_io) THEN
         allocate(mtf%datasets(di)%forcn_LB(DEF_Tracer_Forcings_NL(di)%NVAR))
         allocate(mtf%datasets(di)%forcn_UB(DEF_Tracer_Forcings_NL(di)%NVAR))
         
         DO i = 1, DEF_Tracer_Forcings_NL(di)%NVAR
            CALL allocate_block_data (mtf%datasets(di)%tracegrid, &
                                     mtf%datasets(di)%forcn_LB(i))
            CALL allocate_block_data (mtf%datasets(di)%tracegrid, &
                                     mtf%datasets(di)%forcn_UB(i))
         ENDDO
         
         ! 分配辅助数据数组
         CALL allocate_block_data (mtf%datasets(di)%tracegrid, &
                                  mtf%datasets(di)%metdata)
         CALL allocate_block_data (mtf%datasets(di)%tracegrid, &
                                  mtf%datasets(di)%avgcos)
      ENDIF
   ENDIF
   
   ! 在worker进程上分配掩码数组
   IF (p_is_worker) THEN
      IF (numpatch > 0) THEN
         allocate (mtf%datasets(di)%forcmask_pch(numpatch))
         mtf%datasets(di)%forcmask_pch(:) = .true.
      ENDIF
   ENDIF

   ! 处理缺失值（如果需要的话，可以在这里添加）
    if (mtf%datasets(di)%has_missing_value) then
   !    ! 缺失值处理逻辑
    endif
   
   IF (p_is_master) THEN
      WRITE(*,*) "  数据集初始化完成:", trim(mtf%datasets(di)%dataset_name)
      WRITE(*,*) "    变量数:", mtf%datasets(di)%NVAR
   ENDIF
   
   END SUBROUTINE initialize_single_dataset

   ! 清理multi_tracer_forcing结构
   SUBROUTINE cleanup_multi_tracer_forcing()
   
   USE MOD_SPMD_Task
   IMPLICIT NONE
   
   integer :: i, j
   
   IF (p_is_master) THEN
      WRITE(*,*) "清理multi_tracer_forcing结构..."
   ENDIF
   
   ! 清理每个数据集
   IF (allocated(mtf%datasets)) THEN
      DO i = 1, mtf%ndataset
         ! 清理可变大小数组
         IF (allocated(mtf%datasets(i)%vnames)) &
            deallocate(mtf%datasets(i)%vnames)
         IF (allocated(mtf%datasets(i)%tintalgo)) &
            deallocate(mtf%datasets(i)%tintalgo)
         IF (allocated(mtf%datasets(i)%timelog)) &
            deallocate(mtf%datasets(i)%timelog)
         IF (allocated(mtf%datasets(i)%fprefix)) &
            deallocate(mtf%datasets(i)%fprefix)
         IF (allocated(mtf%datasets(i)%dtime)) &
            deallocate(mtf%datasets(i)%dtime)
         IF (allocated(mtf%datasets(i)%offset)) &
            deallocate(mtf%datasets(i)%offset)
            
         ! 清理数据存储数组
         IF (allocated(mtf%datasets(i)%forcn)) &
            deallocate(mtf%datasets(i)%forcn)
         IF (allocated(mtf%datasets(i)%forcn_LB)) &
            deallocate(mtf%datasets(i)%forcn_LB)
         IF (allocated(mtf%datasets(i)%forcn_UB)) &
            deallocate(mtf%datasets(i)%forcn_UB)
            
         ! 清理时间戳数组
         IF (allocated(mtf%datasets(i)%tstamp_LB)) &
            deallocate(mtf%datasets(i)%tstamp_LB)
         IF (allocated(mtf%datasets(i)%tstamp_UB)) &
            deallocate(mtf%datasets(i)%tstamp_UB)
            
         ! 清理掩码数组
         IF (allocated(mtf%datasets(i)%forcmask_pch)) &
            deallocate(mtf%datasets(i)%forcmask_pch)
      ENDDO
      
      deallocate(mtf%datasets)
   ENDIF
   
   ! 清理映射数组
   IF (allocated(mtf%var_to_dataset)) &
      deallocate(mtf%var_to_dataset)
   IF (allocated(mtf%var_to_index)) &
      deallocate(mtf%var_to_index)
   
   ! 重置计数器
   mtf%ndataset = 0
   
   IF (p_is_master) THEN
      WRITE(*,*) "multi_tracer_forcing清理完成!"
   ENDIF
   
   END SUBROUTINE cleanup_multi_tracer_forcing

!-----------------------------------------------------------------------
! !DESCRIPTION:
!  Convert a string to lowercase
!
! !REVISIONS:
!  2024: Added for case-insensitive string comparisons
!
!-----------------------------------------------------------------------
   SUBROUTINE to_lower(str)
      IMPLICIT NONE
      CHARACTER(len=*), INTENT(INOUT) :: str
      INTEGER :: i, str_len
      
      str_len = LEN_TRIM(str)
      DO i = 1, str_len
         IF (str(i:i) >= 'A' .AND. str(i:i) <= 'Z') THEN
            str(i:i) = CHAR(ICHAR(str(i:i)) + 32)
         ENDIF
      ENDDO
      
      ! Clear any remaining characters to avoid buffer issues
      IF (str_len < LEN(str)) THEN
         str(str_len+1:) = ' '
      ENDIF
   END SUBROUTINE to_lower

   ! ---- tracer forcing metreadLBUB ----
   SUBROUTINE tracer_metreadLBUB (idate, di)

   USE MOD_Namelist
   USE MOD_Block
   USE MOD_DataType
   USE MOD_Block
   USE MOD_NetCDFBlock
   USE MOD_RangeCheck
   IMPLICIT NONE

   integer, intent(in) :: idate(3)
   integer, intent(in) :: di  ! dataset index

   ! Local variables
   integer         :: ivar, year, month, day, time_i
   integer         :: iblkme, ib, jb, i, j
   type(timestamp) :: mtstamp
   character(len=256) :: filename

      mtstamp = idate

      DO ivar = 1, mtf%datasets(di)%NVAR

         IF (trim(mtf%datasets(di)%vnames(ivar)) == 'NULL') CYCLE     ! no data, CYCLE

         ! lower and upper boundary data already exist, CYCLE
         IF ( .not.(mtf%datasets(di)%tstamp_LB(ivar)=='NULL') .and. &
              .not.(mtf%datasets(di)%tstamp_UB(ivar)=='NULL') .and. &
            mtf%datasets(di)%tstamp_LB(ivar)<=mtstamp .and. &
            mtstamp<mtf%datasets(di)%tstamp_UB(ivar) ) THEN
            CYCLE
         ENDIF

         ! set lower boundary time stamp and get data
         IF (mtf%datasets(di)%tstamp_LB(ivar) == 'NULL') THEN
            CALL tracer_setstampLB(mtstamp, di, ivar, year, month, day, time_i)

            ! read forcing data
            filename = trim(mtf%datasets(di)%tracer_dir)//&
                      trim(tracerfilename(year, month, day, di, ivar))
            
            CALL ncio_read_block_time (filename, mtf%datasets(di)%vnames(ivar), &
                                      mtf%datasets(di)%tracegrid, time_i, mtf%datasets(di)%metdata)

            CALL block_data_copy (mtf%datasets(di)%metdata, mtf%datasets(di)%forcn_LB(ivar))
         ENDIF

         ! set upper boundary time stamp and get data
         IF (mtf%datasets(di)%tstamp_UB(ivar) == 'NULL' .or. &
             mtf%datasets(di)%tstamp_UB(ivar) <= mtstamp) THEN

            IF ( .not. (mtf%datasets(di)%tstamp_UB(ivar) == 'NULL') ) THEN
               CALL block_data_copy (mtf%datasets(di)%forcn_UB(ivar), mtf%datasets(di)%forcn_LB(ivar))
            ENDIF

            CALL tracer_setstampUB(di, ivar, year, month, day, time_i)

            ! when reaching the END of forcing data, show a Warning but still try to run
            IF ( year>mtf%datasets(di)%endyr .or. &
                 (month>mtf%datasets(di)%endmo .and. year==mtf%datasets(di)%endyr) ) THEN
               write(*,*) 'model year: ', year, 'tracer forcing end year defined: ', mtf%datasets(di)%endyr
               print *, 'Warning: reaching the END of tracer forcing data defined!'
            ENDIF

            ! read forcing data
            filename = trim(mtf%datasets(di)%tracer_dir)//&
                      trim(tracerfilename(year, month, day, di, ivar))
            
            CALL ncio_read_block_time (filename, mtf%datasets(di)%vnames(ivar), &
                                      mtf%datasets(di)%tracegrid, time_i, mtf%datasets(di)%metdata)

            CALL block_data_copy (mtf%datasets(di)%metdata, mtf%datasets(di)%forcn_UB(ivar))

            ! calculate time average coszen, for shortwave radiation (if needed for tracers)
            ! IF (ivar == 7) THEN
            !    CALL tracer_calavgcos(idate, di)
            ! ENDIF
         ENDIF

      ENDDO

   END SUBROUTINE tracer_metreadLBUB

   ! ---- tracer forcing setstampUB ----
   SUBROUTINE tracer_setstampUB(di, var_i, year, month, mday, time_i)

   IMPLICIT NONE
   integer,         intent(in)  :: di     ! dataset index
   integer,         intent(in)  :: var_i  ! variable index
   integer,         intent(out) :: year
   integer,         intent(out) :: month
   integer,         intent(out) :: mday
   integer,         intent(out) :: time_i

   integer :: day, sec
   integer :: months(0:12)

      ! calculate the time stamp
      IF ( mtf%datasets(di)%tstamp_UB(var_i) == 'NULL' ) THEN
         mtf%datasets(di)%tstamp_UB(var_i) = mtf%datasets(di)%tstamp_LB(var_i) + mtf%datasets(di)%dtime(var_i)
      ELSE
         mtf%datasets(di)%tstamp_LB(var_i) = mtf%datasets(di)%tstamp_UB(var_i)
         mtf%datasets(di)%tstamp_UB(var_i) = mtf%datasets(di)%tstamp_UB(var_i) + mtf%datasets(di)%dtime(var_i)
      ENDIF

      ! calculate initial year, day, and second values
      year = mtf%datasets(di)%tstamp_UB(var_i)%year
      day  = mtf%datasets(di)%tstamp_UB(var_i)%day
      sec  = mtf%datasets(di)%tstamp_UB(var_i)%sec

      IF ( trim(mtf%datasets(di)%groupby) == 'year' ) THEN

         ! adjust year value
         IF ( sec==86400 .and. mtf%datasets(di)%offset(var_i).eq.0 ) THEN
            sec = 0
            day = day + 1
            IF( isleapyear(year) .and. day==367) THEN
               year = year + 1; day = 1
            ENDIF
            IF( .not. isleapyear(year) .and. day==366) THEN
               year = year + 1; day = 1
            ENDIF
         ENDIF

         ! in case of leapyear with a non-leapyear calendar
         ! USE the data 1 day before after FEB 28th (Julian day 59).
         IF ( .not. mtf%datasets(di)%leapyear .and. isleapyear(year) .and. day>59 ) THEN
            day = day - 1
         ENDIF

         ! set record index
         sec = 86400*(day-1) + sec
         time_i = floor( (sec-mtf%datasets(di)%offset(var_i)) *1. / mtf%datasets(di)%dtime(var_i) ) + 1
      ENDIF

      IF ( trim(mtf%datasets(di)%groupby) == 'month' ) THEN

         IF ( isleapyear(year) ) THEN
            months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
         ELSE
            months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
         ENDIF

         ! calculate initial month and day values
         CALL julian2monthday(year, day, month, mday)

         ! record in the next day, adjust year, month and second values
         IF ( sec==86400 .and. mtf%datasets(di)%offset(var_i).eq.0 ) THEN
            sec  = 0
            mday = mday + 1
            IF ( mday > (months(month)-months(month-1)) ) THEN
               mday = 1
               IF (month == 12) THEN
                  month = 1
                  year = year + 1
               ELSE
                  month = month + 1
               ENDIF
            ENDIF
         ENDIF

         ! in case of leapyear with a non-leapyear calendar
         ! for day 29th Feb, USE the data 1 day before, i.e., 28th FEB.
         IF ( .not. mtf%datasets(di)%leapyear .and. isleapyear(year) .and. month==2 .and. mday==29 ) THEN
            mday = 28
         ENDIF

         ! set record index
         sec    = 86400*(mday-1) + sec
         time_i = floor( (sec-mtf%datasets(di)%offset(var_i)) *1. / mtf%datasets(di)%dtime(var_i) ) + 1
      ENDIF

      IF ( trim(mtf%datasets(di)%groupby) == 'day' ) THEN
         IF ( isleapyear(year) ) THEN
            months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
         ELSE
            months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
         ENDIF

         ! calculate initial month and day values
         CALL julian2monthday(year, day, month, mday)

         ! record in the next day, adjust year, month and second values
         IF ( sec==86400 .and. mtf%datasets(di)%offset(var_i).eq.0 ) THEN
            sec  = 0
            mday = mday + 1
            IF ( mday > (months(month)-months(month-1)) ) THEN
               mday = 1
               IF (month == 12) THEN
                  month = 1
                  year = year + 1
               ELSE
                  month = month + 1
               ENDIF
            ENDIF
         ENDIF

         ! in case of leapyear with a non-leapyear calendar
         ! for day 29th Feb, USE the data 1 day before, i.e., 28th FEB.
         IF ( .not. mtf%datasets(di)%leapyear .and. isleapyear(year) .and. month==2 .and. mday==29 ) THEN
            mday = 28
         ENDIF

         ! set record index
         time_i = floor( (sec-mtf%datasets(di)%offset(var_i)) *1. / mtf%datasets(di)%dtime(var_i) ) + 1
      ENDIF

      IF(time_i < 0) THEN
         write(6, *) "got the wrong time record of tracer forcing! STOP!"; CALL CoLM_stop()
      ENDIF

      RETURN

   END SUBROUTINE tracer_setstampUB

   ! ---- tracer metpreprocess ----
   SUBROUTINE tracer_metpreprocess (tracegrid, forcn, di)

   USE MOD_Block
   USE MOD_DataType
   USE MOD_SPMD_Task
   IMPLICIT NONE

   type(grid_type), intent(in) :: tracegrid
   type(block_data_real8_2d), intent(inout) :: forcn(:)
   integer, intent(in) :: di  ! dataset index

   ! Local variables
   integer :: ivar, iblkme, ib, jb, i, j
   real(r8) :: min_val, max_val
   character(len=64) :: tracer_name, var_name

   ! Get tracer name for this dataset
   tracer_name = trim(mtf%datasets(di)%tracer_name)

   ! Tracer-specific preprocessing
   SELECT CASE (trim(tracer_name))
   
   CASE ('O18')
      ! Preprocessing for O18 isotope data
      DO ivar = 1, mtf%datasets(di)%NVAR
         IF (trim(mtf%datasets(di)%vnames(ivar)) == 'NULL') CYCLE
      ENDDO
   CASE ('H2')
      ! Preprocessing for H2 isotope data
      DO ivar = 1, mtf%datasets(di)%NVAR
         IF (trim(mtf%datasets(di)%vnames(ivar)) == 'NULL') CYCLE
         ! Add H2-specific preprocessing logic here when needed
      ENDDO
      
   CASE DEFAULT
      ! Default preprocessing for other tracers
      IF (p_is_master) THEN
         WRITE(*,*) 'Warning: No specific preprocessing for tracer: ', trim(tracer_name)
      ENDIF
   END SELECT

   END SUBROUTINE tracer_metpreprocess

  

   ! ---- tracer forcing setstampLB ----
   SUBROUTINE tracer_setstampLB(mtstamp, di, var_i, year, month, mday, time_i)

   IMPLICIT NONE
   type(timestamp), intent(in)  :: mtstamp
   integer,         intent(in)  :: di     ! dataset index
   integer,         intent(in)  :: var_i  ! variable index
   integer,         intent(out) :: year
   integer,         intent(out) :: month
   integer,         intent(out) :: mday
   integer,         intent(out) :: time_i

   integer :: i, day, sec, ntime
   integer :: months(0:12)

      year = mtstamp%year
      day  = mtstamp%day
      sec  = mtstamp%sec

      mtf%datasets(di)%tstamp_LB(var_i)%year = year
      mtf%datasets(di)%tstamp_LB(var_i)%day  = day

      ! in the case of one year one file
      IF ( trim(mtf%datasets(di)%groupby) == 'year' ) THEN

         ! calculate the initial second
         sec    = 86400*(day-1) + sec
         time_i = floor( (sec-mtf%datasets(di)%offset(var_i)) *1. / mtf%datasets(di)%dtime(var_i) ) + 1
         sec    = (time_i-1)*mtf%datasets(di)%dtime(var_i) + mtf%datasets(di)%offset(var_i) - 86400*(day-1)
         mtf%datasets(di)%tstamp_LB(var_i)%sec = sec

         ! set time stamp (ststamp_LB)
         IF (sec < 0) THEN
            mtf%datasets(di)%tstamp_LB(var_i)%sec = 86400 + sec
            mtf%datasets(di)%tstamp_LB(var_i)%day = day - 1
            IF (mtf%datasets(di)%tstamp_LB(var_i)%day == 0) THEN
               mtf%datasets(di)%tstamp_LB(var_i)%year = year - 1
               IF ( isleapyear(year) ) THEN
                  mtf%datasets(di)%tstamp_LB(var_i)%day = 366
               ELSE
                  mtf%datasets(di)%tstamp_LB(var_i)%day = 365
               ENDIF
            ENDIF
         ENDIF

         ! set record info (year, time_i)
         IF ( sec<0 .or. (sec==0 .and. mtf%datasets(di)%offset(var_i).NE.0) ) THEN

            ! IF the required data just behind the first record
            ! -> set to the first record
            IF ( year==mtf%datasets(di)%startyr .and. month==mtf%datasets(di)%startmo .and. day==1 ) THEN
               sec = mtf%datasets(di)%offset(var_i)

               ! ELSE, set to one record backward
            ELSE
               sec = 86400 + sec
               day = day - 1
               IF (day == 0) THEN
                  year = year - 1
                  IF ( isleapyear(year) ) THEN
                     day = 366
                  ELSE
                     day = 365
                  ENDIF
               ENDIF
            ENDIF
         ENDIF ! ENDIF (sec <= 0)

         ! in case of leapyear with a non-leapyear calendar
         ! USE the data 1 day before after FEB 28th (Julian day 59).
         IF ( .not. mtf%datasets(di)%leapyear .and. isleapyear(year) .and. day>59 ) THEN
            day = day - 1
         ENDIF

         ! get record time index
         sec = 86400*(day-1) + sec
         time_i = floor( (sec-mtf%datasets(di)%offset(var_i)) *1. / mtf%datasets(di)%dtime(var_i) ) + 1
      ENDIF

      ! in the case of one month one file
      IF ( trim(mtf%datasets(di)%groupby) == 'month' ) THEN

         IF ( isleapyear(year) ) THEN
            months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
         ELSE
            months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
         ENDIF

         ! calculate initial month and day values
         CALL julian2monthday(year, day, month, mday)

         ! calculate initial second value
         sec    = 86400*(mday-1) + sec
         time_i = floor( (sec-mtf%datasets(di)%offset(var_i)) *1. / mtf%datasets(di)%dtime(var_i) ) + 1
         sec    = (time_i-1)*mtf%datasets(di)%dtime(var_i) + mtf%datasets(di)%offset(var_i) - 86400*(mday-1)
         mtf%datasets(di)%tstamp_LB(var_i)%sec  = sec

         ! set time stamp (ststamp_LB)
         IF (sec < 0) THEN
            mtf%datasets(di)%tstamp_LB(var_i)%sec = 86400 + sec
            mtf%datasets(di)%tstamp_LB(var_i)%day = day - 1
            IF (mtf%datasets(di)%tstamp_LB(var_i)%day == 0) THEN
               mtf%datasets(di)%tstamp_LB(var_i)%year = year - 1
               IF ( isleapyear(mtf%datasets(di)%tstamp_LB(var_i)%year) ) THEN
                  mtf%datasets(di)%tstamp_LB(var_i)%day = 366
               ELSE
                  mtf%datasets(di)%tstamp_LB(var_i)%day = 365
               ENDIF
            ENDIF
         ENDIF

         ! set record info (year, month, time_i)
         IF ( sec<0 .or. (sec==0 .and. mtf%datasets(di)%offset(var_i).ne.0) ) THEN

            ! IF just behind the first record -> set to first record
            IF ( year==mtf%datasets(di)%startyr .and. month==mtf%datasets(di)%startmo .and. mday==1 ) THEN
               sec = mtf%datasets(di)%offset(var_i)

               ! set to one record backward
            ELSE
               sec = 86400 + sec
               mday = mday - 1
               IF (mday == 0) THEN
                  month = month - 1
                  IF (month == 0) THEN
                     month = 12
                     year = year - 1
                     mday = 31
                  ELSE
                     mday = months(month) - months(month-1)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         ! in case of leapyear with a non-leapyear calendar
         ! USE the data 1 day before, i.e., FEB 28th.
         IF ( .not. mtf%datasets(di)%leapyear .and. isleapyear(year) .and. month==2 .and. mday==29 ) THEN
            mday = 28
         ENDIF

         ! get record time index
         sec = 86400*(mday-1) + sec
         time_i = floor( (sec-mtf%datasets(di)%offset(var_i)) *1. / mtf%datasets(di)%dtime(var_i) ) + 1
      ENDIF

      ! in the case of one day one file
      IF ( trim(mtf%datasets(di)%groupby) == 'day' ) THEN

         ! calculate initial month and day values
         CALL julian2monthday(year, day, month, mday)

         ! calculate initial second value
         time_i = floor( (sec-mtf%datasets(di)%offset(var_i)) *1. / mtf%datasets(di)%dtime(var_i) ) + 1
         sec    = (time_i-1)*mtf%datasets(di)%dtime(var_i) + mtf%datasets(di)%offset(var_i)
         mtf%datasets(di)%tstamp_LB(var_i)%sec  = sec

         ! set time stamp (ststamp_LB)
         IF (sec < 0) THEN
            mtf%datasets(di)%tstamp_LB(var_i)%sec = 86400 + sec
            mtf%datasets(di)%tstamp_LB(var_i)%day = day - 1
            IF (mtf%datasets(di)%tstamp_LB(var_i)%day == 0) THEN
               mtf%datasets(di)%tstamp_LB(var_i)%year = year - 1
               IF ( isleapyear(mtf%datasets(di)%tstamp_LB(var_i)%year) ) THEN
                  mtf%datasets(di)%tstamp_LB(var_i)%day = 366
               ELSE
                  mtf%datasets(di)%tstamp_LB(var_i)%day = 365
               ENDIF
            ENDIF

            IF ( year==mtf%datasets(di)%startyr .and. month==mtf%datasets(di)%startmo .and. mday==1 ) THEN
               sec = mtf%datasets(di)%offset(var_i)
            ! set to one record backward
            ELSE
               sec = 86400 + sec
               year = mtf%datasets(di)%tstamp_LB(var_i)%year
               CALL julian2monthday(mtf%datasets(di)%tstamp_LB(var_i)%year, mtf%datasets(di)%tstamp_LB(var_i)%day, month, mday)
            ENDIF
         ENDIF

         ! in case of leapyear with a non-leapyear calendar
         ! USE the data 1 day before, i.e., FEB 28th.
         IF ( .not. mtf%datasets(di)%leapyear .and. isleapyear(year) .and. month==2 .and. mday==29 ) THEN
            mday = 28
         ENDIF

         ! get record time index
         time_i = floor( (sec-mtf%datasets(di)%offset(var_i)) *1. / mtf%datasets(di)%dtime(var_i) ) + 1
      ENDIF

      IF (time_i <= 0) THEN
         write(6, *) "got the wrong time record of tracer forcing! STOP!"; CALL CoLM_stop()
      ENDIF

      RETURN

   END SUBROUTINE tracer_setstampLB

   !-----------------------------------------------------------------------
   ! !DESCRIPTION:
   !    提供对tracer forcing数据的访问接口
   !    返回指向特定tracer和变量的forcing数据的指针
   !
   ! !ARGUMENTS:
   !    tracer_idx: tracer索引 (在DEF_Tracer_Forcings_NL中的索引)
   !    var_idx:    变量索引 (在该tracer的变量列表中的索引)
   !
   ! !REVISIONS:
   !    2024: 为支持tracer accumulation模块而添加
   !-----------------------------------------------------------------------
   FUNCTION get_tracer_forcing_data(tracer_idx, var_idx) RESULT(data_ptr)

   USE MOD_Tracer_Vars_1DForcing, only: forc_prc_O18, forc_prl_O18, forc_q_O18
   IMPLICIT NONE

   integer, intent(in) :: tracer_idx  ! tracer index
   integer, intent(in) :: var_idx     ! variable index within tracer

   real(r8), pointer :: data_ptr(:)

   ! 局部变量
   character(len=64) :: tracer_name, var_name

      ! 初始化指针为空
      data_ptr => null()

      ! 检查索引有效性
      IF (tracer_idx < 1 .OR. tracer_idx > SIZE(DEF_Tracer_Forcings_NL, 1)) THEN
         RETURN
      ENDIF

      IF (var_idx < 1 .OR. var_idx > DEF_Tracer_Forcings_NL(tracer_idx)%NVAR) THEN
         RETURN
      ENDIF

      ! 获取tracer和变量名称
      tracer_name = trim(DEF_Tracer_Forcings_NL(tracer_idx)%tracer_name)
      var_name = trim(DEF_Tracer_Forcings_NL(tracer_idx)%vname(var_idx))

      ! 根据tracer类型和变量名称返回对应的数据指针
      SELECT CASE (trim(tracer_name))
      CASE ('O18')
         SELECT CASE (trim(var_name))
         CASE ('prc_O18', 'prc_o18', 'PRC_O18', 'prate1sfc')
            IF (allocated(forc_prc_O18)) THEN
               data_ptr => forc_prc_O18
            ENDIF
         CASE ('prl_O18', 'prl_o18', 'PRL_O18', 'prate2sfc')
            IF (allocated(forc_prl_O18)) THEN
               data_ptr => forc_prl_O18
            ENDIF
         CASE ('q_O18', 'q_o18', 'Q_O18', 'spfh12m')
            IF (allocated(forc_q_O18)) THEN
               data_ptr => forc_q_O18
            ENDIF
         END SELECT
      CASE ('H2')
         ! H2的数据访问可以在这里添加
         ! 当H2的forcing变量可用时
      CASE DEFAULT
         ! 对于其他tracer，可以扩展这里的逻辑
      END SELECT

   END FUNCTION get_tracer_forcing_data


   

   






END MODULE MOD_Tracer_Forcing

