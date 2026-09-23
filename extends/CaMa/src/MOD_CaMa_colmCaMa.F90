#include <define.h>
MODULE MOD_CaMa_colmCaMa
#if defined(CaMa_Flood) && (defined(GridRiverLakeFlow) || defined(CatchLateralFlow))
#error "CaMa_Flood cannot run alongside another routing solver; choose one routing scheme"
#endif
#if(defined CaMa_Flood)
!DESCRIPTION
!===========
   ! This MODULE is the coupler for the colm and CaMa-Flood model.

!ANCILLARY FUNCTIONS AND SUBROUTINES
!-------------------
   !* :SUBROUTINE:"colm_CaMa_init" :  Initialization of the coupler
   !* :SUBROUTINE:"colm_CaMa_drv"  :  Coupling between colm and CaMa-Flood
   !* :SUBROUTINE:"colm_CaMa_exit" :  Finalization of the coupler
   !* :SUBROUTINE:"get_fldinfo"    :  Get floodplain information from CaMa-Flood model
   !* :SUBROUTINE:"get_fldevp"     :  calculate floodplain evaporation

!REVISION HISTORY
!----------------
   ! 2023.02.21  Zhongwang Wei @ SYSU
   ! 2021.12.02  Zhongwang Wei @ SYSU
   ! 2020.10.01  Zhongwang Wei @ SYSU

   USE CMF_COUPLING_BUDGET_MOD
   USE MOD_Namelist
   USE MOD_CaMa_Vars
   USE PARKIND1,                  only: JPRB, JPRM, JPIM
   USE CMF_DRV_CONTROL_MOD,       only: CMF_DRV_INPUT,   CMF_DRV_INIT,    CMF_DRV_END
   USE CMF_DRV_ADVANCE_MOD,       only: CMF_DRV_ADVANCE
   USE CMF_CTRL_OUTPUT_MOD,       only: CMF_OUTPUT_INIT,CMF_OUTPUT_END,NVARSOUT,VAROUT
   USE YOS_CMF_INPUT,             only: NXIN, NYIN, DT,DTIN,IFRQ_INP,LLEAPYR,NX,NY,RMIS,DMIS
   USE MOD_Precision,             only: r8,r4
   USE YOS_CMF_INPUT ,            only: LROSPLIT,LWEVAP,LWINFILT,LDAMIRR,CSETFILE,LSEDIMENT,LCOLMFEEDBACK
   USE YOS_CMF_MAP,               only: D1LON, D1LAT
   USE YOS_CMF_INPUT,             only: WEST,EAST,NORTH,SOUTH
   USE YOS_CMF_INPUT,             ONLY: LOGNAM

   USE MOD_SPMD_Task
   USE CMF_CTRL_TIME_MOD
   USE MOD_Vars_Global,           only: spval
   USE MOD_Vars_1DFluxes
#ifdef CROP
   USE MOD_LandCrop
#endif
   USE MOD_Qsadv
   USE YOS_CMF_PROG,              only: dirrig_cama,dirrig_cama_orig,dirrig_cama_unmt,&
                                          release_cama,release_cama_riv,release_cama_dam,release_cama_rof

   USE CMF_CTRL_RESTART_MOD
   USE CMF_CTRL_SED_MOD,          only: CMF_SED_RESTART_WRITE,CMF_SED_FORCING_PUT
   USE CMF_CTRL_OUTPUT_MOD,       only: CVNAMES, CVARSOUT
   IMPLICIT NONE
   !----------------------- Dummy argument --------------------------------
   real(r8) :: coupling_elapsed=0._r8
   ! runoff volumes [m3] over the run (master): supplied by the land, and dropped
   ! because no routing cell receives it
   real(r8) :: runoff_supplied = 0._r8, runoff_unrouted = 0._r8
   logical  :: unrouted_reported = .false.
   real(r8), allocatable :: exchange_area(:,:), evap_debit(:), infil_debit(:)
   integer I,J
   real(r8),ALLOCATABLE            :: ZBUFF(:,:,:)        ! Buffer to store forcing runoff
   real(r8),ALLOCATABLE            :: ZBUFF_SED(:,:)        ! Buffer to store forcing sediment
   INTERFACE colm_CaMa_init
      MODULE PROCEDURE colm_CaMa_init
   END INTERFACE

   INTERFACE colm_CaMa_drv
      MODULE PROCEDURE colm_CaMa_drv
   END INTERFACE

   INTERFACE colm_CaMa_exit
      MODULE PROCEDURE colm_CaMa_exit
   END INTERFACE
#endif
#ifdef GridRiverLakeFlow
   ! The open-water surface flux is shared with GridRiverLakeFlow; the CaMa
   ! driver and its state remain excluded from GridRiverLakeFlow builds.
   USE MOD_Namelist, only: DEF_USE_CBL_HEIGHT
   USE MOD_Precision, only: r8
   USE MOD_Qsadv
   IMPLICIT NONE
#endif
#if defined(CaMa_Flood) || defined(GridRiverLakeFlow)
CONTAINS
#endif

#ifdef CaMa_Flood
   SUBROUTINE colm_CaMa_init(start_date)
   USE MOD_LandPatch
   USE YOS_CMF_TIME,          only: YYYY0
   USE MOD_TimeManager, only: julian2monthday, monthday2julian
   USE MOD_Utils, only: areaquad
   USE YOS_CMF_MAP, only: INPX,INPY,INPA,NSEQMAX,NSEQRIV,I1SEQX,I1SEQY,I1NEXT
   USE YOS_CMF_INPUT, only: LPTHOUT
   USE CMF_CTRL_FORCING_MOD, only: LINTERP
   USE YOS_CMF_MAP, only: I1P_OUT,I1P_OUTN,I1P_INF,I1P_INFN,PTH_UPST,PTH_DOWN

   IMPLICIT NONE
   integer, intent(in) :: start_date(3)
   !** local variables
   integer :: i,j,c,k,nextc,ncut_link,ncut_bif,spinup_julian
   integer, allocatable :: ix(:,:),iy(:,:)
   real(r8), allocatable :: w(:,:),coverage(:)

   integer(KIND=JPIM)          :: JF

      IF (p_is_master) THEN
         SELECT CASE(TRIM(DEF_UnitCatchment_file))
         CASE('', 'null', 'NULL', 'none', 'NONE')
            ! Legacy input compatibility: preserve the original CaMa configuration.
            CSETFILE = DEF_CaMa_Namelist
            SELECT CASE(TRIM(CSETFILE))
            CASE('', 'null', 'NULL', 'none', 'NONE')
               ERROR STOP 'CaMa: set DEF_UnitCatchment_file, or DEF_CaMa_Namelist for legacy inputs'
            END SELECT
            CALL CMF_DRV_INPUT
         CASE DEFAULT
            ! Same input file and bifurcation switch as GridRiverLakeFlow.
            IF(DEF_USE_LEVEE.OR.DEF_Reservoir_Method/=0) &
               ERROR STOP 'CaMa NC mode: levee/reservoir operations require CaMa-specific parameters'
            LCOLMFEEDBACK = DEF_CaMa_FloodFeedback
            CALL CMF_DRV_INPUT(DEF_UnitCatchment_file, DEF_CaMa_Restart_file)
            LPTHOUT = DEF_USE_BIFURCATION
            IF(LPTHOUT) CVARSOUT=TRIM(CVARSOUT)//',pthflw,pthout'
         END SELECT
         !get the time information from colm namelist
         DTIN     = IFRQ_INP*3600                                              ! time step of input data [sec]
         SYEAR=start_date(1)
         CALL julian2monthday(start_date(1),start_date(2),SMON,SDAY)
         CALL CMF_TIME_SET_START_SECONDS(start_date(3))
         EYEAR    = DEF_simulation_time%end_year                               ! end year
         EMON     = DEF_simulation_time%end_month                              ! end month
         EDAY     = DEF_simulation_time%end_day                                ! end day
         CALL CMF_TIME_SET_END_SECONDS(DEF_simulation_time%end_sec)
         IF(DT<=0._r8.OR.MOD(NINT(DT),60)/=0.OR.ABS(DT-NINT(DT))>1.e-8_r8) &
            ERROR STOP 'CaMa: DT must be a positive whole number of minutes'
         IF(IFRQ_INP<=0) ERROR STOP 'CaMa: invalid coupling interval'
         LLEAPYR  = DEF_forcing%leapyear                                       ! leap year flag
         LLEGACY_DAILY_HISTORY = TRIM(ADJUSTL(DEF_HIST_FREQ)) == 'DAILY'
         IF (LLEGACY_DAILY_HISTORY) THEN
            CALL monthday2julian(DEF_simulation_time%spinup_year, DEF_simulation_time%spinup_month, &
                                 DEF_simulation_time%spinup_day, spinup_julian)
            IF (start_date(1)<DEF_simulation_time%spinup_year .OR. &
                (start_date(1)==DEF_simulation_time%spinup_year .AND. &
                 (start_date(2)<spinup_julian .OR. &
                  (start_date(2)==spinup_julian .AND. start_date(3)<=DEF_simulation_time%spinup_sec)))) &
               LLEGACY_DAILY_HISTORY=.false.
         ENDIF

         CALL system('mkdir -p ' // trim(DEF_dir_restart)//'/CaMa')


         !----------------------- Dummy argument --------------------------------
         YYYY0    = SYEAR
         RMIS     = spval
         DMIS     = spval

         CALL CMF_DRV_INIT       !INITIALIZATION

         !Initialize varialbes to be outputed from variable list
         DO JF=1,NVARSOUT
            SELECT CASE (CVNAMES(JF))
            CASE ('rivout') ! river discharge [m3/s]
               DEF_hist_cama_vars%rivout=.true.
            CASE ('rivsto') ! river storage  [m3]
               DEF_hist_cama_vars%rivsto=.true.
            CASE ('rivdph') ! river depth   [m]
               DEF_hist_cama_vars%rivdph=.true.
            CASE ('rivvel') ! river velocity [m/s]
               DEF_hist_cama_vars%rivvel=.true.
            CASE ('fldout') ! floodplain discharge [m3/s]
               DEF_hist_cama_vars%fldout=.true.
            CASE ('fldsto') ! floodplain storage [m3]
               DEF_hist_cama_vars%fldsto=.true.
            CASE ('flddph') ! floodplain depth [m]
               DEF_hist_cama_vars%flddph=.true.
            CASE ('fldfrc') ! floodplain fraction (0-1)
               DEF_hist_cama_vars%fldfrc=.true.
            CASE ('fldare') ! floodplain area [m2]
               DEF_hist_cama_vars%fldare=.true.
            CASE ('sfcelv') ! water surface elevation [m]
               DEF_hist_cama_vars%sfcelv=.true.
            CASE ('totout') ! total discharge(river+floodplain) [m3/s]
               DEF_hist_cama_vars%totout=.true.
            CASE ('outflw') ! compatibility for previous file name =totout
               DEF_hist_cama_vars%outflw=.true.
            CASE ('totsto') ! total storage(river+floodplain) [m3]
               DEF_hist_cama_vars%totsto=.true.
            CASE ('storge') ! compatibility for previous file name =totsto
               DEF_hist_cama_vars%storge=.true.
            CASE ('pthflw') ! bifurcation channel discharge [m3/s]
               DEF_hist_cama_vars%pthflw=.true.
            CASE ('pthout') ! net bifurcation discharge [m3/s]
               DEF_hist_cama_vars%pthout=.true.
            CASE ('maxflw') ! daily maximum discharge [m3/s]
               DEF_hist_cama_vars%maxflw=.true.
            CASE ('maxdph') ! daily maximum depth [m]
               DEF_hist_cama_vars%maxdph=.true.
            CASE ('maxsto') ! daily maximum storage [m3]
               DEF_hist_cama_vars%maxsto=.true.
            !TODO: check the difference between gwsto and gdwsto
            CASE ('gwsto')  ! ground water storage [m3]
               DEF_hist_cama_vars%gwsto=.true.
            CASE ('gdwsto') ! ground water storage [m3]
               DEF_hist_cama_vars%gdwsto=.true.
            CASE ('gwout')  ! ground water discharge [m3/s]
               DEF_hist_cama_vars%gwout=.true.
            !TODO: check the difinition of gdwrtn, runoff, rofsfc, rofsub
            CASE ('gdwrtn') ! Ground water return flow [m3/s]
               DEF_hist_cama_vars%gdwrtn=.true.
            CASE ('runoff') ! total runoff [m3/s]                !!  compatibility for previous file name
               DEF_hist_cama_vars%runoff=.true.
            CASE ('runoffsub') ! subsurface runoff [m3/s]        !!  compatibility for previous file name
               DEF_hist_cama_vars%runoffsub=.true.
            CASE ('rofsfc') ! surface runoff [m3/s]              !!  compatibility for previous file name
               DEF_hist_cama_vars%rofsfc=.true.
            CASE ('rofsub')  ! input sub-surface runoff [m3/s]
               DEF_hist_cama_vars%rofsub=.true.
            CASE ('damsto')   ! reservoir storage [m3]
               DEF_hist_cama_vars%damsto=.true.
            CASE ('daminf')   ! reservoir inflow [m3/s]
                  DEF_hist_cama_vars%daminf=.true.
            CASE ('levsto')   !flood storage in protected side (storage betwen river & levee) [m3]
                  DEF_hist_cama_vars%levsto=.true.
            CASE ('levdph')   !flood depth in protected side [m]
                  DEF_hist_cama_vars%levdph=.true.
            CASE ('wevap')    ! input inundation Evaporation [m]
               IF (LWEVAP) THEN
                  DEF_hist_cama_vars%wevap=.true.
               ELSE
                  DEF_hist_cama_vars%wevap=.false.
               ENDIF
            CASE ('winfilt')  ! input inundation re-infiltrition [m]
               IF (LWINFILT) THEN
                  DEF_hist_cama_vars%winfilt=.true.
               ELSE
                  DEF_hist_cama_vars%winfilt=.false.
               ENDIF
            CASE ('outflw_ocean') ! discharge to ocean [m3/s]
               DEF_hist_cama_vars%outflw_ocean=.true.
            CASE ('sedout') ! suspended sediment flow [m3/s]
               IF (LSEDIMENT) THEN
                  DEF_hist_cama_vars%sedout=.true.
               ELSE
                  DEF_hist_cama_vars%sedout=.false.
               ENDIF
            CASE ('sedcon') ! suspended sediment concentration [m3/m3]
               IF (LSEDIMENT) THEN
                  DEF_hist_cama_vars%sedcon=.true.
               ELSE
                  DEF_hist_cama_vars%sedcon=.false.
               ENDIF
            CASE ('sedinp') ! sediment inflow from land [m3/s]
               IF (LSEDIMENT) THEN
                  DEF_hist_cama_vars%sedinp=.true.
               ELSE
                  DEF_hist_cama_vars%sedinp=.false.
               ENDIF
            CASE ('bedout') ! bedload [m3/s]
               IF (LSEDIMENT) THEN
                  DEF_hist_cama_vars%bedout=.true.
               ELSE
                  DEF_hist_cama_vars%bedout=.false.
               ENDIF
            CASE ('netflw') ! net entrainment flow [m3/s]
               IF (LSEDIMENT) THEN
                  DEF_hist_cama_vars%netflw=.true.
               ELSE
                  DEF_hist_cama_vars%netflw=.false.
               ENDIF
            CASE ('layer') ! exchange layer volume [m3]
               IF (LSEDIMENT) THEN
                  DEF_hist_cama_vars%layer=.true.
               ELSE
                  DEF_hist_cama_vars%layer=.false.
               ENDIF
            CASE ('deplyr') ! river bed volume (vertical layer) [m3]
               IF (LSEDIMENT) THEN
                  ! Note: deplyr is a special case for sediment deposition layers
                  ! It requires additional parsing to determine layer number
                  DEF_hist_cama_vars%layer=.true.  ! Use existing layer flag
               ELSE
                  DEF_hist_cama_vars%layer=.false.
               ENDIF
            CASE DEFAULT
               STOP
            END SELECT
         ENDDO

      ENDIF

      !Broadcast the variables to all the processors
      CALL mpi_bcast (NX      ,   1, MPI_INTEGER,   p_address_master, p_comm_glb, p_err) ! number of grid points in x-direction of CaMa-Flood
      CALL mpi_bcast (NY      ,   1, MPI_INTEGER,   p_address_master, p_comm_glb, p_err) ! number of grid points in y-direction of CaMa-Flood
      CALL mpi_bcast (IFRQ_INP,   1, MPI_INTEGER,   p_address_master, p_comm_glb, p_err) ! input frequency of CaMa-Flood (hour)
      CALL mpi_bcast (LWEVAP  ,   1, MPI_LOGICAL,   p_address_master, p_comm_glb, p_err) ! switch for inundation evaporation
      CALL mpi_bcast (LWINFILT,   1, MPI_LOGICAL,   p_address_master, p_comm_glb, p_err) ! switch for inundation re-infiltration
      CALL mpi_bcast (LDAMIRR ,   1, MPI_LOGICAL,   p_address_master, p_comm_glb, p_err) ! switch for deficit irrigation
      CALL mpi_bcast (LSEDIMENT,  1, MPI_LOGICAL,   p_address_master, p_comm_glb, p_err) ! switch for sediment transport

      IF (.not. allocated(D1LAT))  allocate (D1LAT(NY))
      IF (.not. allocated(D1LON))  allocate (D1LON(NX))

      CALL mpi_bcast (D1LAT, NY, MPI_REAL8, p_address_master, p_comm_glb, p_err) !
      CALL mpi_bcast (D1LON, NX, MPI_REAL8, p_address_master, p_comm_glb, p_err) !
      CALL mpi_bcast (SOUTH,  1, MPI_REAL8, p_address_master, p_comm_glb, p_err) !
      CALL mpi_bcast (NORTH,  1, MPI_REAL8, p_address_master, p_comm_glb, p_err) !
      CALL mpi_bcast (WEST ,  1, MPI_REAL8, p_address_master, p_comm_glb, p_err) !
      CALL mpi_bcast (EAST ,  1, MPI_REAL8, p_address_master, p_comm_glb, p_err) !

      !allocate the data structure for cama
      CALL gcama%define_by_center (D1LAT,D1LON,SOUTH,NORTH,WEST,EAST) !define the grid for cama
      CALL mp2g_cama%build_arealweighted (gcama, landpatch) !build the mapping between cama and mpi
      CALL mg2p_cama%build_arealweighted (gcama, landpatch)
      CALL cama_gather%set (gcama)

      !allocate the cama-flood related variable for accumulation
      CALL allocate_2D_cama_Fluxes  (gcama) !allocate the 2D variables
      CALL allocate_acc_cama_Fluxes () !allocate the accumulation variables
      CALL FLUSH_acc_cama_fluxes    () !initialize the accumulation variables

      !Only master processor allocate the 2D variables
      IF (p_is_master) THEN
         allocate (runoff_2d (NX,NY))
         allocate (prcp_2d (NX,NY))
         allocate (fevpg_2d  (NX,NY))
         allocate (finfg_2d  (NX,NY))
         !allocate data buffer for input forcing, flood fraction and flood depth
         allocate (ZBUFF(NX,NY,2))
         allocate (ZBUFF_SED(NX,NY))
         allocate (fldfrc_tmp(NX,NY))
         allocate (flddepth_tmp(NX,NY))
         allocate (dirrig_2d (NX,NY))
         allocate (dirrig_2d_orig (NX,NY))
         allocate (dirrig_cama (NX,NY))
         allocate (dirrig_cama_orig(NX,NY))
         allocate (dirrig_cama_unmt(NX,NY))
         allocate (withdrawal_tmp(NX,NY))
         allocate (withdrawal_riv_tmp(NX,NY))
         allocate (withdrawal_dam_tmp(NX,NY))
         allocate (withdrawal_rof_tmp(NX,NY))
         allocate (release_cama(NX,NY))
         allocate (release_cama_riv(NX,NY))
         allocate (release_cama_dam(NX,NY))
         allocate (release_cama_rof(NX,NY))
         !Initialize the data buffer for input forcing, flood fraction and flood depth
         runoff_2d(:,:)    = 0.0 !runoff in master processor
         prcp_2d(:,:)      = 0.0
         fevpg_2d(:,:)     = 0.0 !evaporation in master processor
         finfg_2d(:,:)     = 0.0 !re-infiltration in master processor
         ZBUFF(:,:,:)      = 0.0 !input forcing in master processor
         ZBUFF_SED(:,:)    = 0.0 !input forcing in master processor

         fldfrc_tmp(:,:)   = 0.0 !flood fraction in master processor
         flddepth_tmp(:,:) = 0.0 !flood depth in master processor
         dirrig_2d(:,:)    = 0.0 !deficit_irrigation in master processor
         dirrig_2d_orig(:,:)     = 0.0 !deficit_irrigation in master processor
         dirrig_cama(:,:)        = 0.0 !deficit_irrigation in master processor
         dirrig_cama_orig(:,:)   = 0.0 !deficit_irrigation in master processor
         dirrig_cama_unmt(:,:)   = 0.0 !deficit_irrigation in master processor
         release_cama(:,:)       = 0.0 !release in master processor
         release_cama_riv(:,:)   = 0.0 !release in master processor
         release_cama_dam(:,:)   = 0.0 !release in master processor
         release_cama_rof(:,:)   = 0.0 !release in master processor
         withdrawal_tmp(:,:)     = 0.0 !withdrawal in master processor
         withdrawal_riv_tmp(:,:) = 0.0 !withdrawal in master processor
         withdrawal_dam_tmp(:,:) = 0.0 !withdrawal in master processor
         withdrawal_rof_tmp(:,:) = 0.0 !withdrawal in master processor
      ENDIF
      !allocate the cama-flood related variable in worker processors
      IF (p_is_worker) THEN
         allocate (flddepth_cama(numpatch)) !flood depth in worker processors
         allocate (fldfrc_cama(numpatch))   !flood fraction in worker processors
         allocate (fevpg_fld(numpatch))     !evaporation in worker processors
         allocate (finfg_fld(numpatch))     !re-infiltration in worker processors
         allocate (dirrig_tmp(numpatch))    ! demand irrigation in worker processors
         allocate (dirrig_day(numpatch))    ! demand irrigation(day) in worker processors
         allocate (withdrawal_cama(numpatch))     !withdrawal in worker processors
         allocate (withdrawal_riv_cama(numpatch)) !withdrawal in worker processors
         allocate (withdrawal_dam_cama(numpatch)) !withdrawal in worker processors
         allocate (withdrawal_rof_cama(numpatch)) !withdrawal in worker processors

         flddepth_cama(:)     =  0.0
         fldfrc_cama(:)       =  0.0
         fevpg_fld(:)         =  0.0
         finfg_fld(:)         =  0.0
         dirrig_tmp(:)        =  0.0
         dirrig_day(:)        =  0.0
         withdrawal_cama(:)   =  0.0
         withdrawal_riv_cama(:) =  0.0
         withdrawal_dam_cama(:) =  0.0
         withdrawal_rof_cama(:) =  0.0
      ENDIF
      CALL mg2p_cama%allocate_part(flood_credit_part)
      CALL mg2p_cama%allocate_part(flood_sink_part)
      IF(p_is_worker) ALLOCATE(flood_credit(numpatch))
      IF(p_is_master) THEN
         ALLOCATE(exchange_area(NX,NY),coverage(NSEQMAX),evap_debit(NSEQMAX),infil_debit(NSEQMAX))
         DO j=1,NY
            DO i=1,NX
               exchange_area(i,j)=areaquad(gcama%lat_s(j),gcama%lat_n(j),gcama%lon_w(i),gcama%lon_e(i))*1.e6_r8
            ENDDO
         ENDDO
         IF(LINTERP) THEN
            IF(NXIN/=NX.OR.NYIN/=NY) ERROR STOP 'CaMa: coupled input grid must match routing grid coordinates'
            ! Legacy binary matrices may carry missing-value area sentinels
            ! behind an inactive x index. Never change the input arrays.
            w=INPA
            WHERE(INPX<=0) w=0._r8
            CALL budget_init(INPX,INPY,w,exchange_area)
         ELSE
            ALLOCATE(ix(NSEQMAX,1),iy(NSEQMAX,1),w(NSEQMAX,1))
            ix(:,1)=I1SEQX; iy(:,1)=I1SEQY; w=1._r8
            CALL budget_init(ix,iy,w,exchange_area)
         ENDIF
      ENDIF
      ! Diagnose regional coverage without extrapolating runoff over missing land.
      IF(p_is_worker) THEN
         a_rnof_cama=1._r8
         nacc=1._r8
      ENDIF
      CALL colm2cama_real8(a_rnof_cama,f_rnof_cama,runoff_2d,integral=.true.)
      IF(p_is_master) THEN
         ! DEF_CaMa_StrictDomain = .false. (default): a regional domain does not have to
         ! contain complete drainage basins.  Runoff is then routed only where a
         ! routing cell receives it (the rest is counted and reported), and water
         ! entering across the domain boundary is absent.  Strict: stop instead.
         CALL budget_set_strict(DEF_CaMa_StrictDomain)
         CALL budget_runoff(runoff_2d,coverage)
         ncut_link = 0
         DO c=1,NSEQRIV
            nextc=I1NEXT(c)
            IF(nextc>0) THEN
               IF((coverage(c)>0._r8).NEQV.(coverage(nextc)>0._r8)) THEN
                  IF(DEF_CaMa_StrictDomain) &
                     ERROR STOP 'CaMa: regional domain cuts a river link; supply complete basins or external boundaries'
                  ncut_link = ncut_link + 1
               ENDIF
            ENDIF
         ENDDO
         ! Bifurcations can reverse: both connected catchments need forcing.
         ncut_bif = 0
         IF(LPTHOUT) THEN
            DO c=1,NSEQMAX
               IF(coverage(c)<=0._r8) CYCLE
               DO k=1,I1P_OUTN(c)
                  IF(coverage(PTH_DOWN(I1P_OUT(c,k)))<=0._r8) THEN
                     IF(DEF_CaMa_StrictDomain) ERROR STOP 'CaMa: regional domain cuts a bidirectional bifurcation'
                     ncut_bif = ncut_bif + 1
                  ENDIF
               ENDDO
               DO k=1,I1P_INFN(c)
                  IF(coverage(PTH_UPST(I1P_INF(c,k)))<=0._r8) THEN
                     IF(DEF_CaMa_StrictDomain) ERROR STOP 'CaMa: regional domain cuts a bidirectional bifurcation'
                     ncut_bif = ncut_bif + 1
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
         IF(ncut_link > 0 .OR. ncut_bif > 0) write(*,'(A,I0,A,I0,A)') &
            ' WARNING: CaMa regional domain cuts ', ncut_link, ' river links and ', ncut_bif, &
            ' bifurcation ends; inflow across the domain boundary is not supplied (DEF_CaMa_StrictDomain = .true. stops here).'
      ENDIF
      CALL flush_acc_cama_fluxes()
      coupling_elapsed=0._r8
      ! Publish immediately, including the first land step after a restart.
      IF(LWINFILT.OR.LWEVAP) CALL publish_flood_credit()
   END SUBROUTINE colm_CaMa_init

!####################################################################
   SUBROUTINE colm_cama_drv(idate_sec,deltim,force_flush)
   USE YOS_CMF_PROG, only: D2RUNOFF,D2ROFSUB,D2WEVAP,D2WINFILT,P2FLDSTO
   USE YOS_CMF_DIAG, only: D2WEVAPEX,D2WINFILTEX
   USE CMF_CTRL_PHYSICS_MOD, only: CMF_PHYSICS_FLDSTG
   IMPLICIT NONE
   integer, intent(in)  :: idate_sec     ! calendar (year, julian day, seconds)
   real(r8), intent(in) :: deltim
   logical, intent(in) :: force_flush
   real(r8) :: dt_saved,dtin_saved,remaining
      IF(deltim<=0._r8.OR.MOD(NINT(deltim),60)/=0.OR.ABS(deltim-NINT(deltim))>1.e-8_r8) &
         ERROR STOP 'CaMa: land time step must be a positive whole number of minutes'
      coupling_elapsed=coupling_elapsed+deltim
      !Accumulate cama-flood related flux variables
      CALL accumulate_cama_fluxes(deltim)
      ! If the time is the same as the input time step of cama-flood
      IF (MOD(idate_sec,3600*int(IFRQ_INP))==0 .OR. coupling_elapsed>=3600._r8*IFRQ_INP .OR. force_flush) THEN
         ! Prepare sending the accumulated runoff flux varilble to cama model (master processor to worker processors)
         CALL colm2cama_real8 (a_rnof_cama, f_rnof_cama, runoff_2d, integral=.true.)
         IF (LSEDIMENT) THEN
            CALL colm2cama_real8 (a_prcp_cama, f_prcp_cama, prcp_2d, valid_time=a_prcp_time)
         ENDIF
         IF (LDAMIRR) THEN
            ! Prepare sending the accumulated deficit irriggation varilble to cama model (master processor to worker processors)
            IF (p_is_worker) THEN
               dirrig_tmp = a_dirrig_cama
               IF (idate_sec == 3600*int(IFRQ_INP)) THEN
                  dirrig_day = 0.0
               ENDIF
               dirrig_day = dirrig_day + dirrig_tmp
            ENDIF
            CALL colmvar2cama_real8 (a_dirrig_cama, f_dirrig_cama, dirrig_2d)
         ENDIF

         ! Prepare sending the accumulated inundation evaporation flux to cama model (master processor to worker processors)
         ! only if the inundation evaporation is turned on
         IF (LWEVAP) THEN
            CALL colm2cama_real8 (a_fevpg_fld, f_fevpg_fld, fevpg_2d, integral=.true., flood_sink=.true.)
         ENDIF
         ! Prepare sending the accumulated inundation re-infiltrition flux to cama model (master processor to worker processors)
         ! only if the inundation re-infiltrition is turned on
         IF (LWINFILT) THEN
            CALL colm2cama_real8 (a_finfg_fld, f_finfg_fld, finfg_2d, integral=.true., flood_sink=.true.)
         ENDIF


      ! Reset the accumulation variables
         CALL flush_acc_cama_fluxes

         ! Initialize the variables unit for cama-flood input
         IF(p_is_master)THEN
            dt_saved=DT; dtin_saved=DTIN
            DTIN=coupling_elapsed
            ZBUFF(:,:,1)=runoff_2d ! already integrated over covered land [m3/s]
            ZBUFF(:,:,2)=0._r8
            IF (LSEDIMENT) THEN
               ZBUFF_SED(:,:) = prcp_2d(:,:) !* 1.0D-3 ! mm/s
            ENDIF

            IF (LDAMIRR) THEN
               DO i = 1,NX
                  DO j = 1, NY
                     IF (dirrig_2d(i,j).LE.0.0) THEN
                        dirrig_cama(i,j) = 0.0
                     ELSEIF (dirrig_2d(i,j).GT.0.0) THEN
                        dirrig_cama(i,j) = dirrig_2d(i,j)/1000.0 ! kg -> m3
                     ENDIF


                  ENDDO
               ENDDO

               !!!! save original dirrig_cama
               IF (idate_sec == 3600*int(IFRQ_INP)) THEN
                     dirrig_cama_orig(:,:) = 0.0
               ENDIF
               !!!! update dirrig_cama
               dirrig_cama(:,:) = dirrig_cama(:,:) + dirrig_cama_unmt(:,:)
               !!!! save initial dirrig_cama
               dirrig_cama_orig(:,:) = dirrig_cama_orig(:,:) + dirrig_cama(:,:)

               !!!! initialize
               release_cama_riv(:,:) = 0.0
               release_cama_dam(:,:) = 0.0
               release_cama_rof(:,:) = 0.0
               release_cama(:,:) = 0.0

               !!!!!withdraw water from runoff
               release_cama_rof(:,:) = min(dirrig_cama(:,:), ZBUFF(:,:,1)*DTIN)
               release_cama_rof(:,:) = max(0.0, release_cama_rof(:,:))
               dirrig_cama(:,:) = dirrig_cama(:,:) - release_cama_rof(:,:)
               ZBUFF(:,:,1) = (ZBUFF(:,:,1)*DTIN - release_cama_rof(:,:))/DTIN
            ENDIF

            ! Simulating the hydrodynamics in continental-scale rivers
            ! ----------------------------------------------------------------------

            ! Land has already consumed these credits. Pay them before any
            ! hydrodynamic outflow can spend the same reservoir water.
            IF(LWINFILT.OR.LWEVAP) THEN
               CALL budget_debit(fevpg_2d*DTIN,finfg_2d*DTIN,P2FLDSTO(:,1),evap_debit,infil_debit)
               IF(LWEVAP) D2WEVAPEX(:,1)=evap_debit/DTIN
               IF(LWINFILT) D2WINFILTEX(:,1)=infil_debit/DTIN
               sinks_prepaid=.true.
               CALL CMF_PHYSICS_FLDSTG()
            ENDIF
            ! Do not reinterpret a regional mean as a full-catchment flux.
            CALL budget_runoff(ZBUFF(:,:,1),D2RUNOFF(:,1))
            runoff_supplied = runoff_supplied + SUM(ZBUFF(:,:,1))*DTIN
            runoff_unrouted = runoff_unrouted + budget_unrouted()*DTIN
            IF(.NOT.unrouted_reported .AND. ABS(budget_unrouted())>1.e-12_r8) THEN
               unrouted_reported = .true.
               write(*,'(A,ES10.3,A)') ' WARNING: CaMa cannot route runoff of grid cells that no routing cell receives (', &
                  budget_unrouted(), ' m3/s in this window); it is dropped and totalled at the end of the run.'
            ENDIF
            D2ROFSUB=0._r8
            IF(LWEVAP) D2WEVAP=0._r8
            IF(LWINFILT) D2WINFILT=0._r8
            IF (LSEDIMENT) CALL CMF_SED_FORCING_PUT(real(ZBUFF_SED,kind=JPRB))
            remaining=DTIN
            DO WHILE(remaining>1.e-8_r8)
               DT=MIN(dt_saved,remaining)
               CALL CMF_DRV_ADVANCE(1, stop_at_end=.false.)
               remaining=remaining-DT
            ENDDO
            DT=dt_saved; DTIN=dtin_saved
            sinks_prepaid=.false.
         ENDIF
         coupling_elapsed=0._r8
         IF(LWINFILT.OR.LWEVAP) CALL publish_flood_credit()

         IF (LDAMIRR) THEN
            IF (p_is_master) THEN
               DO i = 1, NX
                  DO j = 1, NY
                     IF (dirrig_cama_orig(i,j).GT.0) then
                        release_cama(i,j) = release_cama_riv(i,j) + release_cama_dam(i,j) + release_cama_rof(i,j)
                        withdrawal_tmp(i,j) = release_cama(i,j)/dirrig_cama_orig(i,j)
                        withdrawal_riv_tmp(i,j) = release_cama_riv(i,j)/dirrig_cama_orig(i,j)
                        withdrawal_dam_tmp(i,j) = release_cama_dam(i,j)/dirrig_cama_orig(i,j)
                        withdrawal_rof_tmp(i,j) = release_cama_rof(i,j)/dirrig_cama_orig(i,j)
                     ELSE
                        withdrawal_tmp(i,j) = 0.0
                        withdrawal_riv_tmp(i,j) = 0.0
                        withdrawal_dam_tmp(i,j) = 0.0
                        withdrawal_rof_tmp(i,j) = 0.0
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF

            IF (p_is_worker) THEN
               withdrawal_cama = dirrig_day
               withdrawal_riv_cama = dirrig_day
               withdrawal_dam_cama = dirrig_day
               withdrawal_rof_cama = dirrig_day
            ENDIF

            CALL camavar2colm_real8 (withdrawal_tmp, f_withdrawal_cama, withdrawal_cama)
            CALL camavar2colm_real8 (withdrawal_riv_tmp, f_withdrawal_riv_cama, withdrawal_riv_cama)
            CALL camavar2colm_real8 (withdrawal_dam_tmp, f_withdrawal_dam_cama, withdrawal_dam_cama)
            CALL camavar2colm_real8 (withdrawal_rof_tmp, f_withdrawal_rof_cama, withdrawal_rof_cama)
         ENDIF
      ENDIF
   END SUBROUTINE colm_cama_drv

   SUBROUTINE colm_cama_exit
   USE YOS_CMF_MAP, only: I2NEXTX

   IMPLICIT NONE

      IF(p_is_master .AND. runoff_unrouted /= 0._r8) write(*,'(A,ES10.3,A,F7.3,A)') &
         ' CaMa runoff dropped (no routing cell receives it): ', runoff_unrouted, ' m3 = ', &
         100._r8*runoff_unrouted/MAX(ABS(runoff_supplied),1.e-30_r8), '% of the land runoff'
      CALL mg2p_cama%deallocate_part(flood_credit_part)
      CALL mg2p_cama%deallocate_part(flood_sink_part)
      IF(ALLOCATED(flood_credit)) DEALLOCATE(flood_credit)
      IF(ALLOCATED(exchange_area)) DEALLOCATE(exchange_area,evap_debit,infil_debit)
      ! finalize CaMa-Flood
      CALL deallocate_acc_cama_Fluxes ()
      IF(p_is_master)THEN
         ! finalize CaMa-Flood
         deallocate(ZBUFF)
         deallocate(ZBUFF_SED)
         deallocate (runoff_2d)
         deallocate (fevpg_2d)
         deallocate (finfg_2d)
         deallocate ( I2NEXTX )
         deallocate (dirrig_2d)
         deallocate (dirrig_2d_orig)
         deallocate (dirrig_cama)
         deallocate (dirrig_cama_orig)
         deallocate (withdrawal_tmp)
         deallocate (withdrawal_riv_tmp)
         deallocate (withdrawal_dam_tmp)
         deallocate (withdrawal_rof_tmp)
         deallocate (release_cama)
         deallocate (release_cama_riv)
         deallocate (release_cama_dam)
         deallocate (release_cama_rof)
      ENDIF
      IF (p_is_worker) THEN
         deallocate (flddepth_cama)
         deallocate (fldfrc_cama)
         deallocate (fevpg_fld)
         deallocate (finfg_fld)
         deallocate (dirrig_tmp)
         deallocate (dirrig_day)
         deallocate (withdrawal_cama)
         deallocate (withdrawal_riv_cama)
         deallocate (withdrawal_dam_cama)
         deallocate (withdrawal_rof_cama)
      ENDIF
   END SUBROUTINE colm_cama_exit


   SUBROUTINE colm_cama_write_restart(idate, lc_year, site, dir_restart)

   IMPLICIT NONE
   integer, intent(in) :: idate(3)
   integer, intent(in) :: lc_year      !year of land cover type data
   character(LEN=*), intent(in) :: site
   character(LEN=*), intent(in) :: dir_restart
   ! Local variables
   character(LEN=256) :: file_restart
   character(len=14)  :: cdate
   character(len=256) :: cyear         !character for lc_year

      IF(coupling_elapsed>1.e-8_r8) ERROR STOP 'CaMa: restart requires a flushed coupling window'
      ! land cover type year
      write(cyear,'(i4.4)') lc_year
      write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') idate(1), idate(2), idate(3)
      CRESTDIR =    trim(DEF_dir_restart)// '/CaMa'//'/'//trim(cdate)//'/'
      CALL system('mkdir -p ' // trim(CRESTDIR))
      CALL CMF_RESTART_WRITE()
      IF (LSEDIMENT) THEN
         CALL CMF_SED_RESTART_WRITE()
      ENDIF

   END SUBROUTINE colm_cama_write_restart

!####################################################################
   SUBROUTINE get_fldinfo()
   USE YOS_CMF_MAP, only: D2GRAREA
   USE YOS_CMF_DIAG, only: D2FLDFRC
   USE YOS_CMF_PROG, only: P2FLDSTO
      CALL budget_publish(P2FLDSTO(:,1),D2FLDFRC(:,1)*D2GRAREA(:,1),flddepth_tmp,fldfrc_tmp)
   END SUBROUTINE get_fldinfo

   SUBROUTINE publish_flood_credit()
      IF(p_is_master) CALL get_fldinfo()
      CALL cama2colm_real8(flddepth_tmp,f_flddepth_cama,flddepth_cama)
      CALL mg2p_cama%grid2part(f_flddepth_cama,flood_credit_part)
      CALL cama2colm_real8(fldfrc_tmp,f_fldfrc_cama,fldfrc_cama)
      IF(p_is_worker) THEN
         flood_credit=MAX(0._r8,flddepth_cama) ! initial patch-mean credit [m]
         fldfrc_cama=MIN(1._r8,MAX(0._r8,fldfrc_cama))
         WHERE(fldfrc_cama>0._r8)
            flddepth_cama=flood_credit/fldfrc_cama*1000._r8 ! wet-area depth [mm]
         ELSEWHERE
            flddepth_cama=0._r8
         END WHERE
      ENDIF
   END SUBROUTINE publish_flood_credit
#endif

#if defined(CaMa_Flood) || defined(GridRiverLakeFlow)
   SUBROUTINE get_fldevp (hu,ht,hq,us,vs,tm,qm,rhoair,psrf,tssea,&
	   hpbl, &
      taux,tauy,fseng,fevpg,tref,qref,z0m,zol,rib,ustar,qstar,tstar,fm,fh,fq)
!DESCRIPTION
!===========
   ! This subrountine compute surface fluxes, derviatives, and exchange coefficiants
   ! This is the main SUBROUTINE to execute the calculation of thermal processes
   ! and surface fluxes

!ANCILLARY FUNCTIONS AND SUBROUTINES
!-------------------
   !* :SUBROUTINE:"qsadv"         :     !  computes saturation mixing ratio and change in saturation
   !* :SUBROUTINE:" moninobukini" :     !  initialzation of Monin-Obukhov length, see MOD_FrictionVelocity.F90
   !* :SUBROUTINE:" moninobuk"    :     !  calculation of friction velocity, relation for potential temperature
                                        !  and humidity profiles of surface boundary layer,see MOD_FrictionVelocity.F90
!REVISION HISTORY
!----------------
   ! 2023.05.05  Shaofeng Liu @ SYSU:
   !             add option to CALL moninobuk_leddy, the LargeEddy
   !             surface turbulence scheme (LZD2022); make a proper update of um.
   ! 2020.10.01  Zhongwang Wei @ SYSU
   ! 2002.08.30  Yongjiu Dai   @ BNU
   ! 1999.09.15  Yongjiu Dai   @ BNU
   USE MOD_Precision
   USE MOD_Const_Physical, only: cpair,rgas,vonkar,grav,tfrz
   USE MOD_FrictionVelocity
   USE MOD_TurbulenceLEddy
   IMPLICIT NONE
   real(r8), intent(in)  :: hu      ! agcm reference height of wind [m]
   real(r8), intent(in)  :: ht      ! agcm reference height of temperature [m]
   real(r8), intent(in)  :: hq      ! agcm reference height of humidity [m]
   real(r8), intent(in)  :: us      ! wind component in eastward direction [m/s]
   real(r8), intent(in)  :: vs      ! wind component in northward direction [m/s]
   real(r8), intent(in)  :: tm      ! temperature at agcm reference height [kelvin]
   real(r8), intent(in)  :: qm      ! specific humidity at agcm reference height [kg/kg]
   real(r8), intent(in)  :: rhoair  ! density air [kg/m3]
   real(r8), intent(in)  :: psrf    ! atmosphere pressure at the surface [pa] [not used]
   real(r8), intent(in)  :: tssea   ! inundation surface temperature [K]-->set to tgrnd
   real(r8), intent(in)  :: hpbl    ! atmospheric boundary layer height [m]
   real(r8), intent(out) :: taux    ! wind stress: E-W [kg/m/s**2]
   real(r8), intent(out) :: tauy    ! wind stress: N-S [kg/m/s**2]
   real(r8), intent(out) :: fseng   ! sensible heat flux from ground [mm/s]
   real(r8), intent(out) :: fevpg   ! evaporation heat flux from ground [mm/s]
   real(r8), intent(out) :: tref    ! 2 m height air temperature [kelvin]
   real(r8), intent(out) :: qref    ! 2 m height air humidity [?]
   real(r8), intent(out) :: z0m     ! effective roughness [m]
   real(r8), intent(out) :: zol     ! dimensionless height (z/L) used in Monin-Obukhov theory
   real(r8), intent(out) :: rib     ! bulk Richardson number in surface layer
   real(r8), intent(out) :: ustar   ! friction velocity [m/s]
   real(r8), intent(out) :: tstar   ! temperature scaling parameter
   real(r8), intent(out) :: qstar   ! moisture scaling parameter
   real(r8), intent(out) :: fm      ! integral of profile function for momentum
   real(r8), intent(out) :: fh      ! integral of profile function for heat
   real(r8), intent(out) :: fq      ! integral of profile function for moisture

   !----------------------- Dummy argument --------------------------------
   integer i
   integer niters            ! maximum number of iterations for surface temperature
   integer iter              ! iteration index
   integer nmozsgn           ! number of times moz changes sign

   real(r8) :: beta          ! coefficient of conective velocity [-]
   real(r8) :: displax       ! zero-displacement height [m]
   real(r8) :: dth           ! diff of virtual temp. between ref. height and surface
   real(r8) :: dqh           ! diff of humidity between ref. height and surface
   real(r8) :: dthv          ! diff of vir. poten. temp. between ref. height and surface
   real(r8) :: eg            ! water vapor pressure at temperature T [Pa]
   real(r8) :: degdT         ! d(eg)/dT
   real(r8) :: obu           ! monin-obukhov length [m]
   real(r8) :: obuold        ! monin-obukhov length from previous iteration
   real(r8) :: qsatg         ! ground saturated specific humidity [kg/kg]
   real(r8) :: qsatgdT       ! d(qsatg)/dT
   real(r8) :: ram           ! aerodynamical resistance [s/m]
   real(r8) :: rah           ! thermal resistance [s/m]
   real(r8) :: raw           ! moisture resistance [s/m]
   real(r8) :: raih          ! temporary variable [kg/m2/s]
   real(r8) :: raiw          ! temporary variable [kg/m2/s]
   real(r8) :: fh2m          ! relation for temperature at 2m
   real(r8) :: fq2m          ! relation for specific humidity at 2m
   real(r8) :: fm10m         ! integral of profile function for momentum at 10m
   real(r8) :: thm           ! intermediate variable (tm+0.0098*ht)
   real(r8) :: th            ! potential temperature (kelvin)
   real(r8) :: thv           ! virtual potential temperature (kelvin)
   real(r8) :: thvstar       ! virtual potential temperature scaling parameter
   real(r8) :: um            ! wind speed including the stablity effect [m/s]
   real(r8) :: ur            ! wind speed at reference height [m/s]
   real(r8) :: visa          ! kinematic viscosity of dry air [m2/s]
   real(r8) :: wc            ! convective velocity [m/s]
   real(r8) :: wc2           ! wc**2
   real(r8) :: xt            !  ---->  temporary variables
   real(r8) :: xq            !  ---->  temporary variables
   real(r8) :: zii           ! convective boundary height [m]
   real(r8) :: zldis         ! reference height "minus" zero displacement heght [m]
   real(r8) :: z0mg          ! roughness length over ground, momentum [m]
   real(r8) :: z0hg          ! roughness length over ground, sensible heat [m]
   real(r8) :: z0qg          ! roughness length over ground, latent heat [m]

   real, parameter :: zsice = 0.04  ! sea ice aerodynamic roughness length [m]

      !-----------------------------------------------------------------------
      ! Potential temperatur at the reference height
      beta = 1.      ! -  (in computing W_*)
      zii  = 1000.    ! m  (pbl height)

      !-----------------------------------------------------------------------
      ! Compute sensible and latent fluxes and their derivatives with respect
      ! to surface temperature using surface temperatures from previous time step.
      !-----------------------------------------------------------------------
      ! Initialization variables
      nmozsgn = 0
      obuold  = 0.
      ! Calculate saturation mixing ratio and change in saturation
      CALL qsadv(tssea,psrf,eg,degdT,qsatg,qsatgdT)

      ! Potential temperatur at the reference height
      thm = tm + 0.0098*ht                    ! intermediate variable equivalent to
                                              ! tm*(pgcm/psrf)**(rgas/cpair)
      th  = tm*(100000./psrf)**(rgas/cpair)   ! potential T
      thv = th*(1.+0.61*qm)                   ! virtual potential T
      ur  = max(0.1,sqrt(us*us+vs*vs))        ! limit set to 0.1

      dth   = thm-tssea                      ! diff of potential temp. between ref. height and surface
      dqh   = qm-qsatg                       ! diff of humidity between ref. height and surface
      dthv  = dth*(1.+0.61*qm)+0.61*th*dqh   ! diff of vir. poten. temp. between ref. height and surface
      !TODO: check if this is correct, inundation may occur over vegetated surface
      zldis = hu-0.                          ! reference height "minus" zero displacement heght

      ! Kinematic viscosity of dry air (m2/s)- Andreas (1989) CRREL Rep. 89-11
      visa=1.326e-5*(1.+6.542e-3*(tm-tfrz) + 8.301e-6*(tm-tfrz)**2 - 4.84e-9*(tm-tfrz)**3)

      ! Loop to obtain initial and good ustar and zo
      ustar=0.06    ! initial value of ustar
      wc=0.5        ! initial value of wc
      !initial value of um
      IF(dthv.ge.0.) THEN
         um=max(ur,0.1)
      ELSE
         um=sqrt(ur*ur+wc*wc)
      ENDIF
      ! initial value of z0mg and ustar
      DO i=1,5
         z0mg=0.013*ustar*ustar/grav+0.11*visa/ustar
         ustar=vonkar*um/log(zldis/z0mg)
      ENDDO
      !
      CALL moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0mg,um,obu)

      ! Evaluated stability-dependent variables using moz from prior iteration
      niters  = 10
      displax = 0.

      !----------------------------------------------------------------
      ITERATION : DO iter = 1, niters         ! begin stability iteration
         !----------------------------------------------------------------
         ! Compute stability-dependent variables
         z0mg = 0.013*ustar*ustar/grav + 0.11*visa/ustar
         xq   = 2.67*(ustar*z0mg/visa)**0.25 - 2.57
         xt   = xq
         z0qg = z0mg/exp(xq)
         z0hg = z0mg/exp(xt)
         !calculation of friction velocity, relation for potential temperature
         IF (DEF_USE_CBL_HEIGHT) THEN
            CALL moninobuk_leddy(hu,ht,hq,displax,z0mg,z0hg,z0qg,obu,um, hpbl, &
                              ustar,fh2m,fq2m,fm10m,fm,fh,fq)
         ELSE
            CALL moninobuk(hu,ht,hq,displax,z0mg,z0hg,z0qg,obu,um,&
                        ustar,fh2m,fq2m,fm10m,fm,fh,fq)
         ENDIF
         !get qstar and tstar
         tstar = vonkar/fh*dth
         qstar = vonkar/fq*dqh

         thvstar=tstar*(1.+0.61*qm)+0.61*th*qstar
         zol=zldis*vonkar*grav*thvstar/(ustar**2*thv) !z/L
         IF(zol >= 0.) THEN       ! stable
            zol = min(2.,max(zol,1.e-6))
         ELSE                     ! unstable
            zol = max(-100.,min(zol,-1.e-6))
         ENDIF
         obu = zldis/zol

         IF(zol >= 0.)THEN
            um = max(ur,0.1) !wind speed at reference height
         ELSE
            IF (DEF_USE_CBL_HEIGHT) THEN !//TODO: Shaofeng, 2023.05.18
               zii = max(5.*hu,hpbl)
            ENDIF !//TODO: Shaofeng, 2023.05.18
            wc = (-grav*ustar*thvstar*zii/thv)**(1./3.) !convective velocity scale
            wc2 = beta*beta*(wc*wc)                     !convective velocity scale squared
            um = sqrt(ur*ur+wc2)                        !wind speed with convective velocity scale
         ENDIF

         IF (obuold*obu < 0.) nmozsgn = nmozsgn+1
         IF(nmozsgn >= 4) EXIT

         obuold = obu

      !----------------------------------------------------------------
      ENDDO ITERATION                         ! end stability iteration
      !----------------------------------------------------------------

      ! Get derivative of fluxes with repect to ground temperature
      ram     =  1./(ustar*ustar/um)
      rah     =  1./(vonkar/fh*ustar)
      raw     =  1./(vonkar/fq*ustar)

      raih    =  rhoair*cpair/rah
      raiw    =  rhoair/raw
      !cgrnds = raih
      !cgrndl = raiw*qsatgdT

      rib = min(5.,zol*ustar**2/(vonkar**2/fh*um**2))

      ! Surface fluxes of momentum, sensible and latent
      ! using ground temperatures from previous time step
      taux   = -rhoair*us/ram
      tauy   = -rhoair*vs/ram

      fseng  = -raih*dth
      fevpg  = -raiw*dqh
      !fsena  = fseng
      !fevpa  = fevpg

      ! 2 m height air temperature
      tref   = thm + vonkar/fh*dth * (fh2m/vonkar - fh/vonkar)
      qref   = qm + vonkar/fq*dqh * (fq2m/vonkar - fq/vonkar)
      z0m   = z0mg
   END SUBROUTINE get_fldevp

#endif
END MODULE MOD_CaMa_colmCaMa
