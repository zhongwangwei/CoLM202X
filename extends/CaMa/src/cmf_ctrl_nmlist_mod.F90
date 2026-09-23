MODULE CMF_CTRL_NMLIST_MOD
!==========================================================
!* PURPOSE: read CaMa-Flood Model configulations from namelist ("input_flood.nam" as default)
!
!* CONTAINS:
! -- CMF_CONFIG_NAMELIST  : read namelist for CaMa-Flood 
! -- CMF_CONFIG_CHECK     : check config conflict
!
! (C) D.Yamazaki & E. Dutra  (U-Tokyo/FCUL)  29Jul 2019
!             Adapted mostly from CMF v362 CONTROL0.F90

  ! Modified by Zhongwang Wei @ SYSU 2022.11.20: add water re-infiltration calculation 

! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
! shared variables in module
USE PARKIND1,                ONLY: JPIM, JPRB, JPRM
USE YOS_CMF_INPUT,           ONLY: LOGNAM
USE YOS_CMF_MAP,             ONLY: REGIONTHIS, REGIONALL

IMPLICIT NONE
CONTAINS
!####################################################################
! -- CMF_CONFIG_NAMELIST  : read namelist for CaMa-Flood 
! -- CMF_CONFIG_CHECK     : check config conflict
!
!####################################################################
SUBROUTINE CMF_CONFIG_NMLIST
USE YOS_CMF_INPUT,      ONLY: TMPNAM,   NSETFILE,   CSETFILE
! run version
USE YOS_CMF_INPUT,      ONLY: LADPSTP,  LFPLAIN,  LKINE,    LFLDOUT,  LPTHOUT,  LDAMOUT,  LDAMIRR,&
                            & LROSPLIT, LGDWDLY,  LSLPMIX,  LMEANSL,  LSEALEV,  LOUTPUT,  &
                            & LRESTART, LSTOONLY, LGRIDMAP, LLEAPYR,  LMAPEND,  LBITSAFE, &
                            & LSTG_ES,  LLEVEE,   LOUTINS,  LOUTINI,  LSEDIMENT,  LTRACE,   &
                            & LSLOPEMOUTH,LWEVAP,LWINFILT, LWEVAPFIX,LWINFILTFIX,LWEXTRACTRIV,       LSPAMAT, LCOLMFEEDBACK
! dimention & time
USE YOS_CMF_INPUT,      ONLY: CDIMINFO, DT,       NX,NY,    NLFP,     NXIN,NYIN,    INPN, &
                            & IFRQ_INP, DTIN,     WEST,EAST,NORTH,SOUTH
! parameters
USE YOS_CMF_INPUT,      ONLY: PMANRIV,  PMANFLD,  PDSTMTH,  PMINSLP,  PGRV, PCADP, &
                            & IMIS, RMIS, DMIS,   CSUFBIN,  CSUFVEC,  CSUFPTH,  CSUFCDF
USE CMF_UTILS_MOD,      ONLY: INQUIRE_FID
IMPLICIT NONE
!* local
CHARACTER(LEN=8)              :: CREG                 !! 
!
NAMELIST/NRUNVER/  LADPSTP,  LFPLAIN,  LKINE,    LFLDOUT,  LPTHOUT,  LDAMOUT,  LDAMIRR,&
                   LROSPLIT, LGDWDLY,  LSLPMIX,  LMEANSL,  LSEALEV,  LOUTPUT,      &
                   LRESTART, LSTOONLY, LGRIDMAP, LLEAPYR,  LMAPEND,  LBITSAFE,     &
                   LSTG_ES,  LLEVEE,   LSEDIMENT,  LTRACE,   LOUTINS,  LSLOPEMOUTH,  &
                   LWEVAP,  LWINFILT,  LWEVAPFIX,LWINFILTFIX,LWEXTRACTRIV,  LOUTINI,  LSPAMAT

NAMELIST/NDIMTIME/ CDIMINFO, DT, IFRQ_INP

NAMELIST/NPARAM/   PMANRIV, PMANFLD, PGRV,    PDSTMTH, PCADP,   PMINSLP, &
                   IMIS, RMIS, DMIS, CSUFBIN, CSUFVEC, CSUFPTH, CSUFCDF
!================================================
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "!--------------------"

! *** 0. SET INPUT UNIT AND OPEN FILE 
NSETFILE=INQUIRE_FID()               !!  for namelist
IF( CSETFILE/="NONE" ) OPEN(NSETFILE,FILE=CSETFILE,STATUS="OLD")
IF( CSETFILE/="NONE" ) WRITE(LOGNAM,*) "CMF::CONFIG_NMLIST: namelist opened: ", TRIM(CSETFILE), NSETFILE

!============================
!*** 1. basic simulation run version

! * defaults
LADPSTP  = .TRUE.            !! true: use adaptive time step
LFPLAIN  = .TRUE.            !! true: consider floodplain (false: only river channel)
LKINE    = .FALSE.           !! true: use kinematic wave
LFLDOUT  = .TRUE.            !! true: floodplain flow (high-water channel flow) active
LPTHOUT  = .FALSE.           !! true: activate bifurcation scheme
LDAMOUT  = .FALSE.           !! true: activate dam operation
LDAMIRR  = .FALSE.           !! true: activate water withdrawals for irrigation  (under development)
LLEVEE   = .FALSE.           !! true: activate levee scheme  (under development)
LSEDIMENT  = .FALSE.           !! true: activate sediment transport (under development)
LTRACE   = .FALSE.           !! true: activate tracer             (under development)
LOUTINS  = .FALSE.           !! true: diagnose instantaneous discharge
LSPAMAT  = .TRUE.            !! true: use quasi sparse matrix (fast but additional memory req)
!!=== this part is used by ECMWF
LROSPLIT = .FALSE.           !! true: input if surface (Qs) and sub-surface (Qsb) runoff
LWEVAP   = .FALSE.           !! true: input evaporation to extract from river 
LWINFILT   = .FALSE.         !! true: input infiltration to extract from river 
LWEVAPFIX= .FALSE.           !! deprecated compatibility option; no effect
LWINFILTFIX= .FALSE.         !! deprecated compatibility option; no effect
LGDWDLY  = .FALSE.           !! true: Activate ground water reservoir and delay
LSLPMIX  = .FALSE.           !! true: activate mixed kinematic and local inertia based on slope
LWEXTRACTRIV=.FALSE.         !! true: also extract water from rivers 
LSLOPEMOUTH =.FALSE.         !! true: prescribe water level slope == elevation slope on river month
!!===

!! dinamic sea level
LMEANSL  = .FALSE.           !! true : boundary condition for mean sea level
LSEALEV  = .FALSE.           !! true : boundary condition for variable sea level

!! restaer & output
LRESTART = .FALSE.           !! true: initial condition from restart file
LSTOONLY = .FALSE.           !! true: storage only restart (mainly for data assimilation)
LOUTPUT  = .FALSE.           !! true: use standard output (to file)
LOUTINI  = .FALSE.           !! true: output initial storage (netCDF only)

LGRIDMAP = .TRUE.            !! true: for standard XY gridded 2D map
LLEAPYR  = .FALSE.           !! true: neglect leap year (Feb29 skipped)
LMAPEND  = .FALSE.           !! true: for map data endian conversion
LBITSAFE = .FALSE.           !! true: for Bit Identical (not used from v410, set in Mkinclude)
LSTG_ES  = .FALSE.           !! true: for Vector Processor optimization (CMF_OPT_FLDSTG_ES) 

!* change
IF( CSETFILE/="NONE" ) REWIND(NSETFILE)
IF( CSETFILE/="NONE" ) READ(NSETFILE,NML=NRUNVER)
IF( CSETFILE=="NONE" )THEN
  ! CoLM NC coupling: fixed local-inertial floodplain solver and two-way sinks.
  LWEVAP=LCOLMFEEDBACK      ! two-way exchange, unless CoLM switched it off (DEF_CaMa_FloodFeedback)
  LWINFILT=LCOLMFEEDBACK
  LOUTPUT=.TRUE.       ! select diagnostics; CoLM owns history file writing
ENDIF

WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "=== NAMELIST, NRUNVER ==="
WRITE(LOGNAM,*) "LADPSTP ",  LADPSTP
WRITE(LOGNAM,*) "LFPLAIN ",  LFPLAIN
WRITE(LOGNAM,*) "LKINE   ",  LKINE
WRITE(LOGNAM,*) "LFLDOUT ",  LFLDOUT
WRITE(LOGNAM,*) "LPTHOUT ",  LPTHOUT
WRITE(LOGNAM,*) "LDAMOUT ",  LDAMOUT
WRITE(LOGNAM,*) "LDAMIRR ",  LDAMIRR
WRITE(LOGNAM,*) "LLEVEE  ",  LLEVEE
WRITE(LOGNAM,*) "LSEDIMENT ",  LSEDIMENT
WRITE(LOGNAM,*) "LTRACE  ",  LTRACE
WRITE(LOGNAM,*) "LOUTINS ",  LOUTINS
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "LROSPLIT ", LROSPLIT
WRITE(LOGNAM,*) "LWEVAP   ", LWEVAP
WRITE(LOGNAM,*) "LWEVAPFIX", LWEVAPFIX
WRITE(LOGNAM,*) "LWINFILT   ", LWINFILT
WRITE(LOGNAM,*) "LWINFILTFIX", LWINFILTFIX
WRITE(LOGNAM,*) "LWEXTRACTRIV", LWEXTRACTRIV
WRITE(LOGNAM,*) "LGDWDLY  ", LGDWDLY
WRITE(LOGNAM,*) "LSLPMIX  ", LSLPMIX
WRITE(LOGNAM,*) "LSLOPEMOUTH ", LSLOPEMOUTH
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "LMEANSL: ", LSEALEV
WRITE(LOGNAM,*) "LSEALEV: ", LSEALEV
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "LRESTART ", LRESTART
WRITE(LOGNAM,*) "LSTOONLY ", LSTOONLY
WRITE(LOGNAM,*) "LOUTPUT  ", LOUTPUT
WRITE(LOGNAM,*) "LOUTINI  ", LOUTINI
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "LGRIDMAP ", LGRIDMAP
WRITE(LOGNAM,*) "LLEAPYR  ", LLEAPYR
WRITE(LOGNAM,*) "LMAPEND  ", LMAPEND
WRITE(LOGNAM,*) "LBITSAFE ", LBITSAFE
WRITE(LOGNAM,*) "LSTG_ES " , LSTG_ES
WRITE(LOGNAM,*) "LSPAMAT " , LSPAMAT

!============================
!*** 2. set model dimention & time

!* defaults (from namelist)
CDIMINFO ="NONE"
DT       = 60*60             !! hourly outer step (adaptive substeps may be shorter)
IFRQ_INP = 1                 !! hourly input/coupling

!* change
IF( CSETFILE/="NONE" ) REWIND(NSETFILE)
IF( CSETFILE/="NONE" ) READ(NSETFILE,NML=NDIMTIME)

DTIN  = IFRQ_INP*60*60       !! hour -> second

WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "=== NAMELIST, NCONF ==="
WRITE(LOGNAM,*) "CDIMINFO  ", TRIM(CDIMINFO)
WRITE(LOGNAM,*) "DT        ", DT
WRITE(LOGNAM,*) "DTIN      ", DTIN
WRITE(LOGNAM,*) "IFRQ_INP  ", IFRQ_INP

!==========
!* default (from diminfo)
NX    = 1440                 !! 15 minute resolution
NY    = 720
NLFP  = 10                   !! 10 floodplain layer
NXIN  = 360                  !! 1 degree input
NYIN  = 180
INPN  = 1                    !! maximum number of input grids corresponding to one CaMa-Flood grid
WEST  = -180._JPRB              !! west, east, north, south edges of the domain
EAST  =  180._JPRB
NORTH =  90._JPRB
SOUTH = -90._JPRB

! Dimension-file I/O is deferred to CMF_MAPS_NMLIST, after routing NC selection.

!============================
!*** 3. set PARAM: parameters
! * defaults
PMANRIV=0.03_JPRB                              !! manning coefficient river
PMANFLD=0.10_JPRB                              !! manning coefficient floodplain
PGRV   =9.8_JPRB                               !! gravity accerelation
PDSTMTH=10000._JPRB                            !! downstream distance at river mouth [m]
PCADP  =0.7_JPRB                               !! CFL coefficient
PMINSLP=1.E-5                                  !! minimum slope (kinematic wave)

IMIS=-9999_JPIM
RMIS=1.E20_JPRM
DMIS=1.E20_JPRB

CSUFBIN='.bin'
CSUFVEC='.vec'
CSUFPTH='.pth'
CSUFCDF='.nc'

! * change
IF( CSETFILE/="NONE" ) REWIND(NSETFILE)
IF( CSETFILE/="NONE" ) READ(NSETFILE,NML=NPARAM)

WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "=== NAMELIST, NPARAM ==="
WRITE(LOGNAM,*) "PMANRIV  ", PMANRIV
WRITE(LOGNAM,*) "PMANFLD  ", PMANFLD
WRITE(LOGNAM,*) "PGRV     ", PGRV
WRITE(LOGNAM,*) "PDSTMTH  ", PDSTMTH
WRITE(LOGNAM,*) "PCADP    ", PCADP
WRITE(LOGNAM,*) "PMINSLP  ", PMINSLP
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "IMIS     ", IMIS
WRITE(LOGNAM,*) "RMIS     ", RMIS
WRITE(LOGNAM,*) "DMIS     ", DMIS
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "CSUFBIN  ", TRIM(CSUFBIN)
WRITE(LOGNAM,*) "CSUFVEC  ", TRIM(CSUFVEC)
WRITE(LOGNAM,*) "CSUFPTH  ", TRIM(CSUFPTH)
WRITE(LOGNAM,*) "CSUFCDF  ", TRIM(CSUFCDF)

!===============================
!*** CLOSE FILE 
IF( CSETFILE/="NONE" ) CLOSE(NSETFILE)

WRITE(LOGNAM,*) "CMF::CONFIG_NMLIST: end "

WRITE(LOGNAM,*) "--------------------!"
WRITE(LOGNAM,*) ""

IF (REGIONALL>=2 )then 
  WRITE(CREG,'(I0)') REGIONTHIS                !! Regional Output for MPI run
  CSUFVEC=TRIM(CSUFVEC)//'-'//TRIM(CREG)       !! Change suffix of output file for each MPI node (only vector output)
ENDIF


END SUBROUTINE CMF_CONFIG_NMLIST
!####################################################################




!####################################################################
SUBROUTINE CMF_CONFIG_CHECK
USE YOS_CMF_INPUT,      ONLY: LADPSTP,  LFPLAIN,  LKINE,    LPTHOUT,     &
                            & LROSPLIT, LGDWDLY,  LSEALEV, &
                            & LWEVAP,   LWEVAPFIX,LWINFILT,   LWINFILTFIX, LWEXTRACTRIV
USE YOS_CMF_INPUT,      ONLY: DT, DTIN, DTSL
IMPLICIT NONE
!================================================

WRITE(LOGNAM,*) "CMF::CONFIG_CHECK: check setting conflicts"

!*** 1. check for time step
IF ( DT<60 .or. MOD( INT(DT),60 )/=0 ) THEN
  WRITE(LOGNAM,*) "DT= ", DT
  WRITE(LOGNAM,*) "DT should be multiple of 60. CaMa-Flood controls time by MINUTE"
  WRITE(LOGNAM,*) "stop"
  STOP 9
ENDIF

IF ( MOD( INT(DTIN), INT(DT) )/=0 ) THEN
  WRITE(LOGNAM,*) "DTIN, DT= ", DTIN, DT
  WRITE(LOGNAM,*) "DTIN should be multiple of DT"
  WRITE(LOGNAM,*) "stop"
  STOP 9
ENDIF

IF ( LSEALEV .and. MOD( INT(DTSL), INT(DT) )/=0 ) THEN
  WRITE(LOGNAM,*) "DTSL, DT= ", DTIN, DT
  WRITE(LOGNAM,*) "DTSL should be multiple of DT"
  WRITE(LOGNAM,*) "stop"
  STOP 9
ENDIF

!*** 2. check for physics options

IF ( .not.LFPLAIN .AND. .not.LKINE ) THEN
  WRITE(LOGNAM,*) "LFPLAIN=.false. & LKINE=.false."
  WRITE(LOGNAM,*) "CAUTION: NO FLOODPLAIN OPTION reccomended to be used with kinematic wave (LKINE=.true.)"
ENDIF

IF ( LKINE .AND. LADPSTP ) THEN
  WRITE(LOGNAM,*) "LKINE=.true. & LADPSTP=.true."
  WRITE(LOGNAM,*) "adaptive time step reccoomended only with local inertial equation (LKINE=.false.)"
  WRITE(LOGNAM,*) "Set appropriate fixed time step for Kinematic Wave"
ENDIF

IF ( LKINE .AND. LPTHOUT ) THEN
  WRITE(LOGNAM,*) "LKINE=.true. & LPATHOUT=.true."
  WRITE(LOGNAM,*) "bifurcation channel flow only available with local inertial equation (LKINE=.false.)"
  WRITE(LOGNAM,*) "STOP"
  STOP 9
ENDIF

IF ( LGDWDLY .AND. .NOT. LROSPLIT ) THEN
  WRITE(LOGNAM,*) "LGDWDLY=true and LROSPLIT=false"
  WRITE(LOGNAM,*) "Ground water reservoir can only be active when runoff splitting is on"
ENDIF

IF ( LWEVAPFIX .OR. LWINFILTFIX ) THEN
  WRITE(LOGNAM,*) 'WARNING: LWEVAPFIX/LWINFILTFIX are deprecated and have no effect.'
  WRITE(LOGNAM,*) 'CoLM floodwater credit/debit accounting is always applied to enabled sinks.'
ENDIF
IF ( LWEXTRACTRIV .AND. .NOT. LWEVAP ) THEN
  WRITE(LOGNAM,*) "LWEXTRACTRIV=true and LWEVAP=false"
  WRITE(LOGNAM,*) "LWEXTRACTRIV can only be active if LWEVAP is active"
ENDIF 


WRITE(LOGNAM,*) "CMF::CONFIG_CHECK: end"

END SUBROUTINE CMF_CONFIG_CHECK
!####################################################################

END MODULE CMF_CTRL_NMLIST_MOD
