#include <define.h>

#ifdef ConstituentTransport
MODULE MOD_CT_Namelist
!-----------------------------------------------------------------------
! !DESCRIPTION:
!
!  Switches and settings of constituent transport.
!
!  They are kept out of MOD_Namelist on purpose: the main namelist file
!  carries one extra group, read here on the master and broadcast,
!
!     &nl_colm_constituent
!        DEF_USE_ISOTOPE        = .true.
!        DEF_ISOTOPE_SPECIES    = 'H2_18O,HDO'
!        DEF_ISOTOPE_PARAM_FILE = 'run/standard_isotope_parameter.nml'
!        DEF_SOLUTE_NAMES       = 'CL'
!        DEF_SOLUTE_PARAM_FILES = 'CL:run/standard_chloride_parameter.nml'
!     /
!
!  The group may be left out; every setting then keeps its default and no
!  constituent is carried. DEF_USE_SEDIMENT stays in &nl_colm, where
!  upstream has it.
!
!  Isotope physics settings live in &nl_colm_isotope_parameter of
!  DEF_ISOTOPE_PARAM_FILE ('null' = code defaults). Solute parameter files
!  are read per species by MOD_CT_Defs; per-species forcing by
!  MOD_CT_ForcingInput.
!-----------------------------------------------------------------------

   USE MOD_Precision
   IMPLICIT NONE
   SAVE

   ! ----- &nl_colm_constituent -----
   logical            :: DEF_USE_ISOTOPE        = .false.
   ! Only species with compiled physics: H2_18O, HDO.
   character(len=256) :: DEF_ISOTOPE_SPECIES    = 'H2_18O,HDO'
   character(len=256) :: DEF_ISOTOPE_PARAM_FILE = 'null'
   character(len=256) :: DEF_SOLUTE_NAMES       = ''
   ! 'NAME:path' entries separated by ',' or ';', or plain paths in the
   ! order of DEF_SOLUTE_NAMES. A solute without a file uses code defaults.
   character(len=2048):: DEF_SOLUTE_PARAM_FILES = 'null'

   ! Abort after this many cells fail a check (0 = report only).
   integer            :: DEF_CT_BALANCE_ABORT_NBAD = 0
   integer            :: DEF_CT_RESID_ABORT_NBAD   = 0
   integer            :: DEF_CT_LULCC_ABORT_NBAD   = 0
   ! Water depth that carries the aquifer signature when the
   ! variably saturated flow scheme is on (required > 0 for isotopes then).
   real(r8)           :: DEF_CT_AQUIFER_MIXING_WATER_MM = -1._r8
   logical            :: DEF_CT_USE_SOIL_INIT   = .false.
   character(len=256) :: DEF_CT_SOIL_INIT_FILE  = 'null'
   character(len=256) :: DEF_CT_SOIL_INIT_VARS  = ''

   ! ----- &nl_colm_isotope_parameter (in DEF_ISOTOPE_PARAM_FILE) -----
   logical            :: DEF_ISO_USE_FRACTIONATION     = .false.
   character(len=16)  :: DEF_ISO_KINETIC_SCHEME        = 'MERLIVAT1978'
   real(r8)           :: DEF_ISO_ICE_SUPERSAT_SLOPE    = 0.003_r8
   real(r8)           :: DEF_ISO_CG_RELHUM_MAX         = 0.99_r8
   character(len=16)  :: DEF_ISO_OPEN_WATER_KINETIC    = 'MJ79'
   real(r8)           :: DEF_ISO_SUBL_SKIN_MM          = 5.0_r8
   character(len=16)  :: DEF_ISO_SOIL_KINETIC          = 'RESISTANCE'
   logical            :: DEF_ISO_SOIL_DIFFUSION        = .true.
   logical            :: DEF_ISO_SOIL_VAPOR_DIFFUSION  = .true.
   real(r8)           :: DEF_ISO_CANOPY_EQUILIBRATION  = 0.0_r8
   real(r8)           :: DEF_ISO_SNOWMELT_EQUILIBRATION= 0.0_r8
   real(r8)           :: DEF_ISO_NSS_LEAF_WATER_PER_LAI= 0.12_r8
   real(r8)           :: DEF_ISO_NSS_LEAF_PATH_LENGTH  = 0.01_r8
   real(r8)           :: DEF_ISO_NSS_LEAF_RB           = 100._r8
   ! Per-species values; defaults are VSMOW and the former
   ! standard_O18/HDO_parameter.nml.
   real(r8)           :: DEF_ISO_REF_RATIO_H2_18O      = 2.0052e-3_r8
   real(r8)           :: DEF_ISO_INIT_DELTA_H2_18O     = -10._r8
   real(r8)           :: DEF_ISO_REF_RATIO_HDO         = 1.5576e-4_r8
   real(r8)           :: DEF_ISO_INIT_DELTA_HDO        = -70._r8

   logical, private :: loaded = .false.

   PUBLIC :: ct_read_namelist
   PUBLIC :: ct_any_enabled
   PUBLIC :: ct_namelist_group_present

CONTAINS

   SUBROUTINE ct_read_namelist (nlfile)
   !-----------------------------------------------------------------------
   !  Read the constituent groups once. nlfile is the main namelist file;
   !  without it the first command-line argument is used, as in CoLM.F90.
   !-----------------------------------------------------------------------
   USE MOD_SPMD_Task
   IMPLICIT NONE

   character(len=*), intent(in), optional :: nlfile

   character(len=256) :: mainfile
   integer :: ierr, unit_nml
   character(len=512) :: iomsg

      namelist /nl_colm_constituent/ DEF_USE_ISOTOPE, DEF_ISOTOPE_SPECIES, &
         DEF_ISOTOPE_PARAM_FILE, DEF_SOLUTE_NAMES, DEF_SOLUTE_PARAM_FILES, &
         DEF_CT_BALANCE_ABORT_NBAD, DEF_CT_RESID_ABORT_NBAD, DEF_CT_LULCC_ABORT_NBAD, &
         DEF_CT_AQUIFER_MIXING_WATER_MM, &
         DEF_CT_USE_SOIL_INIT, DEF_CT_SOIL_INIT_FILE, DEF_CT_SOIL_INIT_VARS

      namelist /nl_colm_isotope_parameter/ DEF_ISO_USE_FRACTIONATION, &
         DEF_ISO_KINETIC_SCHEME, DEF_ISO_ICE_SUPERSAT_SLOPE, DEF_ISO_CG_RELHUM_MAX, &
         DEF_ISO_OPEN_WATER_KINETIC, DEF_ISO_SUBL_SKIN_MM, DEF_ISO_SOIL_KINETIC, &
         DEF_ISO_SOIL_DIFFUSION, DEF_ISO_SOIL_VAPOR_DIFFUSION, &
         DEF_ISO_CANOPY_EQUILIBRATION, DEF_ISO_SNOWMELT_EQUILIBRATION, &
         DEF_ISO_NSS_LEAF_WATER_PER_LAI, DEF_ISO_NSS_LEAF_PATH_LENGTH, DEF_ISO_NSS_LEAF_RB, &
         DEF_ISO_REF_RATIO_H2_18O, DEF_ISO_INIT_DELTA_H2_18O, &
         DEF_ISO_REF_RATIO_HDO, DEF_ISO_INIT_DELTA_HDO

      IF (loaded) RETURN

      IF (p_is_master) THEN

         IF (present(nlfile)) THEN
            mainfile = nlfile
         ELSE
            CALL get_command_argument (1, mainfile)
         ENDIF

         IF (ct_namelist_group_present (mainfile, 'nl_colm_constituent')) THEN
            open (newunit=unit_nml, status='OLD', file=trim(mainfile), form='FORMATTED')
            iomsg = ''
            read (unit_nml, nml=nl_colm_constituent, iostat=ierr, iomsg=iomsg)
            close (unit_nml)
            IF (ierr /= 0) THEN
               write(*,*) 'ERROR in ', trim(mainfile), ' : ', trim(iomsg)
               CALL CoLM_Stop (' ***** ERROR: Problem reading &nl_colm_constituent')
            ENDIF
         ENDIF

         IF (len_trim(DEF_SOLUTE_PARAM_FILES) == len(DEF_SOLUTE_PARAM_FILES)) &
            CALL CoLM_Stop ('DEF_SOLUTE_PARAM_FILES is too long.')

         IF (DEF_USE_ISOTOPE .and. ct_file_given (DEF_ISOTOPE_PARAM_FILE)) THEN
            IF (.not. ct_namelist_group_present (DEF_ISOTOPE_PARAM_FILE, &
                                                 'nl_colm_isotope_parameter')) THEN
               CALL CoLM_Stop ('No &nl_colm_isotope_parameter in '//trim(DEF_ISOTOPE_PARAM_FILE))
            ENDIF
            open (newunit=unit_nml, status='OLD', file=trim(DEF_ISOTOPE_PARAM_FILE), form='FORMATTED')
            iomsg = ''
            read (unit_nml, nml=nl_colm_isotope_parameter, iostat=ierr, iomsg=iomsg)
            close (unit_nml)
            IF (ierr /= 0) THEN
               write(*,*) 'ERROR in ', trim(DEF_ISOTOPE_PARAM_FILE), ' : ', trim(iomsg)
               CALL CoLM_Stop (' ***** ERROR: Problem reading &nl_colm_isotope_parameter')
            ENDIF
         ENDIF

         CALL check_settings ()

      ENDIF

#ifdef USEMPI
      CALL mpi_bcast (DEF_USE_ISOTOPE               ,1   ,mpi_logical   ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISOTOPE_SPECIES           ,256 ,mpi_character ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISOTOPE_PARAM_FILE        ,256 ,mpi_character ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_SOLUTE_NAMES              ,256 ,mpi_character ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_SOLUTE_PARAM_FILES        ,len(DEF_SOLUTE_PARAM_FILES) ,mpi_character ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_CT_BALANCE_ABORT_NBAD     ,1   ,mpi_integer   ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_CT_RESID_ABORT_NBAD       ,1   ,mpi_integer   ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_CT_LULCC_ABORT_NBAD       ,1   ,mpi_integer   ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_CT_AQUIFER_MIXING_WATER_MM,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_CT_USE_SOIL_INIT          ,1   ,mpi_logical   ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_CT_SOIL_INIT_FILE         ,256 ,mpi_character ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_CT_SOIL_INIT_VARS         ,256 ,mpi_character ,p_address_master ,p_comm_glb ,p_err)

      CALL mpi_bcast (DEF_ISO_USE_FRACTIONATION     ,1   ,mpi_logical   ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_KINETIC_SCHEME        ,16  ,mpi_character ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_ICE_SUPERSAT_SLOPE    ,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_CG_RELHUM_MAX         ,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_OPEN_WATER_KINETIC    ,16  ,mpi_character ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_SUBL_SKIN_MM          ,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_SOIL_KINETIC          ,16  ,mpi_character ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_SOIL_DIFFUSION        ,1   ,mpi_logical   ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_SOIL_VAPOR_DIFFUSION  ,1   ,mpi_logical   ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_CANOPY_EQUILIBRATION  ,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_SNOWMELT_EQUILIBRATION,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_NSS_LEAF_WATER_PER_LAI,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_NSS_LEAF_PATH_LENGTH  ,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_NSS_LEAF_RB           ,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_REF_RATIO_H2_18O      ,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_INIT_DELTA_H2_18O     ,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_REF_RATIO_HDO         ,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
      CALL mpi_bcast (DEF_ISO_INIT_DELTA_HDO        ,1   ,mpi_real8     ,p_address_master ,p_comm_glb ,p_err)
#endif

      loaded = .true.

   END SUBROUTINE ct_read_namelist


   SUBROUTINE check_settings ()

   USE MOD_SPMD_Task, only: CoLM_Stop
   IMPLICIT NONE

      SELECT CASE (trim(adjustl(DEF_ISO_KINETIC_SCHEME)))
      CASE ('CAPPA2003', 'Cappa2003', 'cappa2003')
         DEF_ISO_KINETIC_SCHEME = 'CAPPA2003'
      CASE ('MERLIVAT1978', 'Merlivat1978', 'merlivat1978')
         DEF_ISO_KINETIC_SCHEME = 'MERLIVAT1978'
      CASE DEFAULT
         CALL CoLM_Stop ('Invalid DEF_ISO_KINETIC_SCHEME.')
      END SELECT

      SELECT CASE (trim(adjustl(DEF_ISO_OPEN_WATER_KINETIC)))
      CASE ('MJ79', 'mj79')
         DEF_ISO_OPEN_WATER_KINETIC = 'MJ79'
      CASE ('EXPONENT', 'exponent')
         DEF_ISO_OPEN_WATER_KINETIC = 'EXPONENT'
      CASE DEFAULT
         CALL CoLM_Stop ('Invalid DEF_ISO_OPEN_WATER_KINETIC.')
      END SELECT

      SELECT CASE (trim(adjustl(DEF_ISO_SOIL_KINETIC)))
      CASE ('RESISTANCE', 'resistance')
         DEF_ISO_SOIL_KINETIC = 'RESISTANCE'
      CASE ('EXPONENT', 'exponent')
         DEF_ISO_SOIL_KINETIC = 'EXPONENT'
      CASE DEFAULT
         CALL CoLM_Stop ('Invalid DEF_ISO_SOIL_KINETIC.')
      END SELECT

      IF (DEF_ISO_CG_RELHUM_MAX <= 0._r8 .or. DEF_ISO_CG_RELHUM_MAX >= 1._r8 .or. &
          DEF_ISO_CANOPY_EQUILIBRATION < 0._r8 .or. DEF_ISO_CANOPY_EQUILIBRATION > 1._r8 .or. &
          DEF_ISO_SNOWMELT_EQUILIBRATION < 0._r8 .or. DEF_ISO_SNOWMELT_EQUILIBRATION > 1._r8 .or. &
          DEF_ISO_SUBL_SKIN_MM < 0._r8 .or. DEF_ISO_ICE_SUPERSAT_SLOPE < 0._r8) &
         CALL CoLM_Stop ('Invalid isotope fractionation parameter.')

      IF (DEF_CT_BALANCE_ABORT_NBAD < 0 .or. DEF_CT_RESID_ABORT_NBAD < 0 .or. &
          DEF_CT_LULCC_ABORT_NBAD < 0) CALL CoLM_Stop ('Invalid DEF_CT_*_ABORT_NBAD.')

   END SUBROUTINE check_settings


   logical FUNCTION ct_any_enabled ()
   ! True when at least one constituent is switched on.
   IMPLICIT NONE
      ct_any_enabled = DEF_USE_ISOTOPE .or. len_trim(DEF_SOLUTE_NAMES) > 0
   END FUNCTION ct_any_enabled


   logical FUNCTION ct_file_given (fname)
   IMPLICIT NONE
   character(len=*), intent(in) :: fname
      ct_file_given = len_trim(fname) > 0
      IF (ct_file_given) ct_file_given = trim(adjustl(fname)) /= 'null' &
                                   .and. trim(adjustl(fname)) /= 'NULL'
   END FUNCTION ct_file_given


   logical FUNCTION ct_namelist_group_present (fname, group)
   !-----------------------------------------------------------------------
   !  Whether a file holds "&group" at the start of a line, so that an
   !  optional group can be skipped rather than read to end of file.
   !-----------------------------------------------------------------------
   IMPLICIT NONE
   character(len=*), intent(in) :: fname, group

   character(len=1024) :: line
   character(len=128)  :: key
   logical :: fexists
   integer :: ierr, unit_nml, i, ia

      ct_namelist_group_present = .false.
      INQUIRE (file=trim(fname), exist=fexists)
      IF (.not. fexists) THEN
         IF (trim(group) == 'nl_colm_constituent') RETURN
         CALL stop_missing (fname)
      ENDIF

      key = '&'//trim(group)
      open (newunit=unit_nml, status='OLD', file=trim(fname), form='FORMATTED')
      DO
         read (unit_nml, '(A)', iostat=ierr) line
         IF (ierr /= 0) EXIT
         line = adjustl(line)
         DO i = 1, len_trim(line)
            ia = iachar(line(i:i))
            IF (ia >= iachar('A') .and. ia <= iachar('Z')) line(i:i) = achar(ia+32)
         ENDDO
         IF (line(1:1) == '$') line(1:1) = '&'
         IF (index(line, trim(key)) == 1) THEN
            i = len_trim(key) + 1
            IF (i > len_trim(line) .or. line(i:i) == ' ' .or. line(i:i) == '!') THEN
               ct_namelist_group_present = .true.
               EXIT
            ENDIF
         ENDIF
      ENDDO
      close (unit_nml)

   CONTAINS

      SUBROUTINE stop_missing (f)
      USE MOD_SPMD_Task, only: CoLM_Stop
      character(len=*), intent(in) :: f
         CALL CoLM_Stop ('Constituent parameter file not found: '//trim(f))
      END SUBROUTINE stop_missing

   END FUNCTION ct_namelist_group_present

END MODULE MOD_CT_Namelist
#endif
