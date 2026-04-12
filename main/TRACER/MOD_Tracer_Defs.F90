#include <define.h>

MODULE MOD_Tracer_Defs

   USE MOD_Precision
   USE MOD_Namelist, only: DEF_TRACER_NUM, DEF_TRACER_NAMES, DEF_TRACER_TYPES, &
      DEF_TRACER_MRAT, DEF_TRACER_REF_RATIO, DEF_TRACER_INIT_DELTA

   IMPLICIT NONE
   SAVE

   type :: tracer_info_type
      character(len=32) :: name
      character(len=16) :: category
      real(r8)          :: mol_weight
      real(r8)          :: ref_ratio
      real(r8)          :: init_delta
      logical           :: has_fractionation
   end type tracer_info_type

   integer :: ntracers = 0
   type(tracer_info_type), allocatable :: tracers(:)

   real(r8), parameter :: Rsmow_18O = 2.0052e-3_r8
   real(r8), parameter :: Rsmow_D   = 1.5576e-4_r8
   real(r8), parameter :: trc_tiny = 1.0e-30_r8

   PUBLIC :: tracer_defs_init, tracer_defs_final
   PUBLIC :: mass_to_delta, delta_to_R, R_to_mass
   PUBLIC :: ntracers, tracers
   PUBLIC :: tracer_info_type
   PUBLIC :: Rsmow_18O, Rsmow_D, trc_tiny

CONTAINS

   SUBROUTINE tracer_defs_init ()
      IMPLICIT NONE
      integer :: i
      character(len=32) :: tokens(100)
      integer :: ntokens

      ntracers = DEF_TRACER_NUM
      allocate(tracers(ntracers))

      CALL parse_csv(DEF_TRACER_NAMES, tokens, ntokens)
      DO i = 1, min(ntracers, ntokens)
         tracers(i)%name = ADJUSTL(TRIM(tokens(i)))
      ENDDO

      CALL parse_csv(DEF_TRACER_TYPES, tokens, ntokens)
      DO i = 1, min(ntracers, ntokens)
         tracers(i)%category = ADJUSTL(TRIM(tokens(i)))
      ENDDO
      DO i = ntokens+1, ntracers
         tracers(i)%category = 'isotope'
      ENDDO

      CALL parse_csv_real(DEF_TRACER_MRAT, ntracers, tracers(:)%mol_weight, 18.0_r8)
      CALL parse_csv_real(DEF_TRACER_REF_RATIO, ntracers, tracers(:)%ref_ratio, 1.0_r8)
      CALL parse_csv_real(DEF_TRACER_INIT_DELTA, ntracers, tracers(:)%init_delta, 0.0_r8)

      DO i = 1, ntracers
         tracers(i)%has_fractionation = (tracers(i)%category == 'isotope')
      ENDDO
   END SUBROUTINE tracer_defs_init

   SUBROUTINE tracer_defs_final ()
      IF (allocated(tracers)) deallocate(tracers)
      ntracers = 0
   END SUBROUTINE tracer_defs_final

   pure real(r8) FUNCTION mass_to_delta (trc_mass, water_mass, ref_ratio)
      real(r8), intent(in) :: trc_mass, water_mass, ref_ratio
      real(r8) :: R_sample
      IF (water_mass > trc_tiny .and. ref_ratio > trc_tiny) THEN
         R_sample = trc_mass / water_mass
         mass_to_delta = (R_sample / ref_ratio - 1.0_r8) * 1000.0_r8
      ELSE
         mass_to_delta = 0.0_r8
      ENDIF
   END FUNCTION mass_to_delta

   pure real(r8) FUNCTION delta_to_R (delta, ref_ratio)
      real(r8), intent(in) :: delta, ref_ratio
      delta_to_R = ref_ratio * (1.0_r8 + delta / 1000.0_r8)
   END FUNCTION delta_to_R

   pure real(r8) FUNCTION R_to_mass (delta, water_mass, ref_ratio)
      real(r8), intent(in) :: delta, water_mass, ref_ratio
      R_to_mass = water_mass * delta_to_R(delta, ref_ratio)
   END FUNCTION R_to_mass

   SUBROUTINE parse_csv (csvstr, tokens, ntokens)
      character(len=*), intent(in)  :: csvstr
      character(len=32), intent(out) :: tokens(:)
      integer, intent(out) :: ntokens
      integer :: i, j, slen
      character(len=256) :: buf
      buf = ADJUSTL(TRIM(csvstr))
      slen = LEN_TRIM(buf)
      ntokens = 0; j = 1
      DO i = 1, slen
         IF (buf(i:i) == ',') THEN
            ntokens = ntokens + 1
            tokens(ntokens) = ADJUSTL(TRIM(buf(j:i-1)))
            j = i + 1
         ENDIF
      ENDDO
      IF (j <= slen) THEN
         ntokens = ntokens + 1
         tokens(ntokens) = ADJUSTL(TRIM(buf(j:slen)))
      ENDIF
   END SUBROUTINE parse_csv

   SUBROUTINE parse_csv_real (csvstr, n, vals, default_val)
      character(len=*), intent(in)  :: csvstr
      integer, intent(in) :: n
      real(r8), intent(out) :: vals(n)
      real(r8), intent(in) :: default_val
      character(len=32) :: tokens(100)
      integer :: ntokens, i, ios
      CALL parse_csv(csvstr, tokens, ntokens)
      DO i = 1, n
         IF (i <= ntokens) THEN
            READ(tokens(i), *, iostat=ios) vals(i)
            IF (ios /= 0) vals(i) = default_val
         ELSE
            vals(i) = default_val
         ENDIF
      ENDDO
   END SUBROUTINE parse_csv_real

END MODULE MOD_Tracer_Defs
