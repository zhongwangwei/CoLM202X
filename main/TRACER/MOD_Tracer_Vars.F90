#include <define.h>

MODULE MOD_Tracer_Vars

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers

   IMPLICIT NONE
   SAVE

   real(r8), allocatable :: trc_ldew_rain  (:,:)
   real(r8), allocatable :: trc_ldew_snow  (:,:)
   real(r8), allocatable :: trc_wliq_soisno(:,:,:)
   real(r8), allocatable :: trc_wice_soisno(:,:,:)
   real(r8), allocatable :: trc_wa         (:,:)
   real(r8), allocatable :: trc_wdsrf      (:,:)
   real(r8), allocatable :: trc_wetwat     (:,:)

   real(r8), allocatable :: a_trc_precip   (:,:)
   real(r8), allocatable :: a_trc_evap     (:,:)
   real(r8), allocatable :: a_trc_trans    (:,:)
   real(r8), allocatable :: a_trc_rsur     (:,:)
   real(r8), allocatable :: a_trc_rsub     (:,:)
   real(r8), allocatable :: a_trc_rnof     (:,:)
   real(r8), allocatable :: a_trc_qinfl    (:,:)
   real(r8), allocatable :: a_trc_qcharge  (:,:)

   real(r8), allocatable :: trc_storage_beg(:,:)
   real(r8), allocatable :: trc_balance_err(:,:)

   integer :: trc_hist_nac = 0
   real(r8), allocatable :: a_trc_ldew_mass (:,:)
   real(r8), allocatable :: a_water_ldew    (:)
   real(r8), allocatable :: a_trc_soil_mass (:,:,:)
   real(r8), allocatable :: a_water_soil    (:,:)

   PUBLIC :: allocate_Tracer_Vars, deallocate_Tracer_Vars, flush_Tracer_Acc

CONTAINS

   SUBROUTINE allocate_Tracer_Vars (numpatch, maxsnl, nl_soil)
      IMPLICIT NONE
      integer, intent(in) :: numpatch, maxsnl, nl_soil
      IF (ntracers <= 0 .or. numpatch <= 0) RETURN

      allocate(trc_ldew_rain   (ntracers, numpatch));           trc_ldew_rain   = 0._r8
      allocate(trc_ldew_snow   (ntracers, numpatch));           trc_ldew_snow   = 0._r8
      allocate(trc_wliq_soisno (ntracers, maxsnl+1:nl_soil, numpatch)); trc_wliq_soisno = 0._r8
      allocate(trc_wice_soisno (ntracers, maxsnl+1:nl_soil, numpatch)); trc_wice_soisno = 0._r8
      allocate(trc_wa          (ntracers, numpatch));           trc_wa          = 0._r8
      allocate(trc_wdsrf       (ntracers, numpatch));           trc_wdsrf       = 0._r8
      allocate(trc_wetwat      (ntracers, numpatch));           trc_wetwat      = 0._r8

      allocate(a_trc_precip    (ntracers, numpatch));           a_trc_precip    = 0._r8
      allocate(a_trc_evap      (ntracers, numpatch));           a_trc_evap      = 0._r8
      allocate(a_trc_trans     (ntracers, numpatch));           a_trc_trans     = 0._r8
      allocate(a_trc_rsur      (ntracers, numpatch));           a_trc_rsur      = 0._r8
      allocate(a_trc_rsub      (ntracers, numpatch));           a_trc_rsub      = 0._r8
      allocate(a_trc_rnof      (ntracers, numpatch));           a_trc_rnof      = 0._r8
      allocate(a_trc_qinfl     (ntracers, numpatch));           a_trc_qinfl     = 0._r8
      allocate(a_trc_qcharge   (ntracers, numpatch));           a_trc_qcharge   = 0._r8

      allocate(trc_storage_beg (ntracers, numpatch));           trc_storage_beg = 0._r8
      allocate(trc_balance_err (ntracers, numpatch));           trc_balance_err = 0._r8

      allocate(a_trc_ldew_mass (ntracers, numpatch));           a_trc_ldew_mass = 0._r8
      allocate(a_water_ldew    (numpatch));                     a_water_ldew    = 0._r8
      allocate(a_trc_soil_mass (ntracers, nl_soil, numpatch));  a_trc_soil_mass = 0._r8
      allocate(a_water_soil    (nl_soil, numpatch));             a_water_soil    = 0._r8
      trc_hist_nac = 0
   END SUBROUTINE allocate_Tracer_Vars

   SUBROUTINE deallocate_Tracer_Vars ()
      IMPLICIT NONE
      IF (allocated(trc_ldew_rain  )) deallocate(trc_ldew_rain  )
      IF (allocated(trc_ldew_snow  )) deallocate(trc_ldew_snow  )
      IF (allocated(trc_wliq_soisno)) deallocate(trc_wliq_soisno)
      IF (allocated(trc_wice_soisno)) deallocate(trc_wice_soisno)
      IF (allocated(trc_wa         )) deallocate(trc_wa         )
      IF (allocated(trc_wdsrf      )) deallocate(trc_wdsrf      )
      IF (allocated(trc_wetwat     )) deallocate(trc_wetwat     )
      IF (allocated(a_trc_precip   )) deallocate(a_trc_precip   )
      IF (allocated(a_trc_evap     )) deallocate(a_trc_evap     )
      IF (allocated(a_trc_trans    )) deallocate(a_trc_trans    )
      IF (allocated(a_trc_rsur     )) deallocate(a_trc_rsur     )
      IF (allocated(a_trc_rsub     )) deallocate(a_trc_rsub     )
      IF (allocated(a_trc_rnof     )) deallocate(a_trc_rnof     )
      IF (allocated(a_trc_qinfl    )) deallocate(a_trc_qinfl    )
      IF (allocated(a_trc_qcharge  )) deallocate(a_trc_qcharge  )
      IF (allocated(trc_storage_beg)) deallocate(trc_storage_beg)
      IF (allocated(trc_balance_err)) deallocate(trc_balance_err)
      IF (allocated(a_trc_ldew_mass)) deallocate(a_trc_ldew_mass)
      IF (allocated(a_water_ldew   )) deallocate(a_water_ldew   )
      IF (allocated(a_trc_soil_mass)) deallocate(a_trc_soil_mass)
      IF (allocated(a_water_soil   )) deallocate(a_water_soil   )
   END SUBROUTINE deallocate_Tracer_Vars

   SUBROUTINE flush_Tracer_Acc ()
      IMPLICIT NONE
      IF (allocated(a_trc_precip )) a_trc_precip  = 0._r8
      IF (allocated(a_trc_evap   )) a_trc_evap    = 0._r8
      IF (allocated(a_trc_trans  )) a_trc_trans   = 0._r8
      IF (allocated(a_trc_rsur   )) a_trc_rsur    = 0._r8
      IF (allocated(a_trc_rsub   )) a_trc_rsub    = 0._r8
      IF (allocated(a_trc_rnof   )) a_trc_rnof    = 0._r8
      IF (allocated(a_trc_qinfl  )) a_trc_qinfl   = 0._r8
      IF (allocated(a_trc_qcharge)) a_trc_qcharge  = 0._r8
      IF (allocated(a_trc_ldew_mass)) a_trc_ldew_mass = 0._r8
      IF (allocated(a_water_ldew   )) a_water_ldew    = 0._r8
      IF (allocated(a_trc_soil_mass)) a_trc_soil_mass = 0._r8
      IF (allocated(a_water_soil   )) a_water_soil    = 0._r8
      trc_hist_nac = 0
   END SUBROUTINE flush_Tracer_Acc

END MODULE MOD_Tracer_Vars
