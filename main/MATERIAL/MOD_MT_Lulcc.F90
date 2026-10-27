#include <define.h>

#if (defined MaterialTransport) && (defined LULCC)
MODULE MOD_MT_Lulcc
!-----------------------------------------------------------------------
! Material transport across a land cover change.
!
! LulccDriver calls two routines, so the host LULCC files keep their
! upstream layout:
!
!  mt_lulcc_save   before LulccInitialize, while landpatch/mesh/pixel still
!                  describe the old year. It snapshots the material pools
!                  and, for the MEC scheme, the physical area of every old
!                  patch.
!
!  mt_lulcc_remap  after the host variables are recovered (SAT or MEC) and
!                  before the transfer trace is freed. It moves the material
!                  pools, the provider pools and the forcing ratios onto the
!                  new patches and reloads provider inputs for the new year.
!
! With the MEC scheme the transfer trace says where each new patch's area
! came from. Before any pool is moved, the trace is checked against the old
! and new patch areas element by element, so a remap that could not conserve
! mass stops the run instead of silently creating or losing material.
!
! Moved here from the fork's MOD_Lulcc_Driver (#ifdef TRACER blocks).
!-----------------------------------------------------------------------

   USE MOD_Precision
   IMPLICIT NONE
   SAVE

   PUBLIC :: mt_lulcc_save
   PUBLIC :: mt_lulcc_remap

   PRIVATE

   real(r8), allocatable :: old_patch_area(:)
   logical :: inventory_active = .false.

CONTAINS

   SUBROUTINE mt_lulcc_save ()

   USE MOD_SPMD_Task
   USE MOD_Namelist,  only: DEF_LULCC_SCHEME
   USE MOD_LandPatch, only: landpatch, numpatch
   USE MOD_MT_Defs,   only: ntracers
   USE MOD_MT_Vars,   only: save_land_tracer_lulcc_state
   USE MOD_MT_Lifecycle, only: mt_lifecycle_land_save_lulcc_state
   IMPLICIT NONE

      CALL save_land_tracer_lulcc_state ()
      CALL mt_lifecycle_land_save_lulcc_state ()

      IF (allocated(old_patch_area)) deallocate(old_patch_area)
      inventory_active = p_is_worker .and. ntracers > 0 .and. DEF_LULCC_SCHEME == 2
      IF (inventory_active) THEN
         IF (numpatch > 0) THEN
            IF (.not. allocated(landpatch%eindex)) &
               CALL CoLM_stop('MT LULCC old patch map unavailable')
            IF (landpatch%nset /= numpatch) &
               CALL CoLM_stop('MT LULCC old patch count mismatch')
            CALL lulcc_patch_areas(landpatch, old_patch_area)
         ELSE
            allocate(old_patch_area(0))
         ENDIF
      ENDIF

   END SUBROUTINE mt_lulcc_save


   SUBROUTINE mt_lulcc_remap (lc_year, dir_landdata)

   USE MOD_SPMD_Task
   USE MOD_LandPatch, only: landpatch, numpatch
   USE MOD_Vars_TimeInvariants, only: patchclass
   USE MOD_Lulcc_Vars_TimeInvariants, only: landpatch_, patchclass_
   USE MOD_Lulcc_TransferTraceReadin, only: lccpct_patches
   USE MOD_MT_Vars,   only: remap_land_tracer_lulcc_state
   USE MOD_MT_Lifecycle, only: mt_lifecycle_land_remap_lulcc_state, &
                               mt_lifecycle_land_reload_lulcc_inputs
   USE MOD_MT_Conservation, only: deallocate_tracer_conservation
   USE MOD_MT_Forcing, only: mt_forcing_lulcc_remap
   IMPLICIT NONE

   integer, intent(in) :: lc_year
   character(len=*), intent(in) :: dir_landdata

   real(r8), allocatable :: new_patch_area(:), inventory_trace(:,:)

      IF (inventory_active) THEN
         IF ((size(old_patch_area) > 0) .neqv. (numpatch > 0)) &
            CALL CoLM_stop('MT LULCC worker element footprint changed')
         IF (size(old_patch_area) > 0) THEN
            IF (.not. allocated(patchclass_) .or. .not. allocated(landpatch_%eindex)) &
               CALL CoLM_stop('MT LULCC old patch map unavailable')
         ENDIF
         IF (.not. allocated(lccpct_patches)) &
            CALL CoLM_stop('MT LULCC MEC requires the transfer trace')
         IF (numpatch > 0) THEN
            IF (.not. allocated(landpatch%eindex)) &
               CALL CoLM_stop('MT LULCC new patch map unavailable')
            IF (landpatch%nset /= numpatch) &
               CALL CoLM_stop('MT LULCC new patch count mismatch')
            CALL lulcc_patch_areas(landpatch, new_patch_area)
         ELSE
            allocate(new_patch_area(0))
         ENDIF
         IF (allocated(patchclass_)) THEN
            IF (.not. allocated(patchclass)) &
               CALL CoLM_stop('MT LULCC MEC patch class map unavailable')
         ELSE
            IF (size(old_patch_area) > 0 .or. size(new_patch_area) > 0) &
               CALL CoLM_stop('MT LULCC MEC patch class map unavailable')
         ENDIF
         CALL lulcc_inventory_trace(lccpct_patches, inventory_trace)
         IF (size(old_patch_area) > 0 .or. size(new_patch_area) > 0) &
            CALL lulcc_check_inventory_transfer(patchclass_, landpatch_%eindex, old_patch_area, &
               patchclass, landpatch%eindex, new_patch_area, inventory_trace)
      ENDIF

      IF (p_is_worker .and. allocated(patchclass)) THEN
         IF (size(patchclass) > 0 .and. .not. allocated(patchclass_)) &
            CALL CoLM_stop('MT LULCC cannot remap a worker from zero to nonzero patches')
      ENDIF
      IF (p_is_worker .and. allocated(patchclass) .and. allocated(patchclass_) .and. &
          allocated(landpatch%eindex) .and. allocated(landpatch_%eindex)) THEN
         CALL deallocate_tracer_conservation ()
         IF (inventory_active) THEN
            CALL remap_land_tracer_lulcc_state (patchclass, landpatch%eindex, &
               patchclass_, landpatch_%eindex, inventory_trace, new_patch_area, old_patch_area)
            CALL mt_lifecycle_land_remap_lulcc_state (patchclass, landpatch%eindex, &
               patchclass_, landpatch_%eindex, inventory_trace, new_patch_area, old_patch_area)
            CALL mt_forcing_lulcc_remap (patchclass, landpatch%eindex, &
               patchclass_, landpatch_%eindex, lccpct_patches)
         ELSE
            CALL remap_land_tracer_lulcc_state (patchclass, landpatch%eindex, &
               patchclass_, landpatch_%eindex)
            CALL mt_lifecycle_land_remap_lulcc_state (patchclass, landpatch%eindex, &
               patchclass_, landpatch_%eindex)
            CALL mt_forcing_lulcc_remap (patchclass, landpatch%eindex, &
               patchclass_, landpatch_%eindex)
         ENDIF
      ENDIF

      IF (allocated(old_patch_area))  deallocate(old_patch_area)
      IF (allocated(new_patch_area))  deallocate(new_patch_area)
      IF (allocated(inventory_trace)) deallocate(inventory_trace)
      inventory_active = .false.

      CALL mt_lifecycle_land_reload_lulcc_inputs (lc_year, dir_landdata)

   END SUBROUTINE mt_lulcc_remap


   ! Physical area (m2) of every patch, from its pixels and shared fraction.
   SUBROUTINE lulcc_patch_areas(patches, areas)

   USE MOD_SPMD_Task, only: CoLM_stop
   USE MOD_Mesh,      only: mesh
   USE MOD_Pixel,     only: pixel
   USE MOD_Pixelset,  only: pixelset_type
   USE MOD_Utils,     only: areaquad
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   IMPLICIT NONE

   TYPE(pixelset_type), intent(in) :: patches
   real(r8), allocatable, intent(out) :: areas(:)
   integer :: ip, ie, ix, first, last

      allocate(areas(patches%nset))
      areas = 0._r8
      IF (patches%nset <= 0) RETURN
      IF (.not. allocated(mesh) .or. .not. allocated(patches%ielm) .or. &
          .not. allocated(patches%ipxstt) .or. .not. allocated(patches%ipxend) .or. &
          .not. allocated(pixel%lat_s) .or. .not. allocated(pixel%lat_n) .or. &
          .not. allocated(pixel%lon_w) .or. .not. allocated(pixel%lon_e)) &
         CALL CoLM_stop('MT LULCC patch geometry unavailable')
      IF (size(patches%ielm) /= patches%nset .or. size(patches%ipxstt) /= patches%nset .or. &
          size(patches%ipxend) /= patches%nset) &
         CALL CoLM_stop('MT LULCC patch geometry shape mismatch')
      IF (patches%has_shared) THEN
         IF (.not. allocated(patches%pctshared)) &
            CALL CoLM_stop('MT LULCC missing shared patch fractions')
         IF (size(patches%pctshared) /= patches%nset) &
            CALL CoLM_stop('MT LULCC shared patch area shape mismatch')
      ENDIF
      DO ip = 1, patches%nset
         ie = patches%ielm(ip)
         IF (ie < 1 .or. ie > size(mesh)) CALL CoLM_stop('MT LULCC invalid patch element')
         IF (.not. allocated(mesh(ie)%ilat) .or. .not. allocated(mesh(ie)%ilon)) &
            CALL CoLM_stop('MT LULCC element pixels unavailable')
         IF (size(mesh(ie)%ilat) < mesh(ie)%npxl .or. size(mesh(ie)%ilon) < mesh(ie)%npxl) &
            CALL CoLM_stop('MT LULCC element pixel shape mismatch')
         first = patches%ipxstt(ip)
         last = patches%ipxend(ip)
         IF (first == -1 .and. last == -1) CYCLE
         IF (first < 1 .or. last > mesh(ie)%npxl .or. last < first) &
            CALL CoLM_stop('MT LULCC invalid patch pixel range')
         DO ix = first, last
            IF (mesh(ie)%ilat(ix) < 1 .or. mesh(ie)%ilat(ix) > size(pixel%lat_s) .or. &
                mesh(ie)%ilat(ix) > size(pixel%lat_n) .or. mesh(ie)%ilon(ix) < 1 .or. &
                mesh(ie)%ilon(ix) > size(pixel%lon_w) .or. &
                mesh(ie)%ilon(ix) > size(pixel%lon_e)) &
               CALL CoLM_stop('MT LULCC invalid pixel coordinate index')
            areas(ip) = areas(ip) + 1.e6_r8 * areaquad( &
               pixel%lat_s(mesh(ie)%ilat(ix)), pixel%lat_n(mesh(ie)%ilat(ix)), &
               pixel%lon_w(mesh(ie)%ilon(ix)), pixel%lon_e(mesh(ie)%ilon(ix)))
         ENDDO
         IF (patches%has_shared) THEN
            IF (.not. ieee_is_finite(patches%pctshared(ip))) &
               CALL CoLM_stop('MT LULCC invalid shared patch fraction')
            IF (patches%pctshared(ip) < 0._r8 .or. patches%pctshared(ip) > 1._r8) &
               CALL CoLM_stop('MT LULCC invalid shared patch fraction')
            areas(ip) = areas(ip) * patches%pctshared(ip)
         ENDIF
         IF (.not. ieee_is_finite(areas(ip))) &
            CALL CoLM_stop('MT LULCC invalid physical patch area')
         IF (areas(ip) < 0._r8) &
            CALL CoLM_stop('MT LULCC invalid physical patch area')
      ENDDO

   END SUBROUTINE lulcc_patch_areas


   ! The transfer trace is read per land cover class. With PFT/PC all
   ! vegetated classes share one soil patch (class 1 here), and with the fast
   ! PC option cropland and class 14 share the crop patch, so fold the trace
   ! onto the classes the patches actually carry.
   SUBROUTINE lulcc_inventory_trace(raw, mapped)

   USE MOD_SPMD_Task, only: CoLM_stop
   USE MOD_Const_LC,  only: patchtypes
   USE MOD_Vars_Global, only: N_land_classification, CROPLAND
   USE MOD_Namelist,  only: DEF_USE_PFT, DEF_SOLO_PFT, DEF_FAST_PC
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   IMPLICIT NONE

   real(r8), intent(in) :: raw(:,0:)
   real(r8), allocatable, intent(out) :: mapped(:,:)
   integer :: c, dest

      IF (ubound(raw,2) /= N_land_classification) &
         CALL CoLM_stop('MT LULCC transfer trace class count mismatch')
      IF (any(.not. ieee_is_finite(raw))) &
         CALL CoLM_stop('MT LULCC invalid transfer trace')
      IF (any(raw < 0._r8)) &
         CALL CoLM_stop('MT LULCC invalid transfer trace')
      allocate(mapped(size(raw,1),0:N_land_classification))
      mapped = raw
      IF (.not. ((DEF_USE_PFT .and. .not. DEF_SOLO_PFT) .or. DEF_FAST_PC)) RETURN
      mapped = 0._r8
      mapped(:,0) = raw(:,0)
      DO c = 1, N_land_classification
         dest = c
         IF (patchtypes(c) == 0) dest = 1
         IF (DEF_FAST_PC .and. (c == CROPLAND .or. c == 14)) dest = CROPLAND
         mapped(:,dest) = mapped(:,dest) + raw(:,c)
      ENDDO

   END SUBROUTINE lulcc_inventory_trace


   ! Per element: the old area of each class must equal the new area that the
   ! trace says came from that class, and the element's total area must not
   ! change.
   SUBROUTINE lulcc_check_inventory_transfer(old_class, old_element, old_area, &
      new_class, new_element, new_area, trace)

   USE MOD_SPMD_Task, only: CoLM_stop
   USE MOD_Utils,     only: quicksort
   USE MOD_Vars_Global, only: N_land_classification
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   IMPLICIT NONE

   integer, intent(in) :: old_class(:), new_class(:)
   integer*8, intent(in) :: old_element(:), new_element(:)
   real(r8), intent(in) :: old_area(:), new_area(:), trace(:,0:)
   real(r8) :: source(0:N_land_classification), target(0:N_land_classification)
   real(r8) :: old_total, new_total, tol, rowsum
   integer*8, allocatable :: old_sorted(:), new_sorted(:)
   integer, allocatable :: old_order(:), new_order(:)
   integer :: os, oe, ns, ne, op, np, c, k

      IF (size(old_class) /= size(old_area) .or. size(old_class) /= size(old_element) .or. &
          size(new_class) /= size(new_area) .or. size(new_class) /= size(new_element) .or. &
          size(trace,1) /= size(new_class) .or. ubound(trace,2) /= N_land_classification) &
         CALL CoLM_stop('MT LULCC inventory map shape mismatch')
      IF (any(.not. ieee_is_finite(old_area)) .or. any(.not. ieee_is_finite(new_area)) .or. &
          any(.not. ieee_is_finite(trace))) &
         CALL CoLM_stop('MT LULCC invalid inventory area or transfer trace')
      IF (any(old_area < 0._r8) .or. any(new_area < 0._r8) .or. any(trace < 0._r8)) &
         CALL CoLM_stop('MT LULCC invalid inventory area or transfer trace')
      old_sorted = old_element
      new_sorted = new_element
      old_order = [(k, k=1,size(old_class))]
      new_order = [(k, k=1,size(new_class))]
      CALL quicksort(size(old_class), old_sorted, old_order)
      CALL quicksort(size(new_class), new_sorted, new_order)
      os = 1
      ns = 1
      DO WHILE (os <= size(old_class) .or. ns <= size(new_class))
         IF (os > size(old_class) .or. ns > size(new_class)) &
            CALL CoLM_stop('MT LULCC element footprint changed')
         IF (old_sorted(os) /= new_sorted(ns)) &
            CALL CoLM_stop('MT LULCC element footprint changed')
         oe = os
         DO WHILE (oe < size(old_class))
            IF (old_sorted(oe+1) /= old_sorted(os)) EXIT
            oe = oe + 1
         ENDDO
         ne = ns
         DO WHILE (ne < size(new_class))
            IF (new_sorted(ne+1) /= new_sorted(ns)) EXIT
            ne = ne + 1
         ENDDO
         source = 0._r8
         target = 0._r8
         DO k = os, oe
            op = old_order(k)
            IF (old_area(op) <= tiny(1._r8)) CYCLE
            c = old_class(op)
            IF (c < 0 .or. c > N_land_classification) &
               CALL CoLM_stop('MT LULCC invalid old patch class')
            source(c) = source(c) + old_area(op)
         ENDDO
         DO k = ns, ne
            np = new_order(k)
            IF (new_area(np) <= tiny(1._r8)) CYCLE
            c = new_class(np)
            IF (c < 0 .or. c > N_land_classification) &
               CALL CoLM_stop('MT LULCC invalid new patch class')
            rowsum = sum(trace(np,:))
            IF (abs(rowsum - 1._r8) > 1.e-10_r8) &
               CALL CoLM_stop('MT LULCC incomplete source trace')
            target = target + new_area(np) * trace(np,:) / rowsum
         ENDDO
         old_total = sum(source)
         new_total = sum(new_area(new_order(ns:ne)))
         tol = 1.e-10_r8 * max(old_total, new_total, 1._r8)
         IF (abs(old_total - new_total) > tol) THEN
            WRITE(*,'(A,I0,2(A,ES16.8))') 'MT LULCC element=', old_sorted(os), &
               ' old_area_m2=', old_total, ' new_area_m2=', new_total
            CALL CoLM_stop('MT LULCC element physical footprint changed')
         ENDIF
         DO c = 0, N_land_classification
            IF (((source(c) > 0._r8) .neqv. (target(c) > 0._r8)) .or. &
                abs(source(c) - target(c)) > 1.e-10_r8 * max(source(c), target(c))) THEN
               WRITE(*,'(A,I0,A,I0,2(A,ES16.8))') 'MT LULCC element=', old_sorted(os), &
                  ' class=', c, ' source_area_m2=', source(c), ' target_area_m2=', target(c)
               CALL CoLM_stop('MT LULCC source-class physical area mismatch')
            ENDIF
         ENDDO
         os = oe + 1
         ns = ne + 1
      ENDDO

   END SUBROUTINE lulcc_check_inventory_transfer

END MODULE MOD_MT_Lulcc
#endif
