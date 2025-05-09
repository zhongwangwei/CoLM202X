#include <define.h>
#ifdef BGC

MODULE MOD_BGC_CNBalanceCheck

!--------------------------------------------------------------------------------------------
! !DESCRIPTION:
! C and N balance check module.
! run sequential: BeginCNBalance(i) -> all CN cycle processes ->CBalanceCheck & NBalanceCheck
!
! !ORIGINAL:
! The Community Land Model version 5.0 (CLM5.0)
!
! !REVISION:
! Xingjie Lu, 2022, modify original CLM5 to be compatible with CoLM code structure.

   USE MOD_Precision
   USE MOD_Namelist, only: DEF_USE_NITRIF
   USE MOD_BGC_Vars_TimeVariables, only: &
       sminn, col_endcb, col_begcb, totcolc, col_endnb, col_begnb, totcoln, &
       col_vegbegcb, totvegc, ctrunc_veg, col_vegbegnb, totvegn, ntrunc_veg, &
       col_soilbegcb, totsomc, totlitc, totcwdc, ctrunc_soil, &
       col_soilbegnb, totsomn, totlitn, totcwdn, ntrunc_soil, col_sminnbegnb, &
       col_vegendcb, col_vegendnb, col_soilendcb, col_soilendnb, col_sminnendnb 
   USE MOD_BGC_Vars_1DFluxes, only: &
       gpp, er, ar, decomp_hr, fire_closs, hrv_xsmrpool_to_atm, wood_harvestc, grainc_to_cropprodc, &
       som_c_leached, ndep_to_sminn, nfix_to_sminn, supplement_to_sminn, ffix_to_sminn, &
       fert_to_sminn, soyfixn_to_sminn, denit, fire_nloss, wood_harvestn, grainn_to_cropprodn, &
       sminn_leached, f_n2o_nit, smin_no3_leached, smin_no3_runoff, som_n_leached, sminn_to_plant
   USE MOD_BGC_Vars_PFTimeVariables, only: &
       leafc_p, frootc_p, livestemc_p, deadstemc_p, livecrootc_p, deadcrootc_p, &
       leafc_storage_p, frootc_storage_p, livestemc_storage_p, &
       deadstemc_storage_p, livecrootc_storage_p, deadcrootc_storage_p, gresp_storage_p, &
       leafc_xfer_p, frootc_xfer_p, livestemc_xfer_p, &
       deadstemc_xfer_p, livecrootc_xfer_p, deadcrootc_xfer_p, gresp_xfer_p, xsmrpool_p, &
       grainc_p, grainc_storage_p, grainc_xfer_p, ctrunc_p, totvegc_p, cropseedc_deficit_p
   USE MOD_BGC_Vars_1DPFTFluxes, only: &
       grainc_to_food_p
   USE MOD_SPMD_Task
   USE MOD_Vars_PFTimeInvariants, only: pftclass
   
   IMPLICIT NONE
   
   PUBLIC BeginCNBalance
   PUBLIC CBalanceCheck
   PUBLIC NBalanceCheck

CONTAINS

   SUBROUTINE BeginCNBalance(i)
 
   ! !DESCRIPTION:
   ! BeginCNBalance SUBROUTINE stores initial C and N pool size at begining of each time step, which is
   ! further used in CN balance check.
   !
   ! !ORIGINAL:
   ! The Community Land Model version 5.0 (CLM5.0)
   !
   ! !REVISION:
   ! Xingjie Lu, 2022, modify original CLM5 to be compatible with CoLM code structure.
 
   integer, intent(in) :: i ! patch index
 
      col_begcb(i) = totcolc(i)
      col_begnb(i) = totcoln(i)
  
      col_vegbegcb(i) = totvegc(i) + ctrunc_veg(i)
      col_vegbegnb(i) = totvegn(i) + ntrunc_veg(i)
  
      col_soilbegcb(i) = totsomc(i) + totlitc(i) + totcwdc(i) + ctrunc_soil(i) 
      col_soilbegnb(i) = totsomn(i) + totlitn(i) + totcwdn(i) + ntrunc_soil(i)
  
      col_sminnbegnb(i) = sminn(i)
 
   END SUBROUTINE BeginCNBalance
 
   SUBROUTINE CBalanceCheck(i,ps,pe,deltim,dlat,dlon)
 
   ! !DESCRIPTION:
   ! CBalanceCheck tests the carbon balance of each time step, which meet C balance equation:
   ! col_endcb - col_begcb = (col_cinputs - col_coutputs)*deltim
   !
   ! !ORIGINAL:
   ! The Community Land Model version 5.0 (CLM5.0)
   !
   ! !REVISION:
   ! Xingjie Lu, 2022, modify original CLM5 to be compatible with CoLM code structure.
 
   integer, intent(in) :: i      ! patch index
   integer, intent(in) :: ps     ! start pft index
   integer, intent(in) :: pe     ! end pft index
   real(r8),intent(in) :: deltim ! time step in seconds
   real(r8),intent(in) :: dlat   ! latitude (degree)
   real(r8),intent(in) :: dlon   ! longitude (degree)
 
!Local variables
   real(r8),parameter :: cerror = 1.e-7_r8
   real(r8) :: col_cinputs, col_coutputs, col_errcb
   integer m
 
      col_endcb(i)     = totcolc(i)
      col_vegendcb(i)  = totvegc(i) + ctrunc_veg(i)
      col_soilendcb(i) = totsomc(i) + totlitc(i) + totcwdc(i) + ctrunc_soil(i)
  
      col_cinputs = gpp(i)
  
      col_coutputs = er(i) + fire_closs(i) + hrv_xsmrpool_to_atm(i) &
               + wood_harvestc(i) + grainc_to_cropprodc(i) - som_c_leached(i)
  
      col_errcb = (col_cinputs - col_coutputs)*deltim - &
                 (col_endcb(i) - col_begcb(i))
  
      IF(abs(col_errcb) > cerror) THEN
         write(*,*)'column cbalance error    = ', col_errcb, i, p_iam_glb
         write(*,*)'Latdeg,Londeg='             , dlat, dlon
         write(*,*)'begcb                    = ',col_begcb(i)
         write(*,*)'endcb                    = ',col_endcb(i)
         write(*,*)'delta store              = ',col_endcb(i)-col_begcb(i)
         write(*,*)'delta veg                = ',col_vegendcb(i) - col_vegbegcb(i),totvegc(i),col_vegendcb(i),col_vegbegcb(i)
         write(*,*)'delta soil               = ',col_soilendcb(i) - col_soilbegcb(i),totsomc(i),totlitc(i),totcwdc(i),col_soilendcb(i),col_soilbegcb(i)
         DO m = ps, pe
            write(*,*)'m=',m,pftclass(m)
            write(*,*)'vegc,leafc              = ',leafc_p(m)+leafc_storage_p(m)+leafc_xfer_p(m)
            write(*,*)'vegc,frootc             = ',frootc_p(m)+frootc_storage_p(m)+frootc_xfer_p(m)
            write(*,*)'vegc,livestemc          = ',livestemc_p(m)+livestemc_storage_p(m)+livestemc_xfer_p(m)
            write(*,*)'vegc,deadstemc          = ',deadstemc_p(m)+deadstemc_storage_p(m)+deadstemc_xfer_p(m)
            write(*,*)'vegc,livecrootc         = ',livecrootc_p(m)+livecrootc_storage_p(m)+livecrootc_xfer_p(m)
            write(*,*)'vegc,deadcrootc         = ',deadcrootc_p(m)+deadcrootc_storage_p(m)+deadcrootc_xfer_p(m)
            write(*,*)'grainc                  = ',grainc_p(m)+grainc_storage_p(m)+grainc_xfer_p(m)+cropseedc_deficit_p(m)
            write(*,*)'grwoth respiration c    = ',gresp_storage_p(m)+gresp_xfer_p(m)+xsmrpool_p(m)
         ENDDO
         write(*,*)'--- Inputs ---'
         write(*,*)'gpp                      = ',gpp(i)*deltim
         write(*,*)'--- Outputs ---'
         write(*,*)'er                       = ',er(i)*deltim
         write(*,*)'ar                       = ',ar(i)*deltim
         write(*,*)'decomp_hr                = ',decomp_hr(i)*deltim
         write(*,*)'fire_closs           = ',fire_closs(i)*deltim
         write(*,*)'col_hrv_xsmrpool_to_atm  = ',hrv_xsmrpool_to_atm(i)*deltim
         write(*,*)'wood_harvestc            = ',wood_harvestc(i)*deltim
         write(*,*)'grainc_to_cropprodc      = ',grainc_to_cropprodc(i)*deltim, grainc_to_food_p(ps)*deltim
         write(*,*)'-1*som_c_leached         = ',som_c_leached(i)*deltim
#ifdef USEMPI
         CALL mpi_abort (p_comm_glb, p_err)
#else
         CALL abort
#endif
      ENDIF

   END SUBROUTINE CBalanceCheck
 
   SUBROUTINE NBalanceCheck(i,deltim,dlat,dlon)
 
! !DESCRIPTION:
! NBalanceCheck tests the carbon balance of each time step, which meet N balance equation:
! col_endnb - col_begnb = (col_ninputs - col_noutputs)*deltim
!
! !ORIGINAL:
! The Community Land Model version 5.0 (CLM5.0)
!
! !REVISION:
! Xingjie Lu, 2022, modify original CLM5 to be compatible with CoLM code structure.
   
   integer, intent(in) :: i      ! patch index
   real(r8),intent(in) :: deltim ! time step in seconds
   real(r8),intent(in) :: dlat   ! latitude (degree)
   real(r8),intent(in) :: dlon   ! longitude (degree)
 
!Local variables
   real(r8),parameter :: nerror = 1.e-7_r8
   real(r8) :: col_ninputs, col_noutputs, col_errnb
 
      col_endnb(i) = totcoln(i)
      col_vegendnb(i)  = totvegn(i) + ntrunc_veg(i)
      col_soilendnb(i) = totsomn(i) + totlitn(i) + totcwdn(i) + ntrunc_soil(i)
      col_sminnendnb(i) = sminn(i)
  
      col_ninputs = ndep_to_sminn(i) + nfix_to_sminn(i) + supplement_to_sminn(i)
  
      col_ninputs = col_ninputs + fert_to_sminn(i) + soyfixn_to_sminn(i)
  
      col_noutputs = denit(i) + fire_nloss(i) + wood_harvestn(i) + grainn_to_cropprodn(i)
  
      IF(DEF_USE_NITRIF)THEN
         col_noutputs = col_noutputs + f_n2o_nit(i) + smin_no3_leached(i) + smin_no3_runoff(i)
      ELSE
         col_noutputs = col_noutputs + sminn_leached(i)
      ENDIF
  
      col_noutputs = col_noutputs - som_n_leached(i)
      col_errnb    =(col_ninputs - col_noutputs)*deltim - (col_endnb(i) - col_begnb(i))
  
      IF (abs(col_errnb) > nerror) THEN !
         write(*,*)'column nbalance error    = ',col_errnb, i, p_iam_glb
         write(*,*)'Latdeg,Londeg            = ',dlat, dlon
         write(*,*)'begnb                    = ',col_begnb(i)
         write(*,*)'endnb                    = ',col_endnb(i)
         write(*,*)'delta store              = ',col_endnb(i)-col_begnb(i)
         write(*,*)'delta veg                = ',col_vegendnb(i)-col_vegbegnb(i)
         write(*,*)'delta soil               = ',col_soilendnb(i)-col_soilbegnb(i)
         write(*,*)'delta sminn              = ',col_sminnendnb(i)-col_sminnbegnb(i)
         write(*,*)'smin_to_plant            = ',sminn_to_plant(i)*deltim
         write(*,*)'input mass               = ',col_ninputs*deltim
         write(*,*)'output mass              = ',col_noutputs*deltim,f_n2o_nit(i)*deltim,smin_no3_leached(i)*deltim,&
                                                 smin_no3_runoff(i)*deltim, denit(i)*deltim,fire_nloss(i)*deltim,&
                                                 ( wood_harvestn(i) + grainn_to_cropprodn(i))*deltim
         write(*,*)'net flux                 = ',(col_ninputs-col_noutputs)*deltim
         write(*,*)'inputs,ffix,nfix,ndep    = ',ffix_to_sminn(i)*deltim,nfix_to_sminn(i)*deltim,ndep_to_sminn(i)*deltim,&
                                            fert_to_sminn(i)*deltim,soyfixn_to_sminn(i)*deltim
         IF(DEF_USE_NITRIF)THEN
            write(*,*)'outputs,leached,runoff,denit = ',smin_no3_leached(i)*deltim, smin_no3_runoff(i)*deltim,f_n2o_nit(i)*deltim
         ELSE
            write(*,*)'outputs,leached,denit,fire,harvest,som_n_leached',&
                      sminn_leached(i)*deltim,denit(i)*deltim,fire_nloss(i)*deltim,&
                      (wood_harvestn(i)+grainn_to_cropprodn(i))*deltim, - som_n_leached(i)
         ENDIF
#ifdef USEMPI
         CALL mpi_abort (p_comm_glb, p_err)
#else
         CALL abort
#endif
      ENDIF

   END SUBROUTINE NBalanceCheck

END MODULE MOD_BGC_CNBalanceCheck
#endif
