#include <define.h>

SUBROUTINE CoLMDRIVER (idate,deltim,dolai,doalb,dosst,oro)


!=======================================================================
!
!  CoLM MODEL DRIVER
!
!  Initial : Yongjiu Dai, 1999-2014
!  Revised : Hua Yuan, Shupeng Zhang, Nan Wei, Xingjie Lu, Zhongwang Wei, Yongjiu Dai
!            2014-2024
!
!=======================================================================

   USE MOD_Precision
   USE MOD_Const_Physical, only: tfrz, rgas, vonkar
   USE MOD_Const_LC
   USE MOD_Vars_Global
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_TimeVariables
   USE MOD_Vars_1DForcing
   USE MOD_Vars_1DFluxes
   USE MOD_LandPatch, only: numpatch,landpatch
   USE MOD_LandUrban, only: patch2urban
   USE MOD_Namelist, only: DEF_forcing, DEF_URBAN_RUN, DEF_USE_TRACER
   USE MOD_Forcing, only: forcmask_pch
   USE omp_lib
#ifdef HYPERSPECTRAL
  USE MOD_HighRes_Parameters
#endif
#ifdef CaMa_Flood
   ! get flood variables: inundation depth[mm], inundation fraction [0-1],
   ! inundation evaporation [mm/s], inundation re-infiltration[mm/s]
   USE MOD_CaMa_Vars, only: flddepth_cama,fldfrc_cama,fevpg_fld,finfg_fld
#endif

   IMPLICIT NONE

   integer,  intent(in) :: idate(3) ! model calendar for next time step (year, julian day, seconds)
   real(r8), intent(in) :: deltim   ! seconds in a time-step

   logical,  intent(in) :: dolai    ! true if time for time-varying vegetation parameter
   logical,  intent(in) :: doalb    ! true if time for surface albedo calculation
   logical,  intent(in) :: dosst    ! true if time for update sst/ice/snow

   real(r8), intent(inout) :: oro(numpatch)  ! ocean(0)/seaice(2)/ flag

   real(r8) :: deltim_phy
   integer  :: steps_in_one_deltim
   integer  :: i, m, u, k

! ======================================================================

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i, m, u, k, steps_in_one_deltim, deltim_phy) &
!$OMP SCHEDULE(STATIC, 1)
#endif

      DO i = 1, numpatch

         ! Apply forcing mask
         IF (DEF_forcing%has_missing_value) THEN
            IF (.not. forcmask_pch(i)) CYCLE
         ENDIF

         ! Apply patch mask, but still run virtual 2m WMO patch (patch ipxstt=-1)
         IF (DEF_Output_2mWMO) THEN
            IF (.not. patchmask(i) .and. (landpatch%ipxstt(i)>0) ) CYCLE
         ELSE
            IF (.not. patchmask(i)) CYCLE
         ENDIF

         m = patchclass(i)

         steps_in_one_deltim = 1
         ! deltim need to be within 1800s for water body with snow in order to avoid large
         ! temperature fluctuations due to rapid snow heat conductance
         IF(m == WATERBODY) steps_in_one_deltim = ceiling(deltim/1800.)
         deltim_phy = deltim/steps_in_one_deltim

         ! For non urban patch or slab urban
         IF (.not.DEF_URBAN_RUN .or. m.ne.URBAN) THEN

            DO k = 1, steps_in_one_deltim
               !                ***** Call CoLM main program *****
               !
               CALL CoLMMAIN (  i,idate,         coszen(i),       deltim_phy,      &
               patchlonr(i),    patchlatr(i),    patchclass(i),   patchtype(i),    &
               doalb,           dolai,           dosst,           oro(i),          &

             ! SOIL INFORMATION AND LAKE DEPTH
               soil_s_v_alb(i), soil_d_v_alb(i), soil_s_n_alb(i), soil_d_n_alb(i), &
               vf_quartz(1:,i), vf_gravels(1:,i),vf_om(1:,i),     vf_sand(1:,i),   &
               wf_gravels(1:,i),wf_sand(1:,i),   porsl(1:,i),     psi0(1:,i),      &
               bsw(1:,i),       theta_r(1:,i),   fsatmax(i),      fsatdcf(i),      &
               topoweti(i),     alp_twi(i),      chi_twi(i),      mu_twi(i),       &
#ifdef vanGenuchten_Mualem_SOIL_MODEL
               alpha_vgm(1:,i), n_vgm(1:,i),     L_vgm(1:,i),                      &
               sc_vgm(1:,i),    fc_vgm(1:,i),                                      &
#endif
               hksati(1:,i),    csol(1:,i),      k_solids(1:,i),  dksatu(1:,i),    &
               dksatf(1:,i),    dkdry(1:,i),     BA_alpha(1:,i),  BA_beta(1:,i),   &
               rootfr(1:,m),    lakedepth(i),    dz_lake(1:,i),   elvstd(i),       &
               BVIC(i),                                                            &
#if (defined CaMa_Flood)
             ! flood variables [mm, m2/m2, mm/s, mm/s]
               flddepth_cama(i),fldfrc_cama(i),  fevpg_fld(i),    finfg_fld(i),    &
#endif

             ! VEGETATION INFORMATION
               htop(i),         hbot(i),         sqrtdi(m),                        &
               effcon(m),       vmax25(m),       c3c4(m),                          &
               kmax_sun(m),     kmax_sha(m),     kmax_xyl(m),     kmax_root(m),    &
               psi50_sun(m),    psi50_sha(m),    psi50_xyl(m),    psi50_root(m),   &
               ck(m),           slti(m),         hlti(m),         shti(m),         &
               hhti(m),         trda(m),         trdm(m),         trop(m),         &
               g1(m),           g0(m),           gradm(m),        binter(m),       &
               extkn(m),        chil(m),         rho(1:,1:,m),    tau(1:,1:,m),    &
#ifdef HYPERSPECTRAL
               ! variables for high spectral resolution
               ! note that rho & tau are depend on wavelength
               clr_frac,        cld_frac,                                          &
               reflectance(0:,1:,1:),            transmittance(0:,1:,1:),          &
               soil_alb(1:,i),    kw(1:),          nw(1:),                           &
#endif

             ! ATMOSPHERIC FORCING
               forc_pco2m(i),   forc_po2m(i),    forc_us(i),      forc_vs(i),      &
               forc_t(i),       forc_q(i),       forc_prc(i),     forc_prl(i),     &
               forc_rain(i),    forc_snow(i),    forc_psrf(i),    forc_pbot(i),    &
               forc_sols(i),    forc_soll(i),    forc_solsd(i),   forc_solld(i),   &
               forc_frl(i),     forc_hgt_u(i),   forc_hgt_t(i),   forc_hgt_q(i),   &
               forc_rhoair(i),                                                     &
#ifdef HYPERSPECTRAL
             ! solar forcing
               forc_solarin(i),                                                    &
#endif
             ! CBL height forcing
               forc_hpbl(i),                                                       &
             ! Aerosol deposition
               forc_aerdep(:,i),                                                   &

             ! LAND SURFACE VARIABLES REQUIRED FOR RESTART
               z_sno(maxsnl+1:,i),               dz_sno(maxsnl+1:,i),              &
               t_soisno(maxsnl+1:,i),            wliq_soisno(maxsnl+1:,i),         &
               wice_soisno(maxsnl+1:,i),         smp(1:,i),       hk(1:,i),        &
               t_grnd(i),       tleaf(i),        ldew(i),         ldew_rain(i),    &
               ldew_snow(i),    fwet_snow(i),    sag(i),          scv(i),          &
               snowdp(i),       fveg(i),         fsno(i),         sigf(i),         &
               green(i),        lai(i),          sai(i),          alb(1:,1:,i),    &
               ssun(1:,1:,i),   ssha(1:,1:,i),   ssoi(:,:,i),     ssno(:,:,i),     &
               thermk(i),       extkb(i),        extkd(i),        vegwp(1:,i),     &
               gs0sun(i),       gs0sha(i),       &
#ifdef HYPERSPECTRAL
               ! high-res variables
               alb_hires(1:,1:,i),               &
               sol_dir_ln_hires(1:,i), sol_dif_ln_hires(1:,i) ,&
               sr_dir_ln_hires(1:,i) , sr_dif_ln_hires (1:,i) ,&
               reflectance_out(:,:,i), transmittance_out(:,:,i),&
#endif
             ! Ozone Stress Variables
               o3coefv_sun(i),  o3coefv_sha(i),  o3coefg_sun(i),  o3coefg_sha(i),  &
               lai_old(i),      o3uptakesun(i),  o3uptakesha(i),  forc_ozone(i),   &
             ! End ozone stress variables
             ! WUE stomata model parameter
               lambda(m),                                                          &
             ! End WUE model parameter
               zwt(i),          wdsrf(i),        wa(i),           wetwat(i),       &
               t_lake(1:,i),    lake_icefrac(1:,i),               savedtke1(i),    &

             ! SNICAR snow model related
               snw_rds(:,i),    ssno_lyr(:,:,:,i),                                 &
               mss_bcpho(:,i),  mss_bcphi(:,i),  mss_ocpho(:,i),  mss_ocphi(:,i),  &
               mss_dst1(:,i),   mss_dst2(:,i),   mss_dst3(:,i),   mss_dst4(:,i),   &

             ! additional diagnostic variables for output
               laisun(i),       laisha(i),       rootr(1:,i),     rootflux(1:,i),  &
               rstfacsun_out(i),rstfacsha_out(i),gssun_out(i),    gssha_out(i),    &
               assimsun_out(i), etrsun_out(i),   assimsha_out(i), etrsha_out(i),   &
               h2osoi(1:,i),    wat(i),          rss(i),          &

             ! FLUXES
               taux(i),         tauy(i),         fsena(i),        fevpa(i),        &
               lfevpa(i),       fsenl(i),        fevpl(i),        etr(i),          &
               fseng(i),        fevpg(i),        olrg(i),         fgrnd(i),        &
               trad(i),         tref(i),         qref(i),         t2m_wmo(i),      &
               frcsat(i),       rsur(i),         rsur_se(i),      rsur_ie(i),      &
               rsub(i),                                                            &
               rnof(i),         qintr(i),        qinfl(i),        qlayer(0:,i),    &
               lake_deficit(i), qdrip(i),        rst(i),          assim(i),        &
               respc(i),        sabvsun(i),      sabvsha(i),      sabg(i),         &
               sr(i),           solvd(i),        solvi(i),        solnd(i),        &
               solni(i),        srvd(i),         srvi(i),         srnd(i),         &
               srni(i),         solvdln(i),      solviln(i),      solndln(i),      &
               solniln(i),      srvdln(i),       srviln(i),       srndln(i),       &
               srniln(i),       qcharge(i),      xerr(i),         zerr(i),         &

             ! TUNABLE modle constants
               zlnd,            zsno,            csoilc,          dewmx,           &
               ! 'wtfact' is updated to gridded 'fsatmax' data.
               capr,            cnfac,           ssi,             wimp,            &
               pondmx,          smpmax,          smpmin,          trsmx0,          &
               tcrit,           &

             ! additional variables required by coupling with WRF model
               emis(i),         z0m(i),          zol(i),          rib(i),          &
               ustar(i),        qstar(i),        tstar(i),                         &
               fm(i),           fh(i),           fq(i)                             )

            ENDDO
         ENDIF


#if (defined BGC)
         IF(patchtype(i) .eq. 0)THEN
            !
            !               ***** Call CoLM BGC model *****
            !
            CALL bgc_driver (i,idate(1:3),deltim, patchlatr(i)*180/PI,patchlonr(i)*180/PI)
         ENDIF
#endif


#ifdef URBAN_MODEL
         ! For urban model and urban patches
         IF (DEF_URBAN_RUN .and. m.eq.URBAN) THEN

            u = patch2urban(i)
            !
            !              ***** Call CoLM urban model *****
            !
            CALL CoLMMAIN_Urban ( &
          ! MODEL RUNNING PARAMETERS
            i               ,idate           ,coszen(i)       ,deltim          ,&
            patchlonr(i)    ,patchlatr(i)    ,patchclass(i)   ,patchtype(i)    ,&

          ! URBAN PARAMETERS
            froof(u)        ,flake(u)        ,hroof(u)        ,hlr(u)          ,&
            fgper(u)        ,em_roof(u)      ,em_wall(u)      ,em_gimp(u)      ,&
            em_gper(u)      ,cv_roof(:,u)    ,cv_wall(:,u)    ,cv_gimp(:,u)    ,&
            tk_roof(:,u)    ,tk_wall(:,u)    ,tk_gimp(:,u)    ,z_roof(:,u)     ,&
            z_wall(:,u)     ,dz_roof(:,u)    ,dz_wall(:,u)    ,lakedepth(i)    ,&
            dz_lake(1:,i)   ,elvstd(i)       ,BVIC(i)                          ,&

          ! LUCY INPUT PARAMETERS
            fix_holiday(:,u),week_holiday(:,u),hum_prof(:,u)  ,pop_den(u)      ,&
            vehicle(:,u)    ,weh_prof(:,u)   ,wdh_prof(:,u)                    ,&

          ! SOIL INFORMATION AND LAKE DEPTH
            vf_quartz(1:,i) ,vf_gravels(1:,i),vf_om(1:,i)     ,vf_sand(1:,i)   ,&
            wf_gravels(1:,i),wf_sand(1:,i)   ,porsl(1:,i)     ,psi0(1:,i)      ,&
            bsw(1:,i)       ,theta_r(1:,i)   ,fsatmax(i)      ,fsatdcf(i)      ,&
#ifdef vanGenuchten_Mualem_SOIL_MODEL
            alpha_vgm(1:,i) ,n_vgm(1:,i)     ,L_vgm(1:,i)                      ,&
            sc_vgm(1:,i)    ,fc_vgm(1:,i)                                      ,&
#endif
            hksati(1:,i)    ,csol(1:,i)      ,k_solids(1:,i)  ,dksatu(1:,i)    ,&
            dksatf(1:,i)    ,dkdry(1:,i)     ,BA_alpha(1:,i)  ,BA_beta(1:,i)   ,&
            alb_roof(:,:,u) ,alb_wall(:,:,u) ,alb_gimp(:,:,u) ,alb_gper(:,:,u) ,&

          ! VEGETATION INFORMATION
            htop(i)         ,hbot(i)         ,sqrtdi(m)       ,chil(m)         ,&
            effcon(m)       ,vmax25(m)       ,c3c4(m)         ,slti(m)         ,hlti(m)         ,&
            shti(m)         ,hhti(m)         ,trda(m)         ,trdm(m)         ,&
            trop(m)         ,g1(m)           ,g0(m)           ,gradm(m)        ,&
            binter(m)       ,extkn(m)        ,rho(1:,1:,m)    ,tau(1:,1:,m)    ,&
            rootfr(1:,m)    ,&
          ! WUE model parameter
            lambda(m)                                                          ,&
          ! END WUE model parameter

          ! ATMOSPHERIC FORCING
            forc_pco2m(i)   ,forc_po2m(i)    ,forc_us(i)      ,forc_vs(i)      ,&
            forc_t(i)       ,forc_q(i)       ,forc_prc(i)     ,forc_prl(i)     ,&
            forc_rain(i)    ,forc_snow(i)    ,forc_psrf(i)    ,forc_pbot(i)    ,&
            forc_sols(i)    ,forc_soll(i)    ,forc_solsd(i)   ,forc_solld(i)   ,&
            forc_frl(i)     ,forc_hgt_u(i)   ,forc_hgt_t(i)   ,forc_hgt_q(i)   ,&
            forc_rhoair(i)  ,Fhac(u)         ,Fwst(u)         ,Fach(u)         ,&
            Fahe(u)         ,Fhah(u)         ,vehc(u)         ,meta(u)         ,&

          ! LAND SURFACE VARIABLES REQUIRED FOR RESTART
            z_sno_roof  (maxsnl+1:,u)        ,z_sno_gimp  (maxsnl+1:,u)        ,&
            z_sno_gper  (maxsnl+1:,u)        ,z_sno_lake  (maxsnl+1:,u)        ,&
            dz_sno_roof (maxsnl+1:,u)        ,dz_sno_gimp (maxsnl+1:,u)        ,&
            dz_sno_gper (maxsnl+1:,u)        ,dz_sno_lake (maxsnl+1:,u)        ,&
            t_roofsno   (maxsnl+1:,u)        ,t_gimpsno   (maxsnl+1:,u)        ,&
            t_gpersno   (maxsnl+1:,u)        ,t_lakesno   (maxsnl+1:,u)        ,&
            wliq_roofsno(maxsnl+1:,u)        ,wliq_gimpsno(maxsnl+1:,u)        ,&
            wliq_gpersno(maxsnl+1:,u)        ,wliq_lakesno(maxsnl+1:,u)        ,&
            wice_roofsno(maxsnl+1:,u)        ,wice_gimpsno(maxsnl+1:,u)        ,&
            wice_gpersno(maxsnl+1:,u)        ,wice_lakesno(maxsnl+1:,u)        ,&
            z_sno       (maxsnl+1:,i)        ,dz_sno      (maxsnl+1:,i)        ,&
            wliq_soisno (maxsnl+1:,i)        ,wice_soisno (maxsnl+1:,i)        ,&
            t_soisno    (maxsnl+1:,i)        ,&
            smp         (1:,i)               ,hk          (1:,i)               ,&
            t_wallsun   (1:,u)               ,t_wallsha   (1:,u)               ,&

            lai(i)          ,sai(i)          ,fveg(i)         ,sigf(i)         ,&
            green(i)        ,tleaf(i)        ,ldew(i)         ,ldew_rain(i)    ,&
            ldew_snow(i)    ,fwet_snow(i)    ,t_grnd(i)                        ,&

            sag_roof(u)     ,sag_gimp(u)     ,sag_gper(u)     ,sag_lake(u)     ,&
            scv_roof(u)     ,scv_gimp(u)     ,scv_gper(u)     ,scv_lake(u)     ,&
            snowdp_roof(u)  ,snowdp_gimp(u)  ,snowdp_gper(u)  ,snowdp_lake(u)  ,&
            fsno_roof(u)    ,fsno_gimp(u)    ,fsno_gper(u)    ,fsno_lake(u)    ,&
            sag(i)          ,scv(i)          ,snowdp(i)       ,fsno(i)         ,&
            extkd(i)        ,alb(1:,1:,i)    ,ssun(1:,1:,i)   ,ssha(1:,1:,i)   ,&
            sroof(1:,1:,u)  ,swsun(1:,1:,u)  ,swsha(1:,1:,u)  ,sgimp(1:,1:,u)  ,&
            sgper(1:,1:,u)  ,slake(1:,1:,u)  ,lwsun(u)        ,lwsha(u)        ,&
            lgimp(u)        ,lgper(u)        ,lveg(u)         ,fwsun(u)        ,&
            dfwsun(u)       ,t_room(u)       ,troof_inner(u)  ,twsun_inner(u)  ,&
            twsha_inner(u)  ,t_roommax(u)    ,t_roommin(u)    ,tafu(u)         ,&

            zwt(i)          ,wdsrf(i)        ,wa(i)                            ,&
            t_lake(1:,i)    ,lake_icefrac(1:,i),               savedtke1(i)    ,&

          ! SNICAR snow model related
            snw_rds(:,i)    ,ssno_lyr(:,:,:,i)                                 ,&
            mss_bcpho(:,i)  ,mss_bcphi(:,i)  ,mss_ocpho(:,i)  ,mss_ocphi(:,i)  ,&
            mss_dst1(:,i)   ,mss_dst2(:,i)   ,mss_dst3(:,i)   ,mss_dst4(:,i)   ,&

#if (defined CaMa_Flood)
          ! flood variables [mm, m2/m2, mm/s, mm/s]
            flddepth_cama(i),fldfrc_cama(i)  ,fevpg_fld(i)    ,finfg_fld(i)    ,&
#endif

          ! additional diagnostic variables for output
            laisun(i)       ,laisha(i)       ,rss(i)                           ,&
            rstfacsun_out(i),h2osoi(1:,i)    ,wat(i)                           ,&

          ! FLUXES
            taux(i)         ,tauy(i)         ,fsena(i)        ,fevpa(i)        ,&
            lfevpa(i)       ,fsenl(i)        ,fevpl(i)        ,etr(i)          ,&
            fseng(i)        ,fevpg(i)        ,olrg(i)         ,fgrnd(i)        ,&
            fsen_roof(u)    ,fsen_wsun(u)    ,fsen_wsha(u)    ,fsen_gimp(u)    ,&
            fsen_gper(u)    ,fsen_urbl(u)    ,t_roof(u)       ,t_wall(u)       ,&
            lfevp_roof(u)   ,lfevp_gimp(u)   ,lfevp_gper(u)   ,lfevp_urbl(u)   ,&
            trad(i)         ,tref(i)         ,&
            qref(i)         ,rsur(i)         ,rnof(i)         ,qintr(i)        ,&
            qinfl(i)        ,qdrip(i)        ,rst(i)          ,assim(i)        ,&
            respc(i)        ,sabvsun(i)      ,sabvsha(i)      ,sabg(i)         ,&
            sr(i)           ,solvd(i)        ,solvi(i)        ,solnd(i)        ,&
            solni(i)        ,srvd(i)         ,srvi(i)         ,srnd(i)         ,&
            srni(i)         ,solvdln(i)      ,solviln(i)      ,solndln(i)      ,&
            solniln(i)      ,srvdln(i)       ,srviln(i)       ,srndln(i)       ,&
            srniln(i)       ,qcharge(i)      ,xerr(i)         ,zerr(i)         ,&

          ! TUNABLE model constants
            zlnd            ,zsno            ,csoilc          ,dewmx           ,&
            ! 'wtfact' is updated to gridded 'fsatmax' data.
            capr            ,cnfac           ,ssi             ,wimp            ,&
            pondmx          ,smpmax          ,smpmin          ,trsmx0          ,&
            tcrit                                                              ,&

          ! additional variables required by coupling with WRF model
            emis(i)         ,z0m(i)          ,zol(i)          ,rib(i)          ,&
            ustar(i)        ,qstar(i)        ,tstar(i)        ,fm(i)           ,&
            fh(i)           ,fq(i)           ,forc_hpbl(i)                      )
            rsub(i) = rnof(i) - rsur(i)

            ! ---------------------------------------------------------------
            ! Urban tracer handling (Phase 1 approximation).
            !
            ! Urban patches invoke CoLMMAIN_Urban which bypasses
            ! tracer_precip / _evapo / _soil_water. Under Phase 1
            ! (constant atmospheric R = R_init, no fractionation) the
            ! pool ratio must be R_init at every step by mass balance,
            ! so the prognostic pools can be re-synced from the current
            ! water state (sync_tracer_patch_phase1). This keeps the
            ! per-pool δ diagnostics physical instead of drifting to
            ! ±1000‰ as the water side evolves freely.
            !
            ! Urban has real canopy (ldew_rain/snow) — pass them to
            ! the sync via the optional args so trc_ldew is rebuilt
            ! from current ldew, not zeroed like glacier/lake.
            !
            ! The same save_storage → sync → accumulators →
            ! balance_check → hist_accumulate chain used for glacier/
            ! lake applies here (see CoLMMAIN.F90 glacier/lake blocks
            ! for the equivalent plumbing). Iterate over maxsnl+1:
            ! nl_soil (pass maxsnl to save_storage / balance_check) so
            ! in-step snow-layer count changes inside CoLMMAIN_Urban
            ! do not cause storage_beg / storage_end to disagree on
            ! the iteration range.
            !
            ! Phase 2 (time-varying forcing R, fractionation) will
            ! need a true urban mixed-box instead of the Phase-1
            ! re-sync.
            ! ---------------------------------------------------------------
            IF (DEF_USE_TRACER) THEN
               BLOCK
               USE MOD_Tracer_Defs, only: ntracers_u => ntracers, tracers_u => tracers, &
                  delta_to_R_u => delta_to_R
               USE MOD_Tracer_Vars, only: trc_rnof_step_u => trc_rnof_step, &
                                           a_trc_precip_u => a_trc_precip, &
                                           a_trc_evap_u   => a_trc_evap, &
                                           a_trc_rnof_u   => a_trc_rnof, &
                                           trc_wetwat, trc_waterstorage, &
                                           sync_tracer_patch_phase1
               USE MOD_Tracer_Conservation, only: tracer_save_storage, tracer_balance_check
               USE MOD_Tracer_Hist, only: tracer_hist_accumulate
               integer  :: itrc_u, j_u, snl_u
               real(r8) :: R_input_u, precip_mass_u, rnof_mass_u
               real(r8) :: evap_mass_u, dep_mass_u, xerr_tracer_u

               ! Reconstruct snl from current snow layer water content —
               ! CoLMMAIN_Urban maintains the aggregate arrays but
               ! does not surface snl. Same recipe used by tracer_rest
               ! init_from_water at MOD_Tracer_Rest.F90.
               snl_u = 0
               DO j_u = 0, maxsnl + 1, -1
                  IF (wliq_soisno(j_u, i) + wice_soisno(j_u, i) > 0._r8) THEN
                     snl_u = snl_u - 1
                  ELSE
                     EXIT
                  ENDIF
               ENDDO

               ! Purge pools that don't belong to urban (wetland and
               ! irrigation reservoir) BEFORE the storage snapshot —
               ! mirrors the glacier/lake pattern in CoLMMAIN. This
               ! prevents a phantom conservation spike when a LULCC
               ! class switch to urban leaves foreign wetwat or
               ! waterstorage mass in those slots.
               DO itrc_u = 1, ntracers_u
                  trc_wetwat(itrc_u, i) = 0._r8
                  IF (allocated(trc_waterstorage)) THEN
                     trc_waterstorage(itrc_u, i) = 0._r8
                  ENDIF
               ENDDO

               ! Snap pre-sync state; use `maxsnl` as the iteration
               ! lower bound so the snow-layer count can change
               ! freely inside CoLMMAIN_Urban without making
               ! storage_beg and storage_end disagree.
               CALL tracer_save_storage(i, maxsnl, nl_soil)

               CALL sync_tracer_patch_phase1(i, snl_u, maxsnl, nl_soil, &
                  wliq_soisno(:, i), wice_soisno(:, i), &
                  wa(i), wdsrf(i), scv(i), &
                  ldew_rain = ldew_rain(i), ldew_snow = ldew_snow(i))

               ! Atmospheric input / evap / runoff accumulators.
               ! fevpa>0: evap (loss); fevpa<0: condensation (gain).
               precip_mass_u = (forc_rain(i) + forc_snow(i)) * deltim
               rnof_mass_u   = max(rnof(i), 0._r8) * deltim
               evap_mass_u   = max( fevpa(i), 0._r8) * deltim
               dep_mass_u    = max(-fevpa(i), 0._r8) * deltim

               DO itrc_u = 1, ntracers_u
                  trc_rnof_step_u(itrc_u, i) = 0._r8
                  R_input_u = delta_to_R_u(tracers_u(itrc_u)%init_delta, &
                                           tracers_u(itrc_u)%ref_ratio)
                  IF (precip_mass_u > 0._r8) THEN
                     a_trc_precip_u(itrc_u, i) = a_trc_precip_u(itrc_u, i) &
                        + precip_mass_u * R_input_u
                  ENDIF
                  IF (dep_mass_u > 0._r8) THEN
                     a_trc_precip_u(itrc_u, i) = a_trc_precip_u(itrc_u, i) &
                        + dep_mass_u * R_input_u
                  ENDIF
                  IF (evap_mass_u > 0._r8) THEN
                     a_trc_evap_u(itrc_u, i) = a_trc_evap_u(itrc_u, i) &
                        + evap_mass_u * R_input_u
                  ENDIF
                  IF (rnof_mass_u > 0._r8) THEN
                     trc_rnof_step_u(itrc_u, i) = rnof_mass_u * R_input_u
                     a_trc_rnof_u(itrc_u, i) = a_trc_rnof_u(itrc_u, i) &
                        + rnof_mass_u * R_input_u
                  ENDIF
               ENDDO

               ! Close the per-step balance; same `maxsnl` iteration
               ! bound as save_storage.
               CALL tracer_balance_check(i, maxsnl, nl_soil, deltim, xerr_tracer_u, &
                  patchtype_in = 1)   ! urban

               ! Per-pool δ diagnostic. Urban has real canopy and
               ! soil columns, so pass ldew / wliq / wice / wa /
               ! wdsrf / scv. No wetland — pass 0 for wetwat.
               CALL tracer_hist_accumulate(i, snl_u, maxsnl, nl_soil, &
                  ldew_rain(i), ldew_snow(i), &
                  wliq_soisno(snl_u+1:nl_soil, i), &
                  wice_soisno(snl_u+1:nl_soil, i), &
                  wa(i), wdsrf(i), 0._r8, scv(i))
               END BLOCK
            ENDIF
         ENDIF

#endif
      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      ! Surface worst-patch tracer balance error accumulated across this
      ! step. Previously only ipatch==1 / itrc==1 was printed, hiding
      ! bulk failures elsewhere.
      IF (DEF_USE_TRACER) THEN
         BLOCK
         USE MOD_Tracer_Conservation, only: tracer_balance_report
         CALL tracer_balance_report()
         END BLOCK
      ENDIF

END SUBROUTINE CoLMDRIVER
! ---------- EOP ------------
