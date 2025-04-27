#include <define.h>

MODULE MOD_Isotope_Flx
#ifdef USE_ISOTOPE

!-----------------------------------------------------------------------
   USE MOD_Precision
   IMPLICIT NONE
   SAVE



   real, parameter      :: rd = 287.04, sigma = 5.67e-8,    &
   cph2o = 4.218e+3,cpice = 2.106e+3,            &
   lsubf = 3.335e+5

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: isotope_canopy_transpiration_flx
   PUBLIC :: isotope_canopy_evaporation_flx
   PUBLIC :: isotope_soil_evaporation_flx

CONTAINS

!-----------------------------------------------------------------------

   subroutine caltmp(t1, sfctmp, sfcprs, zlvl, q2, th2, t1v, th2v, rho ) 

      implicit none

      ! input:
      real, intent(in)       :: t1     ! skin temperature (k)
      real, intent(in)       :: sfctmp ! air temperature (k) at level zlvl
      real, intent(in)       :: q2     ! specific humidity (kg/kg) at level zlvl
      real, intent(in)       :: sfcprs ! atmospheric pressure (pa) at level zlvl
      real, intent(in)       :: zlvl   ! height (m agl) where atmospheric fields are valid

      ! output:
      real, intent(out)      :: th2    ! potential temperature (referenced to the surface pressure)
      real, intent(out)      :: t1v    ! virtual skin temperature (k)
      real, intent(out)      :: th2v   ! virtual potential temperature at zlvl
      real, intent(out)      :: rho    ! density

      ! local:
      real                   :: t2v

      th2 = sfctmp + ( 0.0098 * zlvl)
      t1v= t1 * (1. + 0.61*q2) 
      th2v = th2 * (1. + 0.61*q2)
      t2v = sfctmp * ( 1. + 0.61*q2 )
      rho = sfcprs / (rd * t2v)

   end subroutine caltmp

   ! -------------------------------------------------------------------------------

   subroutine calhum(sfctmp, sfcprs, q2sat, dqsdt2)

      implicit none
   
      ! input:
      real, intent(in)       :: sfctmp
      real, intent(in)       :: sfcprs
      
      ! output:
      real, intent(out)      :: q2sat   ! saturated specific humidity
      real, intent(out)      :: dqsdt2 
      
      ! local
      real, parameter        :: a2=17.67,a3=273.15,a4=29.65, elwv=2.501e6,        &
                              a23m4=a2*(a3-a4), e0=611.0, rv=461.0,             &
                              epsilon=0.622 
      real                   :: es

      ! es:  e.g. dutton chapter 8, eq 11
         es = e0 * exp ( elwv/rv*(1./a3 - 1./sfctmp) )

      ! q2sat:  
         q2sat = epsilon * es / (sfcprs - (1.-epsilon)*es)

      ! dqsdt2 is calculated assuming q2sat is a specific humidity

         dqsdt2=q2sat*a23m4/(sfctmp-a4)**2 

   end subroutine calhum


   subroutine metalib(e_a,P_a,T_a,Cp,Psy,delta,e_as)

      real(8), intent(in)       :: e_a
      real(8), intent(in)       :: P_a
      real(8), intent(in)       :: T_a
      real(8), intent(out)       :: Cp
      real(8), intent(out)       :: Psy !psychrometic constant
      real(8), intent(out)       :: delta!slope of the saturated wapor-temperature curve
      real(8), intent(out)       :: e_as
      real(8) ::T_a1
      T_a1=T_a-273.15
      e_as=0.6108*exp(17.27*T_a1/(T_a1 + 237.3))!%Kpa
      Cp = 0.24 * 4185.5 * (1.0 + 0.8 * (0.622 * e_a/ (P_a - e_a)))!;%  [J kg-1 ?C-1],
      Psy=0.000665*P_a
      delta=4098.0*(0.6108*exp(17.27*T_a1/(T_a1+237.3)))/ ((T_a1+237.3)**2.0) !kPaC-1
      return
   end subroutine metalib

   subroutine aerodynamicresistance(Zh,wsp,Rav,Rag)

      real(8), intent(in)       :: Zh
      real(8), intent(in)       :: wsp
      real(8), intent(out)      :: Rav
      real(8), intent(out)      :: Rag
      real(8) ::  Zm
      real(8) ::  d_0
      real(8) ::  Zth
      real(8) ::  Z0mv
      real(8) ::  Z0hv
      real(8) ::  Z0mg
      real(8) ::  Z0hg
      real(8) ::  k
      Zm=10.0 !Zm the height of wind speed measurement (m)
      d_0=0.666*Zh!
      Zth=10.0 !Zth the height of temperature and humidity measurement (m)
      Z0mv=0.123*Zh!Z0mv is the roughness length govering momentum transfer above vegetation canopy (m)
      Z0hv=0.1*Z0mv!Z0hv is the roughness length govering heat and vapor transfer above vegetation canopy (m)
      Z0mg=0.001!Z0mg is the roughness length govering momentum transfer above ground surface(m)
      Z0hg=1.0*Z0mg!Z0hg is the roughness length govering heat and vapor transfer above ground surface(m)
      k=0.41!k von Karman constant
      Rav=(log((Zm-d_0)/Z0mv))*(log((Zth-d_0)/Z0hv))/(k*k*wsp)!Rav the aerodynamic resistance for vegetation canopy(s m-1)
      Rag=(log((Zm)/Z0mg))*(log((Zth)/Z0hg))/(k*k*wsp)!Rag the aerodynamic resistance for ground surface(s m-1)
   end subroutine aerodynamicresistance

   subroutine canopy_temperature(T_a,e_a,e_as,rsc,Rav, &
                  LEc,Cp,Psy,rho_a,T_sk)
 
       real(8), intent(in)       :: T_a
       real(8), intent(in)       :: e_a
       real(8), intent(in)       :: e_as
       real(8), intent(in)       :: rsc
       real(8), intent(in)       :: Rav
       real(8), intent(in)       :: LEc
       real(8), intent(in)       :: Cp
       real(8), intent(in)       :: Psy
       real(8), intent(in)       :: rho_a
       real(8), intent(out)      :: T_sk
       real(8) :: Ta_1
       real(8) :: diff_eas
       real(8) :: gw
       real(8) :: DD,VPD
       Ta_1=T_a-273.16
       VPD=e_as-e_a
       diff_eas=-(1527.0*exp((1727.0*Ta_1)/(100.0*(Ta_1+2373.0/10.0))) &
                   *((1727.0*Ta_1)/(100.*(Ta_1 + 2373./10.)**2)  &
                   - 1727.0/(100.0*(Ta_1 + 2373./10.0))))/2500.0
       gw=1.0/(Rav+rsc)
       DD=LEc*Psy/(Cp*rho_a*gw)
       T_sk=(DD-VPD)/diff_eas+T_a
 
       return
       end subroutine canopy_temperature

   SUBROUTINE isotope_canopy_transpiration_flx(z0mv,rc,LEc,rho_a,wsp,T_a, T_c, &
            q2,P_a1,LAI,Vapor_O,Vapor_D,soil_O,soil_D,Xylem_O,Xylem_D, &
            delta_TO, delta_TD,Rhl)
      real(8), intent(in)       :: z0mv !roughness length for momentum transfer above vegetation canopy (m)
      real(8), intent(in)       :: rc !canopy resistance (s m-1)
      real(8), intent(in)       :: LEc !canopy transpiration (W m-2)
      real(8), intent(in)       :: rho_a !air density (kg m-3)
      real(8), intent(in)       :: wsp !wind speed (m s-1)
      real(8), intent(in)       :: T_a !air temperature (K)
      real(8), intent(inout)       :: q2 !specific humidity (kg kg-1)
      real(8), intent(in)       :: P_a1 !air pressure (Pa)
      real(8), intent(in)       :: LAI !leaf area index
      real(8), intent(in)       :: Vapor_O, Vapor_D  !isotopic composition of water vapor (‰)
      real(8), intent(in)       :: T_c !canopy temperature (K)
      real(8), intent(inout)    :: soil_O, soil_D  !isotopic composition of soil water (‰)
      real(8), intent(inout)    :: Xylem_O, Xylem_D  !isotopic composition of xylem water (‰)
      real(8), intent(out)      :: delta_TO, delta_TD !isotopic composition of transpiration (‰)
      real(8), intent(out)      :: Rhl 

      real(8),  parameter :: Rsmow_O = 2005.2 / 1000000.
      real(8),  parameter :: Rsmow_D = 155.76  / 1000000.
      real(8),  parameter :: ddi_O = 1.03189
      real(8),  parameter :: ddi_D = 1.01636
      real(8)  :: e_a1
      real(8)  :: Rawv_O,Rawv_D,Rsxw_O,Rsxw_D,Rssw_O,Rssw_D
      real(8)  :: wcleaf4,fes,lep,lepen,lepns,lepss
      real(8)  :: Rlwes_O, Rblw_D
      real(8)  :: Zh, d0, z0hv, z0mg, z0hg, cumdew, rhow
      real(8)  :: Cp,Psy,delta,e_as
      real(8)  :: rav,rag,T_sk
      real(8)  :: uv,uca,uc,lw
      real(8)  :: ra,rb,e_a, P_a
      real(8)  ::alphaeql_O,alphaeql_D
      real(8)  ::dlil_O,dlil_D
      real(8)  ::alphakl_O,alphakl_D
      real(8)  ::pn_O,pn_D,esattl,qsattl
      real(8)  ::Rtrf1_O, Rtrf1_D,Rlwes1_O,Rlwes1_D,Rblw1_O,Rblw1_D
      real(8)  ::Rlwes4_O,Rlwes4_D,Rblw4_O,Rblw4_D,Rtrf4_O,Rtrf4_D
      real(8)  ::Revf_O, Revf_D,lhe,RH,dqsdt2,q2sat
      integer  :: j
      real(8)  ::Aeq_O,Aeq_D
      Xylem_O=-10.0
      Xylem_D=-70.0
      soil_O=-10.0
      soil_D=-70.0

      Rawv_O = (Vapor_O / 1000. + 1.0) * Rsmow_O
      Rawv_D = (Vapor_D / 1000. + 1.0) * Rsmow_D
      Rsxw_O = (Xylem_O / 1000. + 1.0) * Rsmow_O
      Rsxw_D = (Xylem_D / 1000. + 1.0) * Rsmow_D
      Rssw_O = (soil_O  / 1000. + 1.0) * Rsmow_O
      Rssw_D = (soil_D  / 1000. + 1.0) * Rsmow_D
      wcleaf4=0.333
      fes=0.9
      lep=0.2
      lepen=0.2
      lepns=0.2
      lepss=0.2
      Rlwes_O = Rsxw_O
      Rblw_D  = Rsxw_D
      Zh=z0mv/0.123
      d0 = 0.666 * Zh
      z0hv = 0.1 * z0mv
      z0mg = 1.0 / 100.
      z0hg = 1.0 * z0mg
      cumdew = 0.0
      call calhum(T_a,P_a1,                  q2sat,dqsdt2)     ! out
      if (q2 < 0.1e-5) q2 = 0.1e-5
      if (q2 >=  q2sat) q2 = q2sat*0.99
      e_a1 =q2 * P_a1 / (0.622+0.378*q2)
      e_a=e_a1/1000.0
      P_a=P_a1/1000.0
      rhow = 1000.
      lhe = 2.50025 * 1000000. - 2.365 * 1000. * (t_a-273.15)
      call metalib(e_a,P_a,T_a,Cp,Psy,delta,e_as)
      RH=e_a/e_as
      call aerodynamicresistance(Zh,wsp,rav,rag)
      uv = wsp * log((Zh - d0) / z0mv) / log((10.0 - d0) / z0mv) !Baldocchi et al., 1988
      uca = -0.03 * LAI * LAI + 0.66 * LAI + 0.7
      uc = uv * (1.0 - exp(-uca)) / uca
      lw =0.02
      rb = 283.0 * ( (lw / uc) ** 0.5) / 2.0 / LAI
      ra = rav - rb

      alphaeql_O=  exp(1137.0 / T_c / T_c - 0.4156 / T_c - 0.002067)! canopy_temperature K
      alphaeql_D = exp(24840.0 / T_c / T_c - 76.25 / T_c + 0.0526)
      dlil_O = 119.0 * exp(-637.0 / (T_c - 137.0)) / 1000000000.0
      dlil_D = 116.0 * exp(-626.0 / (T_c - 139.0)) / 1000000000.0

      alphakl_O = (ra + rb * (ddi_O ** 0.666) + ddi_O * rc)/ &
                  (ra + rb + rc)
      alphakl_D = (ra + rb * (ddi_D ** 0.666) + ddi_D * rc)/ &
                  (ra + rb + rc)


      pn_O = LEc/lhe * lep / rhow / dlil_O
      pn_D = LEc/lhe * lep / rhow / dlil_D
      esattl = 6.112 * Exp(17.67 * (T_c - 273.15) / &
                            (T_c -273.15 + 243.5))/10.
      qsattl = 0.622 * esattl / (P_a - 0.378 * esattl)


      Rhl=e_a / esattl



      Rlwes1_O = alphaeql_O * (alphakl_O * Rsxw_O * &
                            (1.0 - e_a / esattl) + Rawv_O * e_a / esattl)
      Rlwes1_D = alphaeql_D * (alphakl_D * Rsxw_D * &
                            (1.0 - e_a / esattl) + Rawv_D * e_a / esattl)
      Rblw1_O = Rlwes1_O
      Rblw1_D = Rlwes1_D

      Rlwes4_O= Rsxw_O
      Rlwes4_D= Rsxw_D
      Rblw4_O = Rsxw_O
      Rblw4_D = Rsxw_D

      wcleaf4=0.333

      do j=1,20
      Rlwes4_O = Rlwes4_O + 1.0 * qsattl * pn_O * &
                      (Rlwes1_O - Rlwes4_O) / &
                      (alphakl_O * alphaeql_O * (rav + rc) * &
                      wcleaf4 * (1.0 - Exp(-1.0 * pn_O)))
      Rlwes4_D = Rlwes4_D + 1.0 * qsattl * pn_D * &
                      (Rlwes1_D - Rlwes4_D) / &
                      (alphakl_D * alphaeql_D * (rav + rc) * &
                      wcleaf4 * (1.0 - Exp(-1.0 * pn_D)))
      Rblw4_O = Rblw4_O + 1.0 * qsattl * pn_O * &
                      (Rblw1_O - Rblw4_O) / &
                      (alphakl_O * alphaeql_O * (rav + rc) * &
                      wcleaf4 * (1.0 - Exp(-1.0 * pn_O)))
      Rblw4_D = Rblw4_D + 1.0 * qsattl * pn_D * &
                      (Rblw1_D - Rblw4_D) / &
                      (alphakl_D * alphaeql_D * (rav + rc) * &
                      wcleaf4 * (1.0 - Exp(-1.0 * pn_D)))
      enddo
      Rblw4_O = (1.0 - Exp(-1.0 * pn_O)) * &
                (Rlwes4_O - Rsmow_O) / pn_O + Rsmow_O
      Rblw4_D = (1.0 - Exp(-1.0 * pn_D)) * &
                (Rlwes4_D - Rsmow_D) / pn_D + Rsmow_D
      Rtrf4_O  = Rsxw_O - (Rlwes4_O  - Rlwes1_O ) / &
                      (alphaeql_O  * alphakl_O  * (1.0 - e_a / esattl))
      Rtrf4_D  = Rsxw_D - (Rlwes4_D - Rlwes1_D) / &
                      (alphaeql_D * alphakl_D * (1.0 - e_a / esattl))

       delta_TO  = (Rtrf4_O/Rsmow_O - 1.0) * 1000.0
       delta_TD  = (Rtrf4_D/Rsmow_D - 1.0) * 1000.0



      if (isnan(delta_TO)) then
        delta_TO=-30.0
        delta_TD=-230.0
      endif

      if (isnan(delta_TD)) then
        delta_TO=-30.0
        delta_TD=-230.0
      endif

   END SUBROUTINE isotope_canopy_transpiration_flx

      
   SUBROUTINE isotope_soil_evaporation_flx(dt,z0mv,rc,rss,LEc,LEs,rho_a,wsp,T_a, &
      q2,P_a1,LAI,Vapor_O,Vapor_D,soil_O,soil_D,Xylem_O,Xylem_D,T_g, &
      delta_TO, delta_TD, delta_SO, delta_SD,delta_ET_O,delta_ET_D,Rhs,Rhl)

      real(8), intent(in)       :: dt !time step (s)
      real(8), intent(in)       :: z0mv !roughness length for momentum transfer above vegetation canopy (m)
      real(8), intent(in)       :: rc !canopy resistance (s m-1)
      real(8), intent(in)       :: rss !soil surface resistance (s m-1)
      real(8), intent(in)       :: LEc !canopy transpiration (W m-2)
      real(8), intent(in)       :: LEs !soil evaporation (W m-2)
      real(8), intent(in)       :: rho_a !air density (kg m-3)
      real(8), intent(in)       :: wsp !wind speed (m s-1)
      real(8), intent(in)       :: T_a !air temperature (K)
      real(8), intent(inout)       :: q2 !specific humidity (kg kg-1)
      real(8), intent(in)       :: P_a1 !air pressure (Pa)
      real(8), intent(in)       :: LAI !leaf area index
      real(8), intent(in)       :: Vapor_O, Vapor_D  !isotopic composition of water vapor (‰)
      real(8), intent(inout)    :: soil_O, soil_D  !isotopic composition of soil water (‰)
      real(8), intent(inout)    :: Xylem_O, Xylem_D  !isotopic composition of xylem water (‰)
      real(8), intent(in)       :: T_g !ground temperature (K)
      real(8), intent(out)      :: delta_TO, delta_TD !isotopic composition of transpiration (‰)
      real(8), intent(out)      :: delta_SO, delta_SD !isotopic composition of soil evaporation (‰)
      real(8), intent(out)      :: delta_ET_O, delta_ET_D !isotopic composition of transpiration (‰)
      real(8), intent(out)      :: Rhl !isotopic composition of soil evaporation (‰)
      real(8), intent(out)      :: Rhs !isotopic composition of transpiration (‰)

      real(8),  parameter :: Rsmow_O = 2005.2 / 1000000.
      real(8),  parameter :: Rsmow_D = 155.76  / 1000000.
      real(8),  parameter :: ddi_O = 1.03189
      real(8),  parameter :: ddi_D = 1.01636
      real(8)  :: e_a1
      real(8)  ::Rawv_O,Rawv_D,Rsxw_O,Rsxw_D,Rssw_O,Rssw_D
      real(8)  ::wcleaf4,fes,lep,lepen,lepns,lepss
      real(8)  ::Rlwes_O, Rblw_D
      real(8)  ::Zh, d0, z0hv, z0mg, z0hg, cumdew, rhow
      real(8)  ::Cp,Psy,delta,e_as
      real(8)  ::rav,rag,T_sk
      real(8)  ::uv,uca,uc,lw
      real(8)  ::T_c
      real(8)  ::ra,rb,e_a, P_a
      real(8)  ::alphaeql_O,alphaeql_D,alphaeqg_O,alphaeqg_D
      real(8)  ::dlil_O,dlil_D,dlig_O,dlig_D
      real(8)  ::alphakl_O,alphakl_D,alphakg_O,alphakg_D
      real(8)  ::pn_O,pn_D,esattl,qsattl,esattg,qsattg
      real(8)  ::Rtrf1_O, Rtrf1_D,Rlwes1_O,Rlwes1_D,Rblw1_O,Rblw1_D
      real(8)  ::Rlwes4_O,Rlwes4_D,Rblw4_O,Rblw4_D,Rtrf4_O,Rtrf4_D
      real(8)  ::Revf_O, Revf_D,lhe,RH,dqsdt2,q2sat
      integer  :: j
      real(8)  ::Aeq_O,Aeq_D
      Xylem_O=-10.0
      Xylem_D=-70.0
      soil_O=-10.0
      soil_D=-70.0

      Rawv_O = (Vapor_O / 1000. + 1.0) * Rsmow_O
      Rawv_D = (Vapor_D / 1000. + 1.0) * Rsmow_D
      Rsxw_O = (Xylem_O / 1000. + 1.0) * Rsmow_O
      Rsxw_D = (Xylem_D / 1000. + 1.0) * Rsmow_D
      Rssw_O = (soil_O  / 1000. + 1.0) * Rsmow_O
      Rssw_D = (soil_D  / 1000. + 1.0) * Rsmow_D
      wcleaf4=0.333
      fes=0.9
      lep=0.2
      lepen=0.2
      lepns=0.2
      lepss=0.2
      Rlwes_O = Rsxw_O
      Rblw_D  = Rsxw_D
      Zh=z0mv/0.123
      d0 = 0.666 * Zh
      z0hv = 0.1 * z0mv
      z0mg = 1.0 / 100.
      z0hg = 1.0 * z0mg
      cumdew = 0.0
      call calhum(T_a,P_a1,                  q2sat,dqsdt2)     ! out
      if (q2 < 0.1e-5) q2 = 0.1e-5
      if (q2 >=  q2sat) q2 = q2sat*0.99
      e_a1 =q2 * P_a1 / (0.622+0.378*q2)
      e_a=e_a1/1000.0
      P_a=P_a1/1000.0
      rhow = 1000.
      lhe = 2.50025 * 1000000. - 2.365 * 1000. * (t_a-273.15)
      call metalib(e_a,P_a,T_a,Cp,Psy,delta,e_as)
      RH=e_a/e_as
      call aerodynamicresistance(Zh,wsp,rav,rag)
      call canopy_temperature(T_a,e_a,e_as,rc,rav,                  LEc,Cp,Psy,rho_a,T_c)
      uv = wsp * log((Zh - d0) / z0mv) / log((10.0 - d0) / z0mv) !Baldocchi et al., 1988
      uca = -0.03 * LAI * LAI + 0.66 * LAI + 0.7
      uc = uv * (1.0 - exp(-uca)) / uca
      lw =0.02
      rb = 283.0 * ( (lw / uc) ** 0.5) / 2.0 / LAI
      ra = rav - rb

      alphaeql_O=  exp(1137.0 / T_c / T_c - 0.4156 / T_c - 0.002067)! canopy_temperature K
      alphaeql_D = exp(24840.0 / T_c / T_c - 76.25 / T_c + 0.0526)
      alphaeqg_O = exp(1137.0 / T_g / T_g - 0.4156 / T_g - 0.002067)
      alphaeqg_D = exp(24840.0 / T_g / T_g - 76.25 / T_g + 0.0526)!! soil_temperature K

      dlil_O = 119.0 * exp(-637.0 / (T_c - 137.0)) / 1000000000.0
      dlil_D = 116.0 * exp(-626.0 / (T_c - 139.0)) / 1000000000.0
      dlig_O = 119.0 * exp(-637.0 / (T_g - 137.0)) / 1000000000.0
      dlig_D = 116.0 * exp(-626.0 / (T_g - 139.0)) / 1000000000.0

      alphakl_O = (ra + rb * (ddi_O ** 0.666) + ddi_O * rc)/ (ra + rb + rc)
      alphakl_D = (ra + rb * (ddi_D ** 0.666) + ddi_D * rc)/  (ra + rb + rc)
      alphakg_O = (rag + ddi_O * rss) / (rag + rss)
      alphakg_D = (rag + ddi_D * rss) / (rag + rss)

      pn_O = LEc/lhe * lep / rhow / dlil_O
      pn_D = LEc/lhe * lep / rhow / dlil_D
      esattl = 6.112 * Exp(17.67 * (T_c - 273.15) / (T_c -273.15 + 243.5))/10.
      qsattl = 0.622 * esattl / (P_a - 0.378 * esattl)

      esattg = 6.112 * Exp(17.67 * (T_g - 273.15) / (T_g - 273.15 + 243.5))/10.
      qsattg = 0.622 * esattg / (P_a - 0.378 * esattg)
      Rhs=e_a / esattg
      Rhl=e_a / esattl



      Rlwes1_O = alphaeql_O * (alphakl_O * Rsxw_O * (1.0 - e_a / esattl) + Rawv_O * e_a / esattl)
      Rlwes1_D = alphaeql_D * (alphakl_D * Rsxw_D * (1.0 - e_a / esattl) + Rawv_D * e_a / esattl)
      Rblw1_O = Rlwes1_O
      Rblw1_D = Rlwes1_D

      Rlwes4_O= Rsxw_O
      Rlwes4_D= Rsxw_D
      Rblw4_O = Rsxw_O
      Rblw4_D = Rsxw_D

      wcleaf4=0.333


      do j=1,200
      Rlwes4_O = Rlwes4_O + 1.0 * qsattl * pn_O *  &
                      (Rlwes1_O - Rlwes4_O) /  &
                      (alphakl_O * alphaeql_O * (rav + rc) &
                      * wcleaf4 * (1.0 - Exp(-1.0 * pn_O)))
      Rlwes4_D = Rlwes4_D + 1.0 * qsattl * pn_D * &
                      (Rlwes1_D - Rlwes4_D) / &
                      (alphakl_D * alphaeql_D * (rav + rc) &
                      * wcleaf4 * (1.0 - Exp(-1.0 * pn_D)))
      Rblw4_O = Rblw4_O + 1.0 * qsattl * pn_O * &
                      (Rblw1_O - Rblw4_O) / &
                      (alphakl_O * alphaeql_O * (rav + rc) &
                      * wcleaf4 * (1.0 - Exp(-1.0 * pn_O)))
      Rblw4_D = Rblw4_D + 1.0 * qsattl * pn_D * &
                      (Rblw1_D - Rblw4_D) / &
                      (alphakl_D * alphaeql_D * (rav + rc) &
                      * wcleaf4 * (1.0 - Exp(-1.0 * pn_D)))
      enddo
      Rblw4_O = (1.0 - Exp(-1.0 * pn_O)) * (Rlwes4_O - Rsmow_O) / pn_O + Rsmow_O
      Rblw4_D = (1.0 - Exp(-1.0 * pn_D)) * (Rlwes4_D - Rsmow_D) / pn_D + Rsmow_D
      Rtrf4_O  = Rsxw_O - (Rlwes4_O  - Rlwes1_O ) / &
                      (alphaeql_O  * alphakl_O  * (1.0 - e_a / esattl))
      Rtrf4_D  = Rsxw_D - (Rlwes4_D - Rlwes1_D) / &
                      (alphaeql_D * alphakl_D * (1.0 - e_a / esattl))

      delta_TO  = (Rtrf4_O/Rsmow_O - 1.0) * 1000.0
      delta_TD  = (Rtrf4_D/Rsmow_D - 1.0) * 1000.0

      Aeq_O = (1.0 - 1.0/alphaeqg_O)*1000.0
      Aeq_D = (1.0 - 1.0/alphaeqg_D)*1000.0

      delta_SO=((1.0/alphaeqg_O) * soil_O - Rhs * Vapor_O - &
                      Aeq_O - (1.0 - Rhs )* 28.0) / &
                      ((1.0-Rhs)+ 0.001*((1.0 - Rhs) * 28.0))

      delta_SD=((1.0/alphaeqg_D) * soil_D - Rhs * Vapor_D - &
                      Aeq_D - (1.0 - Rhs )* 25.0) / &
                      ((1.0-Rhs)+ 0.001*((1.0 - Rhs) * 25.0))

      if (isnan(delta_TO)) then
      delta_TO=-30.0
      delta_TD=-230.0
      endif

      if (isnan(delta_TD)) then
      delta_TO=-30.0
      delta_TD=-230.0
      endif

      if (isnan(delta_SO)) then
      delta_SO=-30.0
      delta_SD=-230.0
      endif
      if (isnan(delta_SD)) then
      delta_SO=-30.0
      delta_SD=-230.0
      endif


      delta_ET_O = (delta_SO*LEs+delta_TO*LEc)/(LEs+LEc)
      delta_ET_D = (delta_SD*LEs+delta_TD*LEc)/(LEs+LEc)


   END SUBROUTINE isotope_soil_evaporation_flx

   SUBROUTINE isotope_canopy_evaporation_flx(dt,z0mv,rc,rss,LEc,LEs,rho_a,wsp,T_a, &
      q2_ref,P_a1,LAI,Vapor_O,Vapor_D,soil_O,soil_D,Xylem_O,Xylem_D,T_g, &
      delta_TO, delta_TD, delta_SO, delta_SD,delta_ET_O,delta_ET_D,Rhs,Rhl)

      real(8), intent(in)       :: dt !time step (s)
      real(8), intent(in)       :: z0mv !roughness length for momentum transfer above vegetation canopy (m)
      real(8), intent(in)       :: rc !canopy resistance (s m-1)
      real(8), intent(in)       :: rss !soil surface resistance (s m-1)
      real(8), intent(in)       :: LEc !canopy transpiration (W m-2)
      real(8), intent(in)       :: LEs !soil evaporation (W m-2)
      real(8), intent(in)       :: rho_a !air density (kg m-3)
      real(8), intent(in)       :: wsp !wind speed (m s-1)
      real(8), intent(in)       :: T_a !air temperature (K)
      real(8), intent(inout)       :: q2_ref !specific humidity (kg kg-1)
      real(8), intent(in)       :: P_a1 !air pressure (Pa)
      real(8), intent(in)       :: LAI !leaf area index
      real(8), intent(in)       :: Vapor_O, Vapor_D  !isotopic composition of water vapor (‰)
      real(8), intent(inout)    :: soil_O, soil_D  !isotopic composition of soil water (‰)
      real(8), intent(inout)    :: Xylem_O, Xylem_D  !isotopic composition of xylem water (‰)
      real(8), intent(in)       :: T_g !ground temperature (K)
      real(8), intent(out)      :: delta_TO, delta_TD !isotopic composition of transpiration (‰)
      real(8), intent(out)      :: delta_SO, delta_SD !isotopic composition of soil evaporation (‰)
      real(8), intent(out)      :: delta_ET_O, delta_ET_D !isotopic composition of transpiration (‰)
      real(8), intent(out)      :: Rhl !isotopic composition of soil evaporation (‰)
      real(8), intent(out)      :: Rhs !isotopic composition of transpiration (‰)

      real(8),  parameter :: Rsmow_O = 2005.2 / 1000000.
      real(8),  parameter :: Rsmow_D = 155.76  / 1000000.
      real(8),  parameter :: ddi_O = 1.03189
      real(8),  parameter :: ddi_D = 1.01636
      real(8)  :: e_a1,q2
      real(8)  ::Rawv_O,Rawv_D,Rsxw_O,Rsxw_D,Rssw_O,Rssw_D
      real(8)  ::wcleaf4,fes,lep,lepen,lepns,lepss
      real(8)  ::Rlwes_O, Rblw_D
      real(8)  ::Zh, d0, z0hv, z0mg, z0hg, cumdew, rhow
      real(8)  ::Cp,Psy,delta,e_as
      real(8)  ::rav,rag,T_sk
      real(8)  ::uv,uca,uc,lw
      real(8)  ::T_c
      real(8)  ::ra,rb,e_a, P_a
      real(8)  ::alphaeql_O,alphaeql_D,alphaeqg_O,alphaeqg_D
      real(8)  ::dlil_O,dlil_D,dlig_O,dlig_D
      real(8)  ::alphakl_O,alphakl_D,alphakg_O,alphakg_D
      real(8)  ::pn_O,pn_D,esattl,qsattl,esattg,qsattg
      real(8)  ::Rtrf1_O, Rtrf1_D,Rlwes1_O,Rlwes1_D,Rblw1_O,Rblw1_D
      real(8)  ::Rlwes4_O,Rlwes4_D,Rblw4_O,Rblw4_D,Rtrf4_O,Rtrf4_D
      real(8)  ::Revf_O, Revf_D,lhe,RH,dqsdt2,q2sat
      integer  :: j
      real(8)  ::Aeq_O,Aeq_D
      Xylem_O=-10.0
      Xylem_D=-70.0
      soil_O=-10.0
      soil_D=-70.0

      Rawv_O = (Vapor_O / 1000. + 1.0) * Rsmow_O
      Rawv_D = (Vapor_D / 1000. + 1.0) * Rsmow_D
      Rsxw_O = (Xylem_O / 1000. + 1.0) * Rsmow_O
      Rsxw_D = (Xylem_D / 1000. + 1.0) * Rsmow_D
      Rssw_O = (soil_O  / 1000. + 1.0) * Rsmow_O
      Rssw_D = (soil_D  / 1000. + 1.0) * Rsmow_D
      wcleaf4=0.333
      fes=0.9
      lep=0.2
      lepen=0.2
      lepns=0.2
      lepss=0.2
      Rlwes_O = Rsxw_O
      Rblw_D  = Rsxw_D
      Zh=z0mv/0.123
      d0 = 0.666 * Zh
      z0hv = 0.1 * z0mv
      z0mg = 1.0 / 100.
      z0hg = 1.0 * z0mg
      cumdew = 0.0
      call calhum(T_a,P_a1,                  q2sat,dqsdt2)     ! out
      q2=q2_ref
      if (q2 < 0.1e-5) q2 = 0.1e-5
      if (q2 >=  q2sat) q2 = q2sat*0.99
      e_a1 =q2 * P_a1 / (0.622+0.378*q2)
      e_a=e_a1/1000.0
      P_a=P_a1/1000.0
      rhow = 1000.
      lhe = 2.50025 * 1000000. - 2.365 * 1000. * (t_a-273.15)
      call metalib(e_a,P_a,T_a,Cp,Psy,delta,e_as)
      RH=e_a/e_as
      call aerodynamicresistance(Zh,wsp,rav,rag)
      call canopy_temperature(T_a,e_a,e_as,rc,rav,               LEc,Cp,Psy,rho_a,T_c)
      uv = wsp * log((Zh - d0) / z0mv) / log((10.0 - d0) / z0mv) !Baldocchi et al., 1988
      uca = -0.03 * LAI * LAI + 0.66 * LAI + 0.7
      uc = uv * (1.0 - exp(-uca)) / uca
      lw =0.02
      rb = 283.0 * ( (lw / uc) ** 0.5) / 2.0 / LAI
      ra = rav - rb

      alphaeql_O=  exp(1137.0 / T_c / T_c - 0.4156 / T_c - 0.002067)! canopy_temperature K
      alphaeql_D = exp(24840.0 / T_c / T_c - 76.25 / T_c + 0.0526)
      alphaeqg_O = exp(1137.0 / T_g / T_g - 0.4156 / T_g - 0.002067)
      alphaeqg_D = exp(24840.0 / T_g / T_g - 76.25 / T_g + 0.0526)!! soil_temperature K

      dlil_O = 119.0 * exp(-637.0 / (T_c - 137.0)) / 1000000000.0
      dlil_D = 116.0 * exp(-626.0 / (T_c - 139.0)) / 1000000000.0
      dlig_O = 119.0 * exp(-637.0 / (T_g - 137.0)) / 1000000000.0
      dlig_D = 116.0 * exp(-626.0 / (T_g - 139.0)) / 1000000000.0

      alphakl_O = (ra + rb * (ddi_O ** 0.666) + ddi_O * rc)/ &
                      (ra + rb + rc)
      alphakl_D = (ra + rb * (ddi_D ** 0.666) + ddi_D * rc)/ &
                      (ra + rb + rc)
      alphakg_O = (rag + ddi_O * rss) / (rag + rss)
      alphakg_D = (rag + ddi_D * rss) / (rag + rss)

      pn_O = LEc/lhe * lep / rhow / dlil_O
      pn_D = LEc/lhe * lep / rhow / dlil_D
      esattl = 6.112 * Exp(17.67 * (T_c - 273.15) / &
                      (T_c -273.15 + 243.5))/10.
      qsattl = 0.622 * esattl / (P_a - 0.378 * esattl)

      esattg = 6.112 * Exp(17.67 * (T_g - 273.15) / &
                      (T_g - 273.15 + 243.5))/10.
      qsattg = 0.622 * esattg / (P_a - 0.378 * esattg)
      Rhs=e_a / esattg
      Rhl=e_a / esattl



      Rlwes1_O = alphaeql_O * (alphakl_O * Rsxw_O * &
                      (1.0 - e_a / esattl) + Rawv_O * e_a / esattl)
      Rlwes1_D = alphaeql_D * (alphakl_D * Rsxw_D * &
                      (1.0 - e_a / esattl) + Rawv_D * e_a / esattl)
      Rblw1_O = Rlwes1_O
      Rblw1_D = Rlwes1_D

      Rlwes4_O= Rsxw_O
      Rlwes4_D= Rsxw_D
      Rblw4_O = Rsxw_O
      Rblw4_D = Rsxw_D

      wcleaf4=0.333


      do j=1,200
      Rlwes4_O = Rlwes4_O + 1.0 * qsattl * pn_O * &
                      (Rlwes1_O - Rlwes4_O) / &
                      (alphakl_O * alphaeql_O * (rav + rc) &
                      * wcleaf4 * (1.0 - Exp(-1.0 * pn_O)))
      Rlwes4_D = Rlwes4_D + 1.0 * qsattl * pn_D * &
                      (Rlwes1_D - Rlwes4_D) / &
                      (alphakl_D * alphaeql_D * (rav + rc) &
                      * wcleaf4 * (1.0 - Exp(-1.0 * pn_D)))
      Rblw4_O = Rblw4_O + 1.0 * qsattl * pn_O * &
                      (Rblw1_O - Rblw4_O)/ &
                      (alphakl_O * alphaeql_O * (rav + rc) &
                      * wcleaf4 * (1.0 - Exp(-1.0 * pn_O)))
      Rblw4_D = Rblw4_D + 1.0 * qsattl * pn_D * &
                      (Rblw1_D - Rblw4_D)/ &
                      (alphakl_D * alphaeql_D * (rav + rc) &
                      * wcleaf4 * (1.0 - Exp(-1.0 * pn_D)))
      enddo
      Rblw4_O = (1.0 - Exp(-1.0 * pn_O)) * (Rlwes4_O - Rsmow_O) / pn_O + Rsmow_O  * (Rlwes4_O - Rsmow_O) / pn_O + Rsmow_O
      Rblw4_D = (1.0 - Exp(-1.0 * pn_D))     * (Rlwes4_D - Rsmow_D) / pn_D + Rsmow_D
      Rtrf4_O  = Rsxw_O - (Rlwes4_O  - Rlwes1_O ) /     (alphaeql_O  * alphakl_O  * (1.0 - e_a / esattl))
      Rtrf4_D  = Rsxw_D - (Rlwes4_D - Rlwes1_D) /     (alphaeql_D * alphakl_D * (1.0 - e_a / esattl))

      delta_TO  = (Rtrf4_O/Rsmow_O - 1.0) * 1000.0
      delta_TD  = (Rtrf4_D/Rsmow_D - 1.0) * 1000.0

      Aeq_O = (1.0 - 1.0/alphaeqg_O)*1000.0
      Aeq_D = (1.0 - 1.0/alphaeqg_D)*1000.0

      delta_SO=((1.0/alphaeqg_O) * soil_O - Rhs * Vapor_O - Aeq_O - (1.0 - Rhs )* 28.0) / ((1.0-Rhs)+ 0.001*((1.0 - Rhs) * 28.0))

      delta_SD=((1.0/alphaeqg_D) * soil_D - Rhs * Vapor_D -      Aeq_D - (1.0 - Rhs )* 25.0) /     ((1.0-Rhs)+ 0.001*((1.0 - Rhs) * 25.0))

      if (isnan(delta_TO)) then
      delta_TO=-30.0
      delta_TD=-230.0
      endif

      if (isnan(delta_TD)) then
      delta_TO=-30.0
      delta_TD=-230.0
      endif

      if (isnan(delta_SO)) then
      delta_SO=-30.0
      delta_SD=-230.0
      endif
      if (isnan(delta_SD)) then
      delta_SO=-30.0
      delta_SD=-230.0
      endif


      delta_ET_O = (delta_SO*LEs+delta_TO*LEc)/(LEs+LEc)
      delta_ET_D = (delta_SD*LEs+delta_TD*LEc)/(LEs+LEc)


   END SUBROUTINE isotope_canopy_evaporation_flx 

#endif
END MODULE MOD_Isotope_Flx
