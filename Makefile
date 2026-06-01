# Makefile for CoLM program

include include/Makeoptions
LINK_FOPTS ?= ${FOPTS}
HEADER = include/define.h
TRACER_ENABLED := $(shell printf '\043include "include/define.h"\n\043ifdef TRACER\nYES\n\043else\nNO\n\043endif\n' | cpp -P -I. -Iinclude - | awk 'NF { last=$$0 } END { print last }')

INCLUDE_DIR = -Iinclude -I.bld/ -I${NETCDF_INC}
MOD_OUT_DIR = .bld
MOD_OUT = $(MOD_CMD)$(if $(filter -J,$(strip $(MOD_CMD))),, )$(MOD_OUT_DIR)
VPATH = include : share : mksrfdata : mkinidata \
	: main : main/TRACER : main/HYDRO : main/BGC : main/URBAN : main/LULCC : main/DA \
	: main/ParaOpt : extends/CaMa/src : postprocess : .bld

# ********** Targets ALL **********
.PHONY: all
all : mkdir_build mksrfdata.x mkinidata.x colm.x postprocess.x lib
	@echo ''
	@echo '*******************************************************'
	@echo '*                                                     *'
	@echo '*        Making all CoLM programs successfully.       *'
	@echo '*                                                     *'
	@echo '*******************************************************'
# ******* End of Targets ALL ******

.PHONY: mkdir_build
mkdir_build :
	mkdir -p .bld

OBJS_SHARED =    \
				  MOD_Precision.o              \
				  MOD_SPMD_Task.o              \
				  MOD_Namelist.o               \
				  MOD_Vars_Global.o            \
				  MOD_Const_Physical.o         \
				  MOD_Const_LC.o               \
				  MOD_Utils.o                  \
				  MOD_IncompleteGamma.o        \
				  MOD_UserDefFun.o             \
				  MOD_TimeManager.o            \
				  MOD_Const_PFT.o              \
				  MOD_NetCDFSerial.o           \
				  MOD_Block.o                  \
				  MOD_Grid.o                   \
				  MOD_Pixel.o                  \
				  MOD_DataType.o               \
				  MOD_NetCDFPoint.o            \
				  MOD_NetCDFBlock.o            \
				  MOD_CatchmentDataReadin.o    \
				  MOD_5x5DataReadin.o          \
				  MOD_Mesh.o                   \
				  MOD_Pixelset.o               \
				  MOD_NetCDFVector.o           \
				  MOD_RangeCheck.o             \
				  MOD_SpatialMapping.o         \
				  MOD_WorkerPushData.o         \
				  MOD_AggregationRequestData.o \
				  MOD_PixelsetShared.o         \
				  MOD_LandElm.o                \
				  MOD_LandHRU.o                \
				  MOD_LandPatch.o              \
				  MOD_Land2mWMO.o              \
				  MOD_LandCrop.o               \
				  MOD_LandPFT.o                \
				  MOD_LandUrban.o              \
				  MOD_Urban_Const_LCZ.o        \
				  MOD_SingleSrfdata.o          \
				  MOD_SrfdataDiag.o            \
				  MOD_SrfdataRestart.o         \
				  MOD_ElmVector.o              \
				  MOD_HRUVector.o              \
				  MOD_MeshFilter.o             \
				  MOD_RegionClip.o

${OBJS_SHARED} : %.o : %.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o .bld/$@ $< $(MOD_OUT)

OBJS_SHARED_T = $(addprefix .bld/,${OBJS_SHARED})

OBJS_MKSRFDATA = \
				  Aggregation_PercentagesPFT.o      \
				  Aggregation_LAI.o                 \
				  Aggregation_SoilHyperAlbedo.o     \
					  Aggregation_SoilBrightness.o      \
					  Aggregation_LakeDepth.o           \
					  Aggregation_LakeSoilC.o          \
					  Aggregation_MethanePH.o          \
					  Aggregation_ForestHeight.o        \
				  Aggregation_SoilParameters.o      \
				  Aggregation_DBedrock.o            \
				  Aggregation_Topography.o          \
				  Aggregation_TopoWetness.o         \
				  Aggregation_TopographyFactors.o   \
				  Aggregation_TopographyFactors_Simple.o \
				  Aggregation_Urban.o               \
				  Aggregation_SoilTexture.o         \
				  MOD_Lulcc_TransferTrace.o         \
				  MKSRFDATA.o

$(OBJS_MKSRFDATA) : %.o : %.F90 ${HEADER} ${OBJS_SHARED}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o .bld/$@ $< $(MOD_OUT)

OBJS_MKSRFDATA_T = $(addprefix .bld/,${OBJS_MKSRFDATA})

# ------- Target 1: mksrfdata --------
mksrfdata.x : mkdir_build ${HEADER} ${OBJS_SHARED} ${OBJS_MKSRFDATA}
	@echo ''
	@echo 'making CoLM surface data start >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	@echo ''
	${FF} ${LINK_FOPTS} ${OBJS_SHARED_T} ${OBJS_MKSRFDATA_T} -o run/mksrfdata.x ${LDFLAGS}
	@echo ''
	@echo '<<<<<<<<<<<<<<<<<<<<<<<<<< making CoLM surface data completed!'
	@echo ''
# ----- End of Target 1 mksrfdata ----

TRACER_ISOTOPE_REGISTERED_SPECIES_OBJS = \
				 MOD_Tracer_Isotope_O18.o         \
				 MOD_Tracer_Isotope_HDO.o

TRACER_ISOTOPE_MAIN_OBJS = \
				 MOD_Tracer_Isotope_Registry.o    \
				 $(TRACER_ISOTOPE_REGISTERED_SPECIES_OBJS) \
				 MOD_Tracer_Isotope_Registrations.o

TRACER_REACTIVE_METHANE_OBJ_SUFFIXES = \
				 Const Registry GIEMS pH VegOverride State Microbes \
				 BgcLink AccFlux Physics Driver Hist Impl

TRACER_REACTIVE_METHANE_OBJS = \
				 $(addprefix MOD_Tracer_Reactive_Methane_, \
				   $(addsuffix .o,$(TRACER_REACTIVE_METHANE_OBJ_SUFFIXES))) \
				 MOD_Tracer_Reactive_Methane.o

TRACER_REACTIVE_REGISTERED_SPECIES_OBJS = \
				  $(TRACER_REACTIVE_METHANE_OBJS)

TRACER_REACTIVE_MAIN_OBJS = \
				MOD_Tracer_Reactive_BgcShim.o             \
				$(TRACER_REACTIVE_REGISTERED_SPECIES_OBJS) \
				MOD_Tracer_Reactive_Registrations.o
TRACER_PARTICLE_REGISTERED_SPECIES_OBJS = \
				MOD_Tracer_Particle_Sediment.o

TRACER_PARTICLE_DISPATCH_OBJS = \
				MOD_Tracer_Particle.o

TRACER_PARTICLE_MAIN_OBJS = \
				$(TRACER_PARTICLE_REGISTERED_SPECIES_OBJS) \
				MOD_Tracer_Particle_Registrations.o
ifeq (${TRACER_ENABLED},YES)
TRACER_BASIC_OBJS = \
					 MOD_Tracer_Defs.o              \
					 $(TRACER_PARTICLE_DISPATCH_OBJS) \
					 $(TRACER_ISOTOPE_MAIN_OBJS)    \
				 MOD_Tracer_Frac.o              \
				 MOD_Tracer_EvapLimit.o         \
				 MOD_Tracer_Vars.o              \
				 MOD_Tracer_RiverLake.o     \
				 MOD_Tracer_Conservation.o      \
				 MOD_Tracer_SoilInit.o          \
				 MOD_Tracer_ForcingInput.o      \
				 MOD_Tracer_Forcing.o           \
				 MOD_Tracer_Precip.o            \
				 MOD_Tracer_Evapo.o             \
				 MOD_Tracer_SoilWater.o         \
				 MOD_Tracer_Snow.o              \
				 MOD_Tracer_Reactive.o          \
				 MOD_Tracer_Rest.o

TRACER_MAIN_OBJS = \
				MOD_Tracer_Main.o                          \
				MOD_Tracer_Hist.o                         \
				$(TRACER_PARTICLE_MAIN_OBJS)              \
				$(TRACER_REACTIVE_MAIN_OBJS)              \
				MOD_Tracer_SpecialPatches.o

TRACER_MKINIDATA_OBJS = \
				MOD_Tracer_Particle_Registrations_Stubs.o \
				MOD_Tracer_Reactive_Registrations_Stubs.o

MOD_Tracer_Isotope_Registrations.o: include/tracer_isotope_species.inc
MOD_Tracer_Frac.o: MOD_Tracer_Isotope_Registry.o MOD_Tracer_Isotope_Registrations.o
MOD_Tracer_ForcingInput.o: MOD_Tracer_Defs.o
MOD_Tracer_Forcing.o: MOD_UserSpecifiedForcing.o MOD_Tracer_Defs.o MOD_Tracer_Vars.o \
					     MOD_Tracer_Isotope_Registry.o MOD_Tracer_Isotope_Registrations.o \
					     MOD_Tracer_ForcingInput.o
MOD_Tracer_Evapo.o MOD_Tracer_SoilWater.o: MOD_Tracer_EvapLimit.o
MOD_Vars_1DAccFluxes.o: MOD_Forcing.o MOD_FrictionVelocity.o MOD_TurbulenceLEddy.o \
					     MOD_Catch_Hist.o MOD_Tracer_Main.o
MOD_HistGridded.o: MOD_HistWriteBack.o MOD_Forcing.o MOD_Vars_1DAccFluxes.o
MOD_HistVector.o MOD_HistSingle.o: MOD_Vars_1DAccFluxes.o
MOD_Tracer_Hist.o: MOD_HistGridded.o MOD_HistVector.o MOD_HistSingle.o MOD_Vars_1DAccFluxes.o
MOD_Tracer_Reactive_Methane.o: MOD_Tracer_Reactive_Methane_Impl.o MOD_Tracer_Reactive_Methane_Hist.o \
					     MOD_Tracer_Reactive_Methane_AccFlux.o MOD_Tracer_Reactive_Methane_Microbes.o \
					     MOD_Tracer_Reactive_Methane_State.o MOD_Tracer_Reactive_Methane_GIEMS.o \
				     MOD_Tracer_Reactive_Methane_pH.o MOD_Tracer_Reactive_Methane_VegOverride.o
MOD_Tracer_Reactive_Methane_Impl.o: MOD_Tracer_Reactive_Methane_Driver.o \
				     MOD_Tracer_Reactive_Methane_State.o MOD_Tracer_Reactive_Methane_Registry.o \
				     MOD_Tracer_Reactive_Methane_Const.o MOD_Tracer_Reactive_Methane_BgcLink.o \
				     MOD_Tracer_Reactive_BgcShim.o MOD_Tracer_Conservation.o
MOD_Tracer_Reactive_Methane_Driver.o: MOD_Tracer_Reactive_Methane_Physics.o \
				     MOD_Tracer_Reactive_Methane_State.o MOD_Tracer_Reactive_Methane_BgcLink.o \
				     MOD_Tracer_Reactive_Methane_VegOverride.o MOD_Tracer_Reactive_Methane_Microbes.o
MOD_Tracer_Reactive_Methane_Physics.o: MOD_Tracer_Reactive_Methane_Const.o \
				     MOD_Tracer_Reactive_Methane_GIEMS.o MOD_Tracer_Reactive_Methane_BgcLink.o \
				     MOD_Tracer_Reactive_Methane_VegOverride.o MOD_Tracer_Reactive_Methane_State.o
MOD_Tracer_Reactive_Methane_BgcLink.o: MOD_Tracer_Reactive_Methane_Const.o \
				     MOD_Tracer_Reactive_Methane_pH.o MOD_Tracer_Reactive_Methane_VegOverride.o
MOD_Tracer_Reactive_Methane_AccFlux.o: MOD_Tracer_Reactive_Methane_Const.o \
					     MOD_Tracer_Reactive_Methane_Microbes.o MOD_Tracer_Reactive_Methane_BgcLink.o
MOD_Tracer_Reactive_Methane_Microbes.o: MOD_Tracer_Reactive_Methane_Const.o
MOD_Tracer_Reactive_Methane_Hist.o: MOD_Tracer_Reactive_Methane_BgcLink.o \
					     MOD_Tracer_Reactive_Methane_Registry.o MOD_Tracer_Reactive_Methane_AccFlux.o \
					     MOD_HistGridded.o
MOD_Tracer_Reactive_Registrations.o: include/tracer_reactive_species.inc \
				     $(TRACER_REACTIVE_REGISTERED_SPECIES_OBJS)
MOD_Tracer_Particle.o: MOD_Tracer_Defs.o
MOD_Tracer_Particle_Registrations.o: include/tracer_particle_species.inc \
				     $(TRACER_PARTICLE_REGISTERED_SPECIES_OBJS)
MOD_Tracer_Particle_Sediment.o: MOD_Tracer_Particle.o MOD_Tracer_Defs.o MOD_Grid_RiverLakeNetwork.o MOD_Vector_ReadWrite.o
MOD_Grid_RiverLakeFlow.o MOD_Grid_RiverLakeHist.o MOD_Grid_RiverLakeTimeVars.o: MOD_Tracer_Particle.o
MOD_Tracer_RiverLake.o MOD_Grid_RiverLakeFlow.o: MOD_Grid_RiverLakeTimeVars.o
else
TRACER_BASIC_OBJS =
TRACER_MAIN_OBJS =
TRACER_MKINIDATA_OBJS =
endif

OBJS_BASIC =    \
				 MOD_Vector_ReadWrite.o         \
				 MOD_dataSpec_PDB.o             \
				 MOD_tav_abs.o                  \
				 MOD_prospect_DB.o              \
				 MOD_Catch_BasinNetwork.o       \
				 MOD_Catch_Vars_TimeVariables.o \
				 MOD_Catch_Vars_1DFluxes.o      \
				 MOD_Grid_RiverLakeNetwork.o    \
				 MOD_Grid_Reservoir.o           \
				 MOD_Grid_RiverLakeLevee.o      \
				 MOD_Grid_RiverLakeBifurcation.o \
				 $(TRACER_BASIC_OBJS)           \
				 MOD_Grid_RiverLakeTimeVars.o   \
				 MOD_Qsadv.o                    \
				 MOD_UserSpecifiedForcing.o     \
				 MOD_BGC_Vars_1DFluxes.o        \
				 MOD_BGC_Vars_1DPFTFluxes.o     \
				 MOD_BGC_Vars_PFTimeVariables.o \
				 MOD_BGC_Vars_TimeInvariants.o  \
				 MOD_BGC_Vars_TimeVariables.o   \
				 MOD_Urban_Vars_1DFluxes.o      \
				 MOD_Urban_Vars_TimeVariables.o \
				 MOD_Urban_Vars_TimeInvariants.o\
				 MOD_DA_Vars_1DFluxes.o         \
				 MOD_Vars_TimeInvariants.o      \
				 MOD_DA_Vars_TimeVariables.o    \
				 MOD_Vars_TimeVariables.o       \
				 MOD_Vars_1DPFTFluxes.o         \
				 MOD_Vars_1DFluxes.o            \
				 MOD_Vars_1DForcing.o           \
				 MOD_Hydro_SoilFunction.o       \
				 MOD_Hydro_SoilWater.o          \
				 MOD_Eroot.o                    \
				 MOD_LAIEmpirical.o             \
				 MOD_LAIReadin.o                \
				 MOD_CropReadin.o               \
				 MOD_NitrifData.o               \
				 MOD_NdepData.o                 \
				 MOD_FireData.o                 \
				 MOD_OrbCoszen.o                \
				 MOD_OrbCosazi.o                \
				 MOD_HighRes_Parameters.o		\
				 MOD_3DCanopyRadiation.o        \
				 MOD_Aerosol.o                  \
				 MOD_SnowSnicar.o               \
				 MOD_Albedo.o                   \
				 MOD_SnowSnicar_HiRes.o			\
				 MOD_Albedo_HiRes.o				\
				 MOD_SnowFraction.o             \
				 MOD_Urban_LAIReadin.o          \
				 MOD_Urban_Shortwave.o          \
				 MOD_Urban_Albedo.o             \
				 MOD_MonthlyinSituCO2MaunaLoa.o \
				 MOD_PercentagesPFTReadin.o     \
				 MOD_LakeDepthReadin.o          \
				 MOD_DBedrockReadin.o           \
				 MOD_SoilColorRefl.o            \
				 MOD_SoilParametersReadin.o     \
				 MOD_SoilTextureReadin.o        \
				 MOD_HtopReadin.o               \
				 MOD_UrbanReadin.o              \
				 MOD_BGC_CNSummary.o            \
				 MOD_IniTimeVariable.o          \
				 MOD_UrbanIniTimeVariable.o     \
				 MOD_ElementNeighbour.o         \
				 MOD_Catch_HillslopeNetwork.o   \
				 MOD_Catch_RiverLakeNetwork.o   \
				 MOD_Catch_Reservoir.o          \
				 MOD_VicParaReadin.o            \
				 MOD_Initialize.o

MOD_UserSpecifiedForcing.o: MOD_Qsadv.o

$(OBJS_BASIC) : %.o : %.F90 ${HEADER} ${OBJS_SHARED}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o .bld/$@ $< $(MOD_OUT)

OBJS_BASIC_T = $(addprefix .bld/,${OBJS_BASIC})

OBJS_MKINIDATA = \
				  $(TRACER_MKINIDATA_OBJS) \
				  CoLMINI.o

$(OBJS_MKINIDATA) : %.o : %.F90 ${HEADER} ${OBJS_SHARED} ${OBJS_BASIC}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o .bld/$@ $< $(MOD_OUT)

OBJS_MKINIDATA_T = $(addprefix .bld/,${OBJS_MKINIDATA})

# -------- Target 2: mkinidata -------
mkinidata.x : mkdir_build ${HEADER} ${OBJS_SHARED} ${OBJS_BASIC} ${OBJS_MKINIDATA}
	@echo ''
	@echo 'making CoLM initial data start >>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	@echo ''
	${FF} ${LINK_FOPTS} ${OBJS_SHARED_T} ${OBJS_BASIC_T} ${OBJS_MKINIDATA_T} -o run/mkinidata.x ${LDFLAGS}
	@echo ''
	@echo '<<<<<<<<<<<<<<<<<<<<<<<<< making CoLM initial data completed!'
	@echo ''
# ----- End of Target 2 mkinidata ----

CaMa = $(shell \
	CAMA_INIT=$$(grep -E '^\#(define|undef)[[:space:]]+CaMa_Flood' $(HEADER) | head -1 | grep -q 'define' && echo 1 || echo 0); \
	SP=$$(grep -E '^\#(define|undef)[[:space:]]+SinglePoint' $(HEADER) | tail -1 | grep -q 'define' && echo 1 || echo 0); \
	MP=$$(if [ $$SP -eq 1 ]; then echo 0; else grep -E '^\#(define|undef)[[:space:]]+USEMPI' $(HEADER) | head -1 | grep -q 'define' && echo 1 || echo 0; fi); \
	if [ $$CAMA_INIT -eq 0 ]; then echo NO; elif [ $$SP -eq 1 ] || [ $$MP -eq 0 ]; then echo NO; else echo YES; fi)
ifeq (${CaMa},YES)# Compile CoLM decoupled with river routing scheme (CaMa-Flood)
#
OBJECTS_CAMA=\
				  parkind1.o              \
				  yos_cmf_input.o         \
				  yos_cmf_time.o          \
				  yos_cmf_map.o           \
				  yos_cmf_prog.o          \
				  yos_cmf_diag.o          \
				  cmf_utils_mod.o         \
				  cmf_calc_outflw_mod.o   \
				  cmf_calc_pthout_mod.o   \
				  cmf_calc_fldstg_mod.o   \
				  cmf_calc_stonxt_mod.o   \
				  cmf_opt_outflw_mod.o    \
				  cmf_ctrl_tracer_mod.o   \
				  cmf_ctrl_mpi_mod.o      \
				  cmf_ctrl_damout_mod.o   \
				  cmf_ctrl_levee_mod.o    \
				  cmf_ctrl_forcing_mod.o  \
				  cmf_ctrl_boundary_mod.o \
				  cmf_ctrl_output_mod.o   \
				  cmf_ctrl_restart_mod.o  \
				  cmf_ctrl_sed_mod.o      \
				  cmf_calc_diag_mod.o     \
				  cmf_ctrl_physics_mod.o  \
				  cmf_ctrl_time_mod.o     \
				  cmf_ctrl_maps_mod.o     \
				  cmf_ctrl_vars_mod.o     \
				  cmf_ctrl_nmlist_mod.o   \
				  cmf_drv_control_mod.o   \
				  cmf_drv_advance_mod.o

$(OBJECTS_CAMA) : %.o : %.F90 ${HEADER}
	$(FCMP)  -c ${FFLAGS} $(MODS) ${CFLAGS} $(INCLUDE_DIR) -o .bld/$@ $< $(MOD_OUT)

OBJS_CAMA_T = $(addprefix .bld/,${OBJECTS_CAMA})

endif

OBJS_MAIN = \
				MOD_Catch_HillslopeFlow.o                 \
				MOD_Catch_SubsurfaceFlow.o                \
				MOD_Catch_RiverLakeFlow.o                 \
				MOD_Catch_Hist.o                          \
				MOD_Catch_WriteParameters.o               \
				MOD_BGC_CNCStateUpdate1.o                 \
				MOD_BGC_CNCStateUpdate2.o                 \
				MOD_BGC_CNCStateUpdate3.o                 \
				MOD_BGC_CNNStateUpdate1.o                 \
				MOD_BGC_CNNStateUpdate2.o                 \
				MOD_BGC_CNNStateUpdate3.o                 \
				MOD_BGC_Soil_BiogeochemNStateUpdate1.o    \
				MOD_BGC_Soil_BiogeochemNitrifDenitrif.o   \
				MOD_BGC_Soil_BiogeochemCompetition.o      \
				MOD_BGC_Soil_BiogeochemDecompCascadeBGC.o \
				MOD_BGC_Soil_BiogeochemDecomp.o           \
				MOD_BGC_Soil_BiogeochemLittVertTransp.o   \
				MOD_BGC_Soil_BiogeochemNLeaching.o        \
				MOD_BGC_Soil_BiogeochemPotential.o        \
				MOD_BGC_Soil_BiogeochemVerticalProfile.o  \
				MOD_BGC_Veg_CNGapMortality.o              \
				MOD_BGC_Veg_CNGResp.o                     \
				MOD_BGC_Veg_CNMResp.o                     \
				MOD_BGC_Daylength.o                       \
				MOD_BGC_Veg_CNPhenology.o                 \
				MOD_BGC_Veg_NutrientCompetition.o         \
				MOD_BGC_Veg_CNVegStructUpdate.o           \
				MOD_BGC_CNAnnualUpdate.o                  \
				MOD_BGC_CNZeroFluxes.o                    \
				MOD_BGC_CNBalanceCheck.o                  \
				MOD_BGC_CNSASU.o                          \
				MOD_BGC_Veg_CNNDynamics.o                 \
				MOD_BGC_Veg_CNFireBase.o                  \
				MOD_BGC_Veg_CNFireLi2016.o                \
				MOD_Vars_2DForcing.o                      \
				MOD_ForcingDownscaling.o                  \
				MOD_Forcing.o                             \
				MOD_DA_TWS.o                              \
				MOD_DA_Const.o                            \
				MOD_DA_RTM.o                              \
				MOD_DA_EnKF.o                             \
				MOD_DA_SM.o                               \
				MOD_DA_Ensemble.o                         \
				MOD_DA_Main.o                             \
				MOD_Opt_Baseflow.o                        \
				MOD_ParameterOptimization.o               \
				MOD_AssimStomataConductance.o             \
				MOD_PlantHydraulic.o                      \
				MOD_FrictionVelocity.o                    \
				MOD_TurbulenceLEddy.o                     \
				MOD_Ozone.o                               \
				MOD_CanopyLayerProfile.o                  \
				MOD_LeafTemperature.o                     \
				MOD_LeafTemperaturePC.o                   \
				MOD_SoilThermalParameters.o               \
				MOD_Hydro_VIC_Variables.o                 \
				MOD_Hydro_VIC.o                           \
				MOD_Runoff.o                              \
				MOD_SoilSnowHydrology.o                   \
				MOD_SnowLayersCombineDivide.o             \
				MOD_PhaseChange.o                         \
				MOD_Glacier.o                             \
				MOD_Lake.o                                \
				MOD_SimpleOcean.o                         \
				MOD_GroundFluxes.o                        \
				MOD_GroundTemperature.o                   \
				MOD_LeafInterception.o                    \
				MOD_NetSolar.o                            \
				MOD_NetSolar_Hyper.o                      \
				MOD_WetBulb.o                             \
				MOD_RainSnowTemp.o                        \
				MOD_SoilSurfaceResistance.o               \
				MOD_NewSnow.o                             \
				MOD_Thermal.o                             \
				$(TRACER_MAIN_OBJS)                       \
				MOD_Vars_1DAccFluxes.o                    \
				MOD_CaMa_Vars.o                           \
				MOD_Irrigation.o                          \
				MOD_BGC_driver.o                          \
				MOD_HistWriteBack.o                       \
				MOD_HistGridded.o                         \
				MOD_HistVector.o                          \
				MOD_HistSingle.o                          \
				MOD_Grid_RiverLakeHist.o                  \
					MOD_Hist.o                                \
					MOD_CheckEquilibrium.o                    \
				MOD_LightningData.o                       \
				MOD_CaMa_colmCaMa.o                       \
				MOD_Catch_LateralFlow.o                   \
				MOD_Grid_RiverLakeFlow.o                  \
				MOD_Urban_Longwave.o                      \
				MOD_Urban_NetSolar.o                      \
				MOD_Urban_Flux.o                          \
				MOD_Urban_GroundFlux.o                    \
				MOD_Urban_RoofFlux.o                      \
				MOD_Urban_RoofTemperature.o               \
				MOD_Urban_WallTemperature.o               \
				MOD_Urban_PerviousTemperature.o           \
				MOD_Urban_ImperviousTemperature.o         \
				MOD_Urban_Hydrology.o                     \
				MOD_Urban_BEM.o                           \
				MOD_Urban_LUCY.o                          \
				MOD_Urban_Thermal.o                       \
				CoLMMAIN_Urban.o                          \
				MOD_Lulcc_Vars_TimeInvariants.o           \
				MOD_Lulcc_Vars_TimeVariables.o            \
				MOD_Lulcc_TransferTraceReadin.o           \
				MOD_Lulcc_MassEnergyConserve.o            \
				MOD_Lulcc_Initialize.o                    \
				MOD_Lulcc_Driver.o                        \
				CoLMDRIVER.o                              \
				CoLMMAIN.o                                \
				CoLM.o

$(OBJS_MAIN) : %.o : %.F90 ${HEADER} ${OBJS_SHARED} ${OBJS_BASIC}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o .bld/$@ $< $(MOD_OUT)




OBJS_MAIN_T = $(addprefix .bld/,${OBJS_MAIN})

# ------ Target 3: main --------

ifneq (${CaMa},YES)# Compile CoLM decoupled without river routing scheme (CaMa-Flood)

colm.x : mkdir_build ${HEADER} ${OBJS_SHARED} ${OBJS_BASIC} ${OBJS_MAIN}
	@echo ''
	@echo 'making CoLM start >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	@echo ''
	${FF} ${LINK_FOPTS} ${OBJS_SHARED_T} ${OBJS_BASIC_T} ${OBJS_MAIN_T} -o run/colm.x ${LDFLAGS}
	@echo ''
	@echo '<<<<<<<<<<<<<<<<<<<<<<<<<<<<< making CoLM completed!'
	@echo ''

else
colm.x : mkdir_build  ${HEADER} ${OBJS_SHARED} ${OBJECTS_CAMA} ${OBJS_BASIC} ${OBJS_MAIN}
	@echo ''
	@echo 'making CoLM with CaMa start >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	@echo ''
	${FF} ${LINK_FOPTS} ${OBJS_SHARED_T} ${OBJS_BASIC_T} ${OBJS_CAMA_T} ${OBJS_MAIN_T} -o run/colm.x ${LDFLAGS}

	@echo ''
	@echo '<<<<<<<<<<<<<<<<<<<<<<<<<<<< making CoLM with CaMa completed!'
	@echo ''

endif

# ----- End of Target 3 main -----

OBJS_POST1 = MOD_Concatenate.o HistConcatenate.o
OBJS_POST2 = MOD_Vector2Grid.o POST_Vector2Grid.o
OBJS_POST3 = SrfDataConcatenate.o
OBJS_POST1_T = $(addprefix .bld/,${OBJS_POST1})
OBJS_POST2_T = $(addprefix .bld/,${OBJS_POST2})
OBJS_POST3_T = $(addprefix .bld/,${OBJS_POST3})

$(OBJS_POST1):%.o:%.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o .bld/$@ $< $(MOD_OUT)

$(OBJS_POST2):%.o:%.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o .bld/$@ $< $(MOD_OUT)

$(OBJS_POST3):%.o:%.F90 ${HEADER}
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o .bld/$@ $< $(MOD_OUT)

hist_concatenate.x : ${HEADER} ${OBJS_SHARED} ${OBJS_POST1}
	${FF} ${LINK_FOPTS} ${OBJS_SHARED_T} ${OBJS_POST1_T} -o run/$@ ${LDFLAGS}

post_vector2grid.x : ${HEADER} ${OBJS_SHARED} ${OBJS_POST2}
	${FF} ${LINK_FOPTS} ${OBJS_SHARED_T} ${OBJS_POST2_T} -o run/$@ ${LDFLAGS}

srfdata_concatenate.x : ${HEADER} ${OBJS_SHARED} ${OBJS_POST3}
	${FF} ${LINK_FOPTS} ${OBJS_SHARED_T} ${OBJS_POST3_T} -o run/$@ ${LDFLAGS}

# ------ Target 4: postprocess --------
DEF = $(shell grep -i CATCHMENT include/define.h)
vector2grid = $(word 1, ${DEF})
ifneq (${vector2grid},\#define)
DEF = $(shell grep -i UNSTRUCTURED include/define.h)
vector2grid = $(word 1, ${DEF})
endif

.PHONY: postprocess.x
ifneq (${vector2grid},\#define)
postprocess.x : mkdir_build hist_concatenate.x srfdata_concatenate.x
	@echo '<<<<<<<<<<<<<<<<<<<<<<<<< making CoLM postprocessing completed!'
else
postprocess.x : mkdir_build hist_concatenate.x srfdata_concatenate.x post_vector2grid.x
	@echo '<<<<<<<<<<<<<<<<<<<<<<<<< making CoLM postprocessing completed!'
endif
# --- End of Target 4 postprocess ------

# ------ Target 5: static libs --------
.PHONY: lib
lib : mksrfdata.x mkinidata.x colm.x postprocess.x
	@echo ''
	@echo 'making CoLM static library >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	mkdir -p lib
	cd lib && find ../.bld -name "*.o" ! -name "CoLM.o" ! -name "MKSRFDATA.o" ! -name "CoLMINI.o" -exec ln -sf {} ./ \;
	cd lib && ar rc libcolm.a *.o && ranlib libcolm.a
	ln -sf lib/libcolm.a ./libcolm.a
# ------End of Target 5: static libs --------

.PHONY: clean
clean :
	rm -rf .bld
	rm -rf lib libcolm.a
	rm -f run/mksrfdata.x run/mkinidata.x run/colm.x
	rm -f run/hist_concatenate.x run/srfdata_concatenate.x run/post_vector2grid.x
	rm -f CaMa/src/*.o CaMa/src/*.mod CaMa/src/*.a
