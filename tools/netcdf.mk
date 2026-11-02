# Only the GNU configuration opts into matching the installed Fortran dependency.
ifeq ($(NETCDF_AUTO),1)
ifneq ($(strip $(MAKECMDGOALS)),clean)
# GNU make otherwise remakes included files even under --dry-run.
ifeq ($(findstring n,$(firstword -$(MAKEFLAGS)))$(filter -n --just-print --dry-run --recon,$(MAKEFLAGS)),)
NETCDF_CACHE ?= .netcdf
NETCDF_PREPARE = env OMPI_FC='$(OMPI_FC)' MPICH_FC='$(MPICH_FC)' I_MPI_F90='$(I_MPI_F90)' I_MPI_FC='$(I_MPI_FC)' python3 tools/prepare_netcdf.py --fc '$(FF)' $(if $(filter YES,$(CaMa)),--check-fc '$(FCMP)') --include '$(NETCDF_INC)' --lib '$(NETCDF_LIB)' --cache-dir '$(NETCDF_CACHE)'
NETCDF_KEY := $(shell $(NETCDF_PREPARE) --key)
ifeq ($(strip $(NETCDF_KEY)),)
$(error Unable to identify NetCDF toolchain; see diagnostics above)
endif
NETCDF_CONFIG := $(NETCDF_CACHE)/$(NETCDF_KEY)/config.mk
.PHONY: check_netcdf
$(NETCDF_CONFIG): check_netcdf
	@$(NETCDF_PREPARE)
include $(NETCDF_CONFIG)
# The matching library must precede system -L paths, including user link flags.
override LDFLAGS := $(NETCDF_F_LDFLAGS) $(LDFLAGS)
endif
endif
endif
