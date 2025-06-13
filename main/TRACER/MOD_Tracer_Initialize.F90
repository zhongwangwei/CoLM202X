#include <define.h>

MODULE MOD_Tracer_Initialize

USE MOD_Precision, only: r8
USE MOD_Tracer_Namelist_Defs, only: DEF_Tracers, DEF_Tracer_Number, DEF_Tracer_Init, DEF_Tracer_Init_file
USE MOD_Tracer_State, only: tracer_soil_concentration, allocate_tracer_state, actual_num_patches_alloc, actual_num_soil_layers_alloc
USE MOD_SPMD_Task, only: p_is_master, numpatch, MPI_REAL8, p_comm_glb, p_err, p_address_master, mpi_bcast, CoLM_stop
USE MOD_NetCDFSerial, only: ncio_open, ncio_inq_varid, ncio_close, ncio_file_exist, ncio_inq_var_ndims, ncio_inq_var_dimid, ncio_inq_dimlen, ncio_inq_dimname, ncio_get_var_real8_1d, ncio_get_var_real8_2d, ncio_get_var_real8_3d, NC_MAX_DIMS, NC_MAX_NAME
USE MOD_Vars_Soil, only: nl_soil
USE MOD_LandPatch, only: landpatch, landpatch_vector_type ! For mapping
USE MOD_Grid, only: grid_type                             ! For defining grid from init file
USE MOD_SpatialMapping, only: spatial_mapping_type        ! For grid2pset mapping
USE MOD_DataType, only: allocate_block_data, flush_block_data ! For temp_grid_layer_data

IMPLICIT NONE
SAVE

PUBLIC :: tracer_initialize_concentrations

CONTAINS

SUBROUTINE tracer_initialize_concentrations()
  integer :: i_tracer, i_patch, i_layer, i
  integer :: nc_file_id, nc_var_id, nc_err
  character(len=256) :: var_in_file_name

  ! Variables for reading NetCDF structure
  integer :: ndims_in_file, dimids_in_file(NC_MAX_DIMS)
  character(len=NC_MAX_NAME) :: dimname_in_file
  integer :: dimlen_in_file
  integer :: file_num_layers, file_num_lat, file_num_lon, file_num_patch
  integer :: layer_dim_idx, lat_dim_idx, lon_dim_idx, patch_dim_idx

  real(r8), allocatable :: temp_data_3d(:,:,:), temp_data_2d(:,:)
  real(r8), allocatable :: file_lat_coords(:), file_lon_coords(:)
  real(r8), allocatable :: mapped_patch_data_1layer(:)

  type(grid_type) :: init_file_grid
  type(spatial_mapping_type) :: init_file_mapper
  logical :: init_file_grid_defined = .false.
  logical :: init_file_mapper_built = .false.

  IF (.NOT. allocated(DEF_Tracers) .OR. DEF_Tracer_Number == 0) THEN
     IF (p_is_master) WRITE(*,*) "Tracer Initialize: DEF_Tracers not allocated or no tracers defined."
     RETURN
  ENDIF

  IF (numpatch <= 0 .OR. nl_soil <= 0) THEN
      IF (p_is_master) WRITE(*,*) "Tracer Initialize: numpatch or nl_soil is zero. Cannot proceed."
      ! Ensure tracer_soil_concentration is not allocated or is deallocated if it was by a previous call
      CALL allocate_tracer_state(0,0,0) ! This will deallocate if already allocated
      RETURN
  ENDIF
  
  CALL allocate_tracer_state(DEF_Tracer_Number, numpatch, nl_soil)
  ! Default initialization to 0.0_r8 is done within allocate_tracer_state

  IF (DEF_Tracer_Init) THEN
     IF (p_is_master) THEN
        IF (TRIM(DEF_Tracer_Init_file) == 'null' .OR. .NOT. ncio_file_exist(TRIM(DEF_Tracer_Init_file))) THEN
           WRITE(*,*) "Tracer Init File: '", TRIM(DEF_Tracer_Init_file), "' is 'null' or does not exist. Using default (zero) concentrations."
        ELSE
           nc_err = ncio_open(TRIM(DEF_Tracer_Init_file), 'READ', nc_file_id)
           IF (nc_err /= 0) THEN
              WRITE(*,*) "Error opening tracer init file: ", TRIM(DEF_Tracer_Init_file)
           ELSE
              DO i_tracer = 1, DEF_Tracer_Number
                 var_in_file_name = TRIM(DEF_Tracers(i_tracer)%name) // "_init"
                 nc_var_id = ncio_inq_varid(nc_file_id, TRIM(var_in_file_name))

                 IF (nc_var_id == -1) THEN
                    WRITE(*,*) "Warning: Var ", TRIM(var_in_file_name), " not found in ", TRIM(DEF_Tracer_Init_file), ". Using default for this tracer."
                    tracer_soil_concentration(i_tracer, :, :) = 0.0_r8
                    CYCLE
                 ENDIF

                 ndims_in_file = ncio_inq_var_ndims(nc_file_id, nc_var_id)
                 CALL ncio_inq_var_dimid(nc_file_id, nc_var_id, dimids_in_file)

                 file_num_layers = -1; file_num_lat = -1; file_num_lon = -1; file_num_patch = -1
                 layer_dim_idx = -1; lat_dim_idx = -1; lon_dim_idx = -1; patch_dim_idx = -1

                 DO i = 1, ndims_in_file
                    dimname_in_file = ncio_inq_dimname(nc_file_id, dimids_in_file(i))
                    dimlen_in_file  = ncio_inq_dimlen(nc_file_id, dimids_in_file(i))
                    SELECT CASE (TRIM(ADJUSTL(dimname_in_file)))
                       CASE ("layer", "soil_layer", "lev", "soil_lev", "soil_layers", "soil_level")
                          file_num_layers = dimlen_in_file; layer_dim_idx = i
                       CASE ("lat", "latitude", "lats", "latitudes")
                          file_num_lat = dimlen_in_file; lat_dim_idx = i
                       CASE ("lon", "longitude", "lons", "longitudes")
                          file_num_lon = dimlen_in_file; lon_dim_idx = i
                       CASE ("patch", "patches", "element", "elements", "pft", "landpatch") 
                          file_num_patch = dimlen_in_file; patch_dim_idx = i
                    END SELECT
                 END DO

                 ! --- Case A: (layer, lat, lon) ---
                 IF (ndims_in_file == 3 .AND. layer_dim_idx /= -1 .AND. lat_dim_idx /= -1 .AND. lon_dim_idx /= -1) THEN
                    WRITE(*,*) "Info: Reading ", TRIM(var_in_file_name), " as (layer,lat,lon)"
                    IF (.NOT. init_file_grid_defined) THEN
                       ALLOCATE(file_lat_coords(file_num_lat))
                       ALLOCATE(file_lon_coords(file_num_lon))
                       CALL ncio_get_var_real8_1d(nc_file_id, ncio_inq_varid(nc_file_id, 'lat'), file_lat_coords) 
                       CALL ncio_get_var_real8_1d(nc_file_id, ncio_inq_varid(nc_file_id, 'lon'), file_lon_coords)
                       CALL init_file_grid%define_by_center(file_lat_coords, file_lon_coords)
                       DEALLOCATE(file_lat_coords, file_lon_coords)
                       init_file_grid_defined = .TRUE.
                       CALL init_file_mapper%build_arealweighted(init_file_grid, landpatch)
                       init_file_mapper_built = .TRUE.
                    ENDIF

                    IF (.NOT. init_file_mapper_built) THEN
                       WRITE(*,*) "Error: Failed to build mapper for init file grid. Skipping ", TRIM(var_in_file_name)
                       tracer_soil_concentration(i_tracer, :, :) = 0.0_r8
                       CYCLE
                    ENDIF
                    
                    ALLOCATE(temp_data_3d(file_num_lon, file_num_lat, file_num_layers)) ! Assuming (lon,lat,layer) order in file
                    CALL ncio_get_var_real8_3d(nc_file_id, nc_var_id, temp_data_3d) 

                    ALLOCATE(mapped_patch_data_1layer(numpatch))
                    DO i_layer = 1, MIN(file_num_layers, nl_soil)
                        type(block_data_real8_2d) :: temp_grid_layer_data
                        CALL allocate_block_data(init_file_grid, temp_grid_layer_data)
                        
                        IF (init_file_grid%xcnt(1) == file_num_lon .AND. init_file_grid%ycnt(1) == file_num_lat) THEN
                             temp_grid_layer_data%blk(1,1)%val(:,:) = temp_data_3d(:,:,i_layer) ! Assuming (lon,lat,layer)
                        ELSE
                             WRITE(*,*) "Grid dimension mismatch for init data layer copy. Skipping layer for: ", TRIM(var_in_file_name)
                             tracer_soil_concentration(i_tracer, :, i_layer) = 0.0_r8
                             CYCLE
                        ENDIF

                        CALL init_file_mapper%grid2pset(temp_grid_layer_data, mapped_patch_data_1layer)
                        tracer_soil_concentration(i_tracer, :, i_layer) = mapped_patch_data_1layer(:)
                        CALL flush_block_data(temp_grid_layer_data, 0.0_r8)
                        ! DEALLOCATE for temp_grid_layer_data is handled by its finalizer if it's a full type
                    END DO
                    DEALLOCATE(temp_data_3d, mapped_patch_data_1layer)
                    IF (file_num_layers /= nl_soil) THEN
                       WRITE(*,*) "Warning: Mismatch soil layers for ", TRIM(var_in_file_name), ". File: ", file_num_layers, ", Model: ", nl_soil
                    ENDIF

                 ! --- Case B: (layer, patch) ---
                 ELSE IF (ndims_in_file == 2 .AND. layer_dim_idx /= -1 .AND. patch_dim_idx /= -1) THEN
                    WRITE(*,*) "Info: Reading ", TRIM(var_in_file_name), " as (layer,patch)"
                    IF (file_num_patch /= numpatch) THEN
                       WRITE(*,*) "Warning: Patch num mismatch for ", TRIM(var_in_file_name), ". File: ", file_num_patch, ", Model: ", numpatch, " Skipping."
                       tracer_soil_concentration(i_tracer, :, :) = 0.0_r8
                       CYCLE
                    ENDIF
                    ALLOCATE(temp_data_2d(file_num_layers, file_num_patch)) ! Or (patch, layer) depending on file
                    CALL ncio_get_var_real8_2d(nc_file_id, nc_var_id, temp_data_2d) 
                    
                    ! Determine if file is (layer,patch) or (patch,layer) based on dimids and known conventions
                    ! For simplicity, assume (layer,patch) if layer_dim_idx < patch_dim_idx (Fortran-like)
                    ! or (patch,layer) if patch_dim_idx < layer_dim_idx.
                    IF (dimids_in_file(layer_dim_idx) < dimids_in_file(patch_dim_idx)) THEN ! (layer, patch)
                        DO i_layer = 1, MIN(file_num_layers, nl_soil)
                           tracer_soil_concentration(i_tracer, :, i_layer) = temp_data_2d(i_layer, :)
                        END DO
                    ELSE ! (patch, layer)
                        DO i_layer = 1, MIN(file_num_layers, nl_soil)
                           tracer_soil_concentration(i_tracer, :, i_layer) = temp_data_2d(:, i_layer)
                        END DO
                    ENDIF
                    DEALLOCATE(temp_data_2d)
                    IF (file_num_layers /= nl_soil) THEN
                       WRITE(*,*) "Warning: Mismatch soil layers for ", TRIM(var_in_file_name), ". File: ", file_num_layers, ", Model: ", nl_soil
                    ENDIF
                 ELSE
                    WRITE(*,*) "Warning: Unsupported dimensions (",ndims_in_file,") for var ", TRIM(var_in_file_name), ". Skipping."
                    tracer_soil_concentration(i_tracer, :, :) = 0.0_r8
                 ENDIF
              END DO ! i_tracer loop
              CALL ncio_close(nc_file_id)
              IF (init_file_mapper_built) THEN
                  CALL init_file_mapper%forc_free_mem() 
              ENDIF
              IF (allocated(init_file_grid%glon)) DEALLOCATE(init_file_grid%glon) ! Manual deallocation if needed
              IF (allocated(init_file_grid%glat)) DEALLOCATE(init_file_grid%glat)
              ! Add other grid component deallocations if necessary
           ENDIF ! file open error
        ENDIF ! file exists
     ENDIF ! p_is_master

#ifdef USEMPI
    IF (DEF_Tracer_Number > 0 .AND. numpatch > 0 .AND. nl_soil > 0) THEN
      IF (allocated(tracer_soil_concentration)) THEN
         CALL mpi_bcast(tracer_soil_concentration, SIZE(tracer_soil_concentration), MPI_REAL8, p_address_master, p_comm_glb, p_err)
      ENDIF
    ENDIF
#endif
  ELSE
     IF (p_is_master) WRITE(*,*) "Tracer Initialize: DEF_Tracer_Init is false. Using default (zero) concentrations."
  ENDIF
END SUBROUTINE tracer_initialize_concentrations

END MODULE MOD_Tracer_Initialize
