#include <define.h>

PROGRAM MKSRFDATA
   ! ======================================================================
   ! Surface grid edges:
   ! The model domain was defined with the north, east, south, west edges:
   !          edgen: northern edge of grid : > -90 and <= 90 (degrees)
   !          edgee: eastern edge of grid  : > western edge and <= 180
   !          edges: southern edge of grid : >= -90  and <  90
   !          edgew: western edge of grid  : >= -180 and < 180
   !
   ! Region (global) latitude grid goes from:
   !                 NORTHERN edge (POLE) to SOUTHERN edge (POLE)
   ! Region (global) longitude grid starts at:
   !                 WESTERN edge (DATELINE with western edge)
   !                 West of Greenwich defined negative for global grids,
   !                 the western edge of the longitude grid starts at the dateline
   !
   ! Land characteristics at the 30 arc-seconds grid resolution (RAW DATA):
   !              1. Global Terrain Dataset (elevation height,...)
   !              2. Global Land Cover Characteristics (land cover TYPE, plant leaf area index, Forest Height, ...)
   !              3. Global Lakes and Wetlands Characteristics (lake and wetlands types, lake coverage and lake depth)
   !              4. Global Glacier Characteristics
   !              5. Global Urban Characteristics (urban extent, ...)
   !              6. Global Soil Characteristics (...)
   !              7. Global Cultural Characteristics (ON-GONG PROJECT)
   !
   ! Land charateristics at the model grid resolution (CREATED):
   !              1. Model grid (longitude, latitude)
   !              2. Fraction (area) of patches of grid (0-1)
   !                 2.1 Fraction of land water bodies (lake, reservoir, river)
   !                 2.2 Fraction of wetland
   !                 2.3 Fraction of glacier
   !                 2.4 Fraction of urban and built-up
   !                 ......
   !              3. Plant leaf area index
   !              4. Tree height
   !              5. Lake depth
   !              6. Soil thermal and hydraulic parameters
   !
   ! Created by Yongjiu Dai, 02/2014
   !
   !
   ! ======================================================================
   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Block
   USE MOD_Pixel
   USE MOD_Grid
   USE MOD_Mesh
   USE MOD_MeshFilter
   USE MOD_LandElm
#ifdef CATCHMENT
   USE MOD_LandHRU
#endif
   USE MOD_LandPatch
   USE MOD_SrfdataRestart
   USE MOD_Const_LC
#ifdef PFT_CLASSIFICATION
   USE MOD_LandPFT
#endif
#ifdef PC_CLASSIFICATION
   USE MOD_LandPC
#endif
#ifdef URBAN_MODEL
   USE MOD_LandUrban
#endif

#ifdef SrfdataDiag
   USE MOD_SrfdataDiag, only : gdiag, srfdata_diag_init
#endif

   USE MOD_RegionClip


   IMPLICIT NONE

   CHARACTER(len=256) :: nlfile

   CHARACTER(LEN=256) :: dir_rawdata
   CHARACTER(LEN=256) :: dir_landdata
   REAL(r8) :: edgen  ! northern edge of grid (degrees)
   REAL(r8) :: edgee  ! eastern edge of grid (degrees)
   REAL(r8) :: edges  ! southern edge of grid (degrees)
   REAL(r8) :: edgew  ! western edge of grid (degrees)

   TYPE (grid_type) :: gridlai, gnitrif, gndep, gfire, gtopo
   TYPE (grid_type) :: grid_urban_5km, grid_urban_100m, grid_urban_500m

   INTEGER*8 :: start_time, end_time, c_per_sec, time_used


#ifdef USEMPI
   CALL spmd_init ()
#endif

   IF (p_is_master) THEN
      CALL system_clock (start_time)
   ENDIF

   CALL getarg(1, nlfile)

   CALL read_namelist (nlfile)

#ifdef SinglePoint
   CALL read_surface_data_single (SITE_fsrfdata, mksrfdata=.true.)
#endif

   IF (USE_srfdata_from_larger_region) THEN

      CALL srfdata_region_clip (DEF_dir_existing_srfdata, DEF_dir_landdata)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
      CALL spmd_exit
#endif
      CALL exit()
   ENDIF

   IF (USE_srfdata_from_3D_gridded_data) THEN

      ! TODO
      ! CALL srfdata_retrieve_from_3D_data (DEF_dir_existing_srfdata, DEF_dir_landdata)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
      CALL spmd_exit
#endif
      CALL exit()
   ENDIF

   dir_rawdata  = DEF_dir_rawdata
   dir_landdata = DEF_dir_landdata
   edges = DEF_domain%edges
   edgen = DEF_domain%edgen
   edgew = DEF_domain%edgew
   edgee = DEF_domain%edgee

   ! define blocks
   CALL gblock%set_by_size (DEF_nx_blocks, DEF_ny_blocks)
   ! CALL gblock%set_by_file (DEF_file_block)

   CALL Init_GlovalVars
   CAll Init_LC_Const

   ! ...........................................................................
   ! 1. Read in or create the modeling grids coordinates and related information
   ! ...........................................................................

   ! define domain in pixel coordinate
   CALL pixel%set_edges (edges, edgen, edgew, edgee)
   CALL pixel%assimilate_gblock ()

   ! define grid coordinates of mesh
#ifdef GRIDBASED
   CALL init_gridbased_mesh_grid ()
#endif

#ifdef CATCHMENT
   CALL gridmesh%define_by_name ('merit_90m')
#endif

#ifdef UNSTRUCTURED
   CALL gridmesh%define_from_file (DEF_file_mesh)
#endif

   ! define grid coordinates of mesh filter
   has_mesh_filter = inquire_mesh_filter ()
   IF (has_mesh_filter) THEN
      CALL grid_filter%define_from_file (DEF_file_mesh_filter)
   ENDIF

   ! define grid coordinates of hydro units in catchment
#ifdef CATCHMENT
   CALL ghru%define_by_name ('merit_90m')
#endif

   ! define grid coordinates of land types
#ifdef USGS_CLASSIFICATION
   CALL gpatch%define_by_name ('colm_1km')
#endif
#ifdef IGBP_CLASSIFICATION
   CALL gpatch%define_by_name ('colm_500m')
#endif
#ifdef PFT_CLASSIFICATION
   CALL gpatch%define_by_name ('colm_500m')
#endif
#ifdef PC_CLASSIFICATION
   CALL gpatch%define_by_name ('colm_500m')
#endif
#ifdef BGC
#if (defined CROP)
   ! define grid for crop parameters
   CALL gcrop%define_by_ndims (720,360)
#endif
#if (defined Fire)
   ! define grid for crop parameters
   CALL gfire%define_by_ndims (720,360)
#endif
#ifdef NITRIF
   CALL gnitrif%define_by_name ('nitrif_2deg')
#endif
   CALL gndep%define_by_name ('nitrif_2deg')
#endif

   ! define grid for land characteristics
   CALL gridlai%define_by_name ('colm_500m')

   ! define grid for topography
   CALL gtopo%define_by_name ('colm_500m')

   ! add by dong, only test for making urban data
#ifdef URBAN_MODEL
   CALL gurban%define_by_name          ('colm_500m')
   CALL grid_urban_500m%define_by_name ('colm_500m')
   CALL grid_urban_5km%define_by_name  ('colm_5km' )
   CALL grid_urban_100m%define_by_name ('colm_100m')

   CALL pixel%assimilate_grid (gurban         )
   CALL pixel%assimilate_grid (grid_urban_500m)
   CALL pixel%assimilate_grid (grid_urban_5km )
   CALL pixel%assimilate_grid (grid_urban_100m)

   CALL pixel%map_to_grid (gurban         )
   CALL pixel%map_to_grid (grid_urban_500m)
   CALL pixel%map_to_grid (grid_urban_5km )
   CALL pixel%map_to_grid (grid_urban_100m)
#endif

   ! assimilate grids to build pixels
#ifndef SinglePoint
   CALL pixel%assimilate_grid (gridmesh)
#endif
   IF (has_mesh_filter) THEN
      CALL pixel%assimilate_grid (grid_filter)
   ENDIF
#ifdef CATCHMENT
   CALL pixel%assimilate_grid (ghru)
#endif
   CALL pixel%assimilate_grid (gpatch)
   CALL pixel%assimilate_grid (gridlai)
#ifdef BGC
#if (defined CROP)
   CALL pixel%assimilate_grid (gcrop )
#endif
#if (defined Fire)
   CALL pixel%assimilate_grid (gfire )
#endif
#ifdef NITRIF
   CALL pixel%assimilate_grid (gnitrif)
#endif
   CALL pixel%assimilate_grid (gndep)
#endif

   CALL pixel%assimilate_grid (gtopo)

   ! map pixels to grid coordinates
#ifndef SinglePoint
   CALL pixel%map_to_grid (gridmesh)
#endif
   IF (has_mesh_filter) THEN
      CALL pixel%map_to_grid (grid_filter)
   ENDIF
#ifdef CATCHMENT
   CALL pixel%map_to_grid (ghru)
#endif
   CALL pixel%map_to_grid (gpatch)
#if (defined CROP)
   CALL pixel%map_to_grid (gcrop )
#endif
#if (defined Fire)
   CALL pixel%map_to_grid (gfire )
#endif
   CALL pixel%map_to_grid (gridlai)
#ifdef NITRIF
   CALL pixel%map_to_grid (gnitrif)
#endif
#ifdef BGC
   CALL pixel%map_to_grid (gndep)
#endif

   CALL pixel%map_to_grid (gtopo)

   ! build land elms
   CALL mesh_build ()
   CALL landelm_build

#ifdef GRIDBASED
   IF (.not. read_mesh_from_file) THEN
      CALL mesh_filter (gpatch, trim(DEF_dir_rawdata)//'/landtype_update.nc', 'landtype')
   ENDIF
#endif

   ! Filtering pixels
   IF (has_mesh_filter) THEN
      CALL mesh_filter (grid_filter, DEF_file_mesh_filter, 'mesh_filter')
   ENDIF

#ifdef CATCHMENT
   CALL landhru_build
#endif

   ! build land patches
   CALL landpatch_build

#ifdef PFT_CLASSIFICATION
   CALL landpft_build
#endif

#ifdef PC_CLASSIFICATION
   CALL landpc_build
#endif

#ifdef URBAN_MODEL
   CALL landurban_build
#endif

   ! ................................................................
   ! 2. SAVE land surface tessellation information
   ! ................................................................

   CALL gblock%save_to_file    (dir_landdata)

   CALL pixel%save_to_file     (dir_landdata)

   CALL mesh_save_to_file      (dir_landdata)

   CALL pixelset_save_to_file  (dir_landdata, 'landelm', landelm)

#ifdef CATCHMENT
   CALL pixelset_save_to_file  (dir_landdata, 'landhru', landhru)
#endif

   !print*, count(landpatch%settyp==13)
   CALL pixelset_save_to_file  (dir_landdata, 'landpatch', landpatch)

#ifdef PFT_CLASSIFICATION
   CALL pixelset_save_to_file  (dir_landdata, 'landpft'  , landpft  )
#endif

#ifdef PC_CLASSIFICATION
   CALL pixelset_save_to_file  (dir_landdata, 'landpc'   , landpc   )
#endif

#ifdef URBAN_MODEL
   CALL pixelset_save_to_file  (dir_landdata, 'landurban', landurban)
#endif

   ! ................................................................
   ! 3. Mapping land characteristic parameters to the model grids
   ! ................................................................
#ifdef SrfdataDiag
#ifdef GRIDBASED
   CALL gdiag%define_by_copy (gridmesh)
#else
   CALL gdiag%define_by_ndims(720,360)
#endif

   CALL srfdata_diag_init (dir_landdata)
#endif

#ifdef BGC
   call Aggregation_NDeposition            (gndep, dir_rawdata, dir_landdata)
#if (defined CROP)
   call Aggregation_CropParameters (gcrop, dir_rawdata, dir_landdata)
#endif
#ifdef Fire
   call Aggregation_Fire            (gfire, dir_rawdata, dir_landdata)
#endif
#if (defined NITRIF)
  call Aggregation_NitrifParameters (gnitrif, dir_rawdata, dir_landdata)
#endif
#endif

   CALL Aggregation_PercentagesPFT     (gpatch,  dir_rawdata, dir_landdata)

   CALL Aggregation_LakeDepth       (gpatch,  dir_rawdata, dir_landdata)

   CALL Aggregation_SoilParameters (gpatch,  dir_rawdata, dir_landdata)

   CALL Aggregation_SoilBrightness (gpatch,  dir_rawdata, dir_landdata)

#ifdef USE_DEPTH_TO_BEDROCK
   CALL Aggregation_DBedrock        (gpatch,  dir_rawdata, dir_landdata)
#endif

   CALL Aggregation_LAI             (gridlai, dir_rawdata, dir_landdata)

   CALL Aggregation_ForestHeight   (gpatch,  dir_rawdata, dir_landdata)

   CALL Aggregation_Topography      (gtopo,   dir_rawdata, dir_landdata)

#ifdef URBAN_MODEL
   CALL Aggregation_Urban (dir_rawdata, dir_landdata, DEF_LC_YEAR, &
                           grid_urban_5km, grid_urban_100m, grid_urban_500m)
#endif


   ! ................................................................
   ! 4. Free memories.
   ! ................................................................

#ifdef SinglePoint
#if (defined PFT_CLASSIFICATION)
   CALL write_surface_data_single (numpatch, numpft)
#else
   CALL write_surface_data_single (numpatch)
#endif
   CALL single_srfdata_final ()
#endif

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

   IF (p_is_master) THEN
      CALL system_clock (end_time, count_rate = c_per_sec)
      time_used = (end_time - start_time) / c_per_sec
      IF (time_used >= 3600) THEN
         write(*,101) time_used/3600, mod(time_used,3600)/60, mod(time_used,60)
         101 format (/, 'Overall system time used:', I4, ' hours', I3, ' minutes', I3, ' seconds.')
      ELSEIF (time_used >= 60) THEN
         write(*,102) time_used/60, mod(time_used,60)
         102 format (/, 'Overall system time used:', I3, ' minutes', I3, ' seconds.')
      ELSE
         write(*,103) time_used
         103 format (/, 'Overall system time used:', I3, ' seconds.')
      ENDIF

      write(*,*)  'Successful in surface data making.'
   ENDIF

#ifdef USEMPI
   CALL spmd_exit
#endif

END PROGRAM MKSRFDATA
! ----------------------------------------------------------------------
! EOP