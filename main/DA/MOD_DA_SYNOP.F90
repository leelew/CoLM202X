#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_SYNOP
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Data assimilation of 2m temperature and 2m humidity from SYNOP stations
!
! AUTHOR:
!   Lu Li, 08/2025
!-----------------------------------------------------------------------------
   USE MOD_DataType
   USE MOD_SpatialMapping
   USE MOD_DA_ObsOperator
   USE MOD_DA_EnKF
   USE MOD_DA_Vars_TimeVariables
   USE MOD_Vars_Global, only: pi, nl_soil
   USE MOD_LandPatch
   USE MOD_Block
   USE MOD_Namelist
   USE MOD_Pixel
   USE MOD_Pixelset
   USE MOD_Mesh
   IMPLICIT NONE
   SAVE

! public functions
   PUBLIC :: init_DA_SYNOP
   PUBLIC :: run_DA_SYNOP
   PUBLIC :: end_DA_SYNOP

   PRIVATE

! local variables
   ! derived types of pixel index in each worker
   ! firstly sorted by pixel latitude index (nlat) and record 
   ! corresponding longitude index and patch id for each pixel 
   type :: idx_type
      integer, allocatable :: ilon (:)
      integer, allocatable :: ipatch (:)
   END type idx_type
   type(idx_type), allocatable :: idx (:)                    ! derived types of pixel index at each worker
   integer, allocatable :: counter (:)                       ! counter of pixel longitude index of each latitude index at each worker

   ! info of all SYNOP sites
   character(len=256) :: file_site                           ! file of location of all SYNOP sites
   real(r8), allocatable :: lat(:), lon(:)                   ! latitude and longitude of all SYNOP sites
   integer :: nsite                                          ! number of all SYNOP sites
   integer, allocatable :: iloc (:,:)                        ! global lat/lon index of pixel cover each site in all SYNOP sites
   integer :: counter_worker_nsite                           ! number of all SYNOP sites located at each worker
   integer, allocatable :: ip_worker (:,:)                   ! patch id of all SYNOP sites located at each worker
   integer, allocatable :: temp (:,:)                        ! temporary array for receiving data from workers
   integer, allocatable :: idx_lat(:), idx_lon(:), pos(:)    ! temporary array save global lat/lon index of each site
   integer, allocatable :: synop_lut (:,:)                   ! look-up-table of all SYNOP sites

   ! logical variables
   logical :: has_file                                       ! whether has file of SMAP at target day
   logical :: has_obs                                        ! whether has obs at current step
   logical :: has_DA                                         ! whether has data assimilation 

   ! time (UTC) at current step
   integer :: month, mday, hour                              ! month, day, hour of current step
   character(len=256) :: yearstr, monthstr, daystr, hourstr  ! string of year, month, day, hour
   integer :: idate_b, idate_e                               ! begin & end seconds since begin of current day (UTC)

   ! time variables used to determine whether has obs
   real(r8), allocatable :: synop_time(:)                    ! seconds of all obs since begin of current day (UTC)
   real(r8), allocatable :: dt_b(:)                          ! delta time between obs and begin seconds of current day
   real(r8), allocatable :: dt_e(:)                          ! delta time between obs and end seconds of current day

   ! observations (dimensions changes with time)
   integer :: num_obs                                        ! number of all obs in current file
   real(r8), allocatable :: synop_lat(:)                     ! latitude of all obs
   real(r8), allocatable :: synop_lon(:)                     ! longitude of all obs
   integer, allocatable :: synop_id(:)                       ! global id of all obs
   real(r8), allocatable :: synop_tref(:)                    ! 2m temperature of all obs ([K])
   integer,  allocatable :: synop_qref(:)                    ! 2m humidity of all obs ([K])

   ! info of SYNOP sites at each step
   character(len=256) :: file_synop                          ! SYNOP file path at each step
   integer, allocatable :: synop_idx (:,:)                   ! index of SYNOP sites (worker/patch id) at each step
   integer, allocatable :: site_id_worker (:)                ! site id of SYNOP sites at each step at each worker
   real(r8), allocatable :: tref_ens_worker (:,:)            ! predicted 2m temperature of SYNOP sites at each step at each worker
   real(r8), allocatable :: qref_ens_worker (:,:)            ! predicted 2m humidity of SYNOP sites at each step at each worker
   real(r8), allocatable :: qref_ens_o (:,:)                 ! predicted 2m temperature of SYNOP sites at each step
   real(r8), allocatable :: tref_ens_o (:,:)                 ! predicted 2m humidity of SYNOP sites at each step

   ! observations around patch (dimensions changes with patch)
   integer :: num_obs_p
   logical, allocatable :: index_p(:)                        ! index of obs around each patch
   real(r8), allocatable :: synop_lat_p(:)                   ! latitude of obs around each patch
   real(r8), allocatable :: synop_lon_p(:)                   ! longitude of obs around each patch
   real(r8), allocatable :: synop_qref_p(:)                  ! 2m temperature of obs around each patch ([K])
   real(r8), allocatable :: synop_tref_p(:)                  ! 2m humidity of obs around each patch ([K])
   real(r8), allocatable :: synop_p(:)                       ! concatenate 2m temperature and humidity of obs around each patch
   real(r8), allocatable :: qref_ens_p(:, :)                 ! predicted 2m temperature around patch
   real(r8), allocatable :: tref_ens_p(:,:)                  ! predicted 2m humidity around patch
   real(r8), allocatable :: pred_synop_ens_p(:,:)            ! predicted 2m temperature and humidity around patch
   real(r8), allocatable :: d_p(:)                           ! distance between obs and patch center

   ! parameters of LETKF
   real(r8), allocatable :: obs_err(:)                       ! observation error
   real(r8), parameter   :: dres = 0.4                       ! search localization radius (deg)
   real(r8), parameter   :: loc_r = 1.0                      ! localization radius
   real(r8), parameter   :: infl = 1.2                       ! inflation factor
   real(r8), parameter   :: static_obs_err_tref = 1.0        ! static observation error (2m temperature)
   real(r8), parameter   :: static_obs_err_qref = 0.1        ! static observation error (2m humidity)

   ! temporary variables for data assimilation
   real(r8), allocatable :: trans(:,:)                       ! transformation matrix on each patch
   real(r8), allocatable :: wice_soi_ens(:,:)                ! soil ice content
   real(r8), allocatable :: wice_soi_ens_da(:,:)             ! soil ice content after data assimilation
   real(r8), allocatable :: wliq_soi_ens(:,:)                ! soil liquid water content
   real(r8), allocatable :: wliq_soi_ens_da(:,:)             ! soil liquid water content after data assimilation
   logical, allocatable ::  filter(:)                        ! to mask the water
   real(r8) :: eff_porsl                                     ! effective porosity of soil layer
   


!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE init_DA_SYNOP()

!-----------------------------------------------------------------------------
      USE MOD_Spmd_Task
      USE MOD_Namelist
      USE MOD_Grid
      USE MOD_NetCDFSerial
      USE MOD_LandPatch
      USE MOD_Pixelset
      USE MOD_RangeCheck
      IMPLICIT NONE

!------------------------ Local Variables ------------------------------------
      integer :: np, ie, ipxstt, ipxend, ipxl, i, ilat, isite, ilon, iwork, mesg(2), isrc, ndata, numpxl_lat, numpxl_lon, numpxl

!-----------------------------------------------------------------------------

#ifndef SinglePoint
!#############################################################################
! Makeup derived types of pixel index for fast access at each worker
!#############################################################################
      IF (p_is_worker) THEN
         allocate (counter (pixel%nlat))
         counter(:) = 0

         ! count the number of pixel lon index for each pixel lat index
         DO np = 1, numpatch
            ie = landpatch%ielm(np)

            ipxstt = landpatch%ipxstt(np)
            ipxend = landpatch%ipxend(np)

            DO ipxl = ipxstt, ipxend
               counter(mesh(ie)%ilat(ipxl)) = counter(mesh(ie)%ilat(ipxl)) + 1
            ENDDO
         ENDDO

         ! allocate derived types of index 
         allocate (idx (pixel%nlat))
         DO i = 1, pixel%nlat
            IF (counter(i) > 0) THEN
               allocate (idx(i)%ilon(counter(i)))
               allocate (idx(i)%ipatch(counter(i)))

               idx(i)%ilon(:) = 0
               idx(i)%ipatch(:) = 0
            ENDIF
         ENDDO

         ! fill the index 
         counter(:) = 0
         DO np = 1, numpatch
            ie = landpatch%ielm(np)

            ipxstt = landpatch%ipxstt(np)
            ipxend = landpatch%ipxend(np)

            DO ipxl = ipxstt, ipxend
               ilat = mesh(ie)%ilat(ipxl)
               counter(mesh(ie)%ilat(ipxl)) = counter(mesh(ie)%ilat(ipxl)) + 1
               idx(ilat)%ilon(counter(mesh(ie)%ilat(ipxl))) = mesh(ie)%ilon(ipxl)
               idx(ilat)%ipatch(counter(mesh(ie)%ilat(ipxl))) = np
            ENDDO
         ENDDO
      ENDIF 


!#############################################################################
! Read the location of SYNOP sites and find the located pixel of each site
!#############################################################################
      ! file of SYNOP sites
      file_site = trim(DEF_DA_obsdir)//'/SYNOP.nc'

      ! read latitude and longitude of sites
      IF (ncio_var_exist(file_site, 'latitude')) THEN
         CALL ncio_read_bcast_serial(file_site, 'latitude', lat)
         CALL ncio_read_bcast_serial(file_site, 'longitude', lon)
      ENDIF
      nsite = size(lat)

      ! find the located pixel of each site & broadcast workers
      IF (p_is_master) THEN
         allocate (iloc (2, nsite))

         DO i = 1, nsite
            numpxl_lat = count(pixel%lat_s <= lat(i) .and. pixel%lat_n >= lat(i))
            numpxl_lon = count(pixel%lon_w <= lon(i) .and. pixel%lon_e >= lon(i))
            IF (allocated (idx_lat)) deallocate (idx_lat)
            IF (allocated (idx_lon)) deallocate (idx_lon)
            allocate (idx_lat(numpxl_lat))
            allocate (idx_lon(numpxl_lon))

            idx_lat = pack([(isite, isite=1, nsite)], pixel%lat_s <= lat(i) .and. pixel%lat_n >= lat(i))
            idx_lon = pack([(isite, isite=1, nsite)], pixel%lon_w <= lon(i) .and. pixel%lon_e >= lon(i))

            iloc(1,i) = idx_lat(1)
            iloc(2,i) = idx_lon(1)
         ENDDO

         CALL mpi_bcast(iloc, 2*nsite, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      ENDIF


!#############################################################################
! Assess the patch id of pixel that cover each observation at each worker
!#############################################################################
      IF (p_is_worker) THEN
         ! count the number of site that located in pixels of each worker
         counter_worker_nsite = 0
         DO i = 1, nsite
            IF (counter(iloc(1, i)) > 0) THEN
               IF (any(idx(iloc(1, i))%ilon == iloc(2, i))) THEN
                  counter_worker_nsite = counter_worker_nsite + 1
               ENDIF
            ENDIF
         ENDDO

         ! assess patch/site id of sites located at each worker
         IF (counter_worker_nsite > 0) THEN
            allocate (ip_worker (2, counter_worker_nsite))
            counter_worker_nsite = 0

            DO i = 1, nsite
               IF (counter(iloc(1, i)) > 0) THEN
                  IF (any(idx(iloc(1, i))%ilon == iloc(2, i))) THEN
                     numpxl = count(idx(iloc(1, i))%ilon == iloc(2, i))
                     allocate (pos(numpxl))

                     pos = pack([(ilon, ilon=1,counter(iloc(1, i)))], (idx(iloc(1, i))%ilon == iloc(2, i))) 
                     counter_worker_nsite = counter_worker_nsite + 1
                     ip_worker(1, counter_worker_nsite) = i
                     ip_worker(2, counter_worker_nsite) = idx(iloc(1, i))%ipatch(pos(1))
                  ENDIF
               ENDIF
            ENDDO
         ENDIF

         ! send the number of site and their patch id to master
#ifdef USEMPI
         IF (counter_worker_nsite > 0) THEN
            mesg = (/p_iam_glb, counter_worker_nsite/)
            CALL mpi_send(mesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
            CALL mpi_send(ip_worker, 2*counter_worker_nsite, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
         ENDIF
#endif

         ! deallocate 
         deallocate (counter)
         DO i = 1, pixel%nlat
            deallocate (idx(i)%ilon)
            deallocate (idx(i)%ipatch)
         ENDDO
         deallocate (idx)   
         deallocate (iloc)
         deallocate (ip_worker)
      ENDIF


!#############################################################################
! Generate look-up-table (contains the worker/patch id of each site) at master
!#############################################################################
      IF (p_is_master) THEN
         allocate (synop_lut (2, nsite))

#ifdef USEMPI
         DO iwork = 1, p_np_worker
            CALL mpi_recv(mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
            isrc = mesg(1)
            ndata = mesg(2)

            IF (ndata > 0) THEN
               allocate(temp(2, ndata))

               CALL mpi_recv(temp, 2*ndata, MPI_INTEGER, isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               DO i = 1, ndata
                  synop_lut(1, temp(1,i)) = isrc
                  synop_lut(2, temp(1,i)) = temp(2,i)
               ENDDO

               deallocate(temp)
            ENDIF
         ENDDO
#endif
      ENDIF

#endif

   END SUBROUTINE init_DA_SYNOP

!-----------------------------------------------------------------------------

   SUBROUTINE run_DA_SYNOP(idate, deltim)

!-----------------------------------------------------------------------------
      USE MOD_Spmd_task
      USE MOD_TimeManager
      USE MOD_NetCDFBlock
      USE MOD_Mesh
      USE MOD_LandElm
      USE MOD_LandPatch
      USE MOD_Vars_Global
      USE MOD_Vars_1DFluxes
      USE MOD_Vars_1DForcing
      USE MOD_Vars_TimeVariables
      USE MOD_Vars_TimeInvariants
      USE MOD_DA_Vars_TimeVariables
      USE MOD_RangeCheck
      USE MOD_UserDefFun
      USE MOD_DA_EnKF
      USE MOD_Const_Physical, only: denice, denh2o
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      integer, intent(in) :: idate(3)
      real(r8), intent(in) :: deltim

!------------------------ Local Variables ------------------------------------
      real(r8) :: lat_p_n, lat_p_s, lon_p_w, lon_p_e
      integer  :: ib, jb, il, jl, ip, iens, np, i, n
      integer  :: sdate(3)
      integer  :: isrc, ndata, iwork, mesg(2)

!-----------------------------------------------------------------------------

!#############################################################################
! Identify if there are observations at this time step 
!#############################################################################
      ! Do not perform DA, only calcuate predict BRT for diagnostic
      IF (DEF_DA_ENS == 1) THEN
         has_file = .false.
         has_obs = .false.
      ELSE
         ! covert local time to UTC for single point
         sdate = idate
         CALL adj2begin(sdate)
#ifdef SinglePoint
         IF (.not. DEF_simulation_time%greenwich) THEN
            CALL localtime2gmt(sdate)
         ENDIF
#endif
         ! calculate year/month/day/hour of current step
         CALL julian2monthday(sdate(1), sdate(2), month, mday)
         hour = int(sdate(3)/3600)

         ! whether has file of synop at target day
         write (yearstr, '(I4.4)') sdate(1)
         write (monthstr, '(I2.2)') month
         write (daystr, '(I2.2)') mday
         write (hourstr, '(I2.2)') hour
         file_synop = trim(DEF_DA_obsdir)//'/pre/'//'/SYNOP_'// &
               trim(yearstr)//'_'//trim(monthstr)//'_'//trim(daystr)//'_'//trim(hourstr)//'.nc'
         inquire (file=trim(file_synop), exist=has_file)

         ! whether have obs at this time interval
         has_obs = .false.
         IF (has_file) THEN
            CALL ncio_read_bcast_serial(file_synop, 'time', synop_time)
            num_obs = size(synop_time)
            idate_b = sdate(3)
            idate_e = sdate(3) + deltim
            allocate (dt_b(num_obs))
            allocate (dt_e(num_obs))
            dt_b = synop_time - idate_b
            dt_e = synop_time - idate_e
            IF (any(dt_b >= 0 .and. dt_e <= 0)) has_obs = .true.
            deallocate (dt_b)
            deallocate (dt_e)
         ELSE
            has_obs = .false.
         ENDIF
      ENDIF


!#############################################################################
! Calculate predicted observations & mapping to world grid & read observations
!#############################################################################
      ! read observations from file
      IF (has_obs) THEN
         CALL ncio_read_bcast_serial(file_synop, 'lat',  synop_lat )
         CALL ncio_read_bcast_serial(file_synop, 'lon',  synop_lon )
         CALL ncio_read_bcast_serial(file_synop, 'id',   synop_id  )
         CALL ncio_read_bcast_serial(file_synop, 'tref',  synop_tref )
         CALL ncio_read_bcast_serial(file_synop, 'qref', synop_qref)
      ENDIF 

      ! allocate memory
      IF (has_obs) THEN
         allocate (synop_idx (2, num_obs))
         allocate (qref_ens_o (num_obs, DEF_DA_ENS))
         allocate (tref_ens_o (num_obs, DEF_DA_ENS))
      ENDIF
      IF (p_is_worker) THEN
         IF (counter_worker_nsite > 0) THEN
            allocate (tref_ens_worker (counter_worker_nsite, DEF_DA_ENS))
            allocate (qref_ens_worker (counter_worker_nsite, DEF_DA_ENS))
            allocate (site_id_worker  (counter_worker_nsite))
         ENDIF
      ENDIF

      ! crop corresponding index (worker id and patch id) of each observation
      IF (p_is_master) THEN
         synop_idx = synop_lut(:, synop_id)
         CALL mpi_bcast(synop_idx, 2*num_obs, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      ENDIF

      ! crop observations at each worker according index & send to master
      IF (p_is_worker) THEN
         counter_worker_nsite = 0
         DO i = 1, num_obs
            IF (synop_idx(1, i) == p_iam_glb) THEN
               counter_worker_nsite = counter_worker_nsite + 1
               tref_ens_worker(counter_worker_nsite,:) = tref_ens(:, synop_idx(2, i))
               qref_ens_worker(counter_worker_nsite,:) = qref_ens(:, synop_idx(2, i))
               site_id_worker (counter_worker_nsite)    = i
            ENDIF
         ENDDO
         
#ifdef USEMPI
         ! send the number of site and their patch id to master
         IF (counter_worker_nsite > 0) THEN
            mesg = (/p_iam_glb, counter_worker_nsite/)
            CALL mpi_send(mesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
            CALL mpi_send(tref_ens_worker, DEF_DA_ENS*counter_worker_nsite, MPI_REAL8, p_address_master, mpi_tag_data, p_comm_glb, p_err)
            CALL mpi_send(qref_ens_worker, DEF_DA_ENS*counter_worker_nsite, MPI_REAL8, p_address_master, mpi_tag_data, p_comm_glb, p_err)
            CALL mpi_send(site_id_worker, counter_worker_nsite, MPI_INTEGER, p_address_master, mpi_tag_data, p_comm_glb, p_err)
         ENDIF
#endif

         deallocate (tref_ens_worker)
         deallocate (qref_ens_worker)
         deallocate (site_id_worker)
      ENDIF

      ! concatenate all predicted observations at master & broadcast to all workers
      IF (p_is_master) THEN
#ifdef USEMPI
         DO iwork = 1, p_np_worker
            CALL mpi_recv(mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
            isrc = mesg(1)
            ndata = mesg(2)

            IF (ndata > 0) THEN
               allocate(tref_ens_worker(ndata, DEF_DA_ENS))
               allocate(qref_ens_worker(ndata, DEF_DA_ENS))
               allocate(site_id_worker (ndata))

               CALL mpi_recv(tref_ens_worker, ndata*DEF_DA_ENS, MPI_REAL8, isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL mpi_recv(qref_ens_worker, ndata*DEF_DA_ENS, MPI_REAL8, isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL mpi_recv(site_id_worker, ndata, MPI_INTEGER, isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               qref_ens_o(site_id_worker,:) = qref_ens_worker
               tref_ens_o(site_id_worker,:) = tref_ens_worker
            ENDIF

            deallocate(tref_ens_worker)
            deallocate(qref_ens_worker)
            deallocate(site_id_worker)
         ENDDO
#endif
         CALL mpi_bcast(qref_ens_o, DEF_DA_ENS*num_obs, MPI_REAL8, p_address_master, p_comm_glb, p_err)
         CALL mpi_bcast(tref_ens_o, DEF_DA_ENS*num_obs, MPI_REAL8, p_address_master, p_comm_glb, p_err)
      ENDIF


!#############################################################################
! Run data assimilation for SYNOP observations
!#############################################################################
      IF (p_is_worker) THEN
         DO np = 1, numpatch
            ! regions info around target patch
            lat_p_n = patchlatr(np)*180/pi + dres
            lat_p_s = patchlatr(np)*180/pi - dres
            lon_p_w = patchlonr(np)*180/pi - dres
            lon_p_e = patchlonr(np)*180/pi + dres

            ! find observations around each patch
            num_obs_p = count( &
               synop_lat(:) < lat_p_n .and. synop_lat(:) > lat_p_s .and. &
               synop_lon(:) > lon_p_w .and. synop_lon(:) < lon_p_e .and. &
               synop_time(:) - idate_b >= 0 .and. synop_time(:) - idate_e <= 0 .and. &
               tref_ens_o(:,1) > 0)

            ! perform data assimilation if have observations
            IF (num_obs_p > 0) THEN
               allocate (index_p         (num_obs_p              ))
               allocate (synop_lat_p     (num_obs_p              ))
               allocate (synop_lon_p     (num_obs_p              ))
               allocate (synop_qref_p    (num_obs_p              ))
               allocate (synop_tref_p    (num_obs_p              ))
               allocate (synop_p         (2*num_obs_p            ))
               allocate (qref_ens_p      (num_obs_p, DEF_DA_ENS  ))
               allocate (tref_ens_p      (num_obs_p, DEF_DA_ENS  ))
               allocate (pred_synop_ens_p(2*num_obs_p, DEF_DA_ENS))
               allocate (d_p             (2*num_obs_p            ))
               allocate (obs_err         (2*num_obs_p            ))
               allocate (trans           (DEF_DA_ENS,DEF_DA_ENS  ))
               allocate (wice_soi_ens    (nl_soil,   DEF_DA_ENS  ))
               allocate (wice_soi_ens_da (nl_soil,   DEF_DA_ENS  ))
               allocate (wliq_soi_ens    (nl_soil,   DEF_DA_ENS  ))
               allocate (wliq_soi_ens_da (nl_soil,   DEF_DA_ENS  ))

               ! index of observations around target patch
               index_p = (synop_lat(:) < lat_p_n .and. synop_lat(:) > lat_p_s .and. &
                  synop_lon(:) > lon_p_w .and. synop_lon(:) < lon_p_e .and. &
                  synop_time(:) - idate_b >= 0 .and. synop_time(:) - idate_e <= 0 .and. &
                  tref_ens_o(:,1) > 0)

               ! crop observations around target patch
               synop_lat_p  = pack(synop_lat , index_p)
               synop_lon_p  = pack(synop_lon , index_p)
               synop_qref_p = pack(synop_qref, index_p)
               synop_tref_p = pack(synop_tref, index_p)

               ! predicted observations around target patch
               DO i = 1, DEF_DA_ENS
                  qref_ens_p(:, i) = pack(qref_ens_o(:, i), index_p)
                  tref_ens_p(:, i) = pack(tref_ens_o(:, i), index_p)
               ENDDO

               ! concatenate predicted observations & observations
               synop_p(1:num_obs_p) = synop_tref_p
               synop_p(num_obs_p+1:2*num_obs_p) = synop_qref_p
               pred_synop_ens_p(1:num_obs_p, 1:DEF_DA_ENS) = tref_ens_p
               pred_synop_ens_p(num_obs_p+1:2*num_obs_p, 1:DEF_DA_ENS) = qref_ens_p

               ! calculate distance between observations and target patch
               d_p(1:num_obs_p) = 2*6.3781e3*asin(sqrt(sin((synop_lat_p*pi/180 - patchlatr(np))/2.0)**2 + &
                  cos(synop_lat_p*pi/180)*cos(patchlatr(np))*sin((synop_lon_p*pi/180 - patchlonr(np))/2.0)**2))
               d_p(num_obs_p+1:2*num_obs_p) = 2*6.3781e3*asin(sqrt(sin((synop_lat_p*pi/180 - patchlatr(np))/2.0)**2 + &
                  cos(synop_lat_p*pi/180)*cos(patchlatr(np))*sin((synop_lon_p*pi/180 - patchlonr(np))/2.0)**2))

               ! setting observation error matrix
               obs_err(1:num_obs_p) = static_obs_err_tref
               obs_err(num_obs_p+1:2*num_obs_p) = static_obs_err_qref

               ! calculate transformation matrix
               CALL letkf(DEF_DA_ENS, 2*num_obs_p, &
                  pred_synop_ens_p, synop_p, obs_err, d_p, loc_r, infl, &
                  trans)

               ! calculate analysis value
               IF (wliq_soisno_ens(1, 1, np) /= spval) THEN
                  has_DA = .TRUE.

                  ! soil layer
                  DO iens = 1, DEF_DA_ENS
                     wliq_soi_ens(:, iens) = wliq_soisno_ens(1:, iens, np)
                     wice_soi_ens(:, iens) = wice_soisno_ens(1:, iens, np)
                  ENDDO

                  ! analysis
                  CALL dgemm('N', 'N', nl_soil, DEF_DA_ENS, DEF_DA_ENS, 1.0_8, wliq_soi_ens, &
                     nl_soil, trans, DEF_DA_ENS, 0.0_8, wliq_soi_ens_da, nl_soil)
                  CALL dgemm('N', 'N', nl_soil, DEF_DA_ENS, DEF_DA_ENS, 1.0_8, wice_soi_ens, &
                     nl_soil, trans, DEF_DA_ENS, 0.0_8, wice_soi_ens_da, nl_soil)

                  DO iens = 1, DEF_DA_ENS
                     wliq_soisno_ens(1:, iens, np) = wliq_soi_ens_da(1:, iens)
                     wice_soisno_ens(1:, iens, np) = wice_soi_ens_da(1:, iens)
                  ENDDO

                  ! limit the soil liquid and ice water in a reasonable range
                  DO i = 1, nl_soil
                     DO iens = 1, DEF_DA_ENS
                        ! lower bound
                        wliq_soisno_ens(i, iens, np) = max(0.0d0, wliq_soisno_ens(i, iens, np))
                        wice_soisno_ens(i, iens, np) = max(0.0d0, wice_soisno_ens(i, iens, np))
                        IF (wliq_soisno_ens(i, iens, np) == 0.0 .and. wice_soisno_ens(i, iens, np) == 0.0) THEN
                           IF (t_soisno_ens(i, iens, np) < -5.0) THEN
                              wice_soisno_ens(i, iens, np) = 1e-10
                           ELSE
                              wliq_soisno_ens(i, iens, np) = 1e-10
                           ENDIF
                        ENDIF

                        ! upper bound
                        wice_soisno_ens(i, iens, np) = min(porsl(i, np)*(dz_soi(i)*denice), wice_soisno_ens(i, iens, np))
                        eff_porsl = max(0.0d0, porsl(i, np) - wice_soisno_ens(i, iens, np)/(dz_soi(i)*denice))
                        wliq_soisno_ens(i, iens, np) = min(eff_porsl*(dz_soi(i)*denh2o), wliq_soisno_ens(i, iens, np))
                     ENDDO
                  ENDDO

                  ! move residual water to water table
                  wa_ens(:, np) = wa_ens(:, np) - sum(wliq_soisno_ens(1:, :, np) + wice_soisno_ens(1:, :, np) - wliq_soi_ens - wice_soi_ens, dim=1)

                  ! update volumetric water content for diagnostic
                  DO iens = 1, DEF_DA_ENS
                     h2osoi_ens(:, iens, np) = wliq_soisno_ens(1:, iens, np)/(dz_soi(:)*denh2o) + wice_soisno_ens(1:, iens, np)/(dz_soi(:)*denice)
                     h2osoi_ens(:, iens, np) = min(1.0d0, h2osoi_ens(:, iens, np))
                     h2osoi_ens(:, iens, np) = max(0.0d0, h2osoi_ens(:, iens, np))
                  ENDDO
               ENDIF
            ENDIF

            ! deallocate memory (cuz dimensions changes with patch)
            IF (allocated(index_p          )) deallocate (index_p           )
            IF (allocated(synop_lat_p      )) deallocate (synop_lat_p       )
            IF (allocated(synop_lon_p      )) deallocate (synop_lon_p       )
            IF (allocated(synop_qref_p     )) deallocate (synop_qref_p      )
            IF (allocated(synop_tref_p     )) deallocate (synop_tref_p      )
            IF (allocated(synop_p          )) deallocate (synop_p           )
            IF (allocated(qref_ens_p       )) deallocate (qref_ens_p        )
            IF (allocated(tref_ens_p       )) deallocate (tref_ens_p        )
            IF (allocated(pred_synop_ens_p )) deallocate (pred_synop_ens_p  )
            IF (allocated(d_p              )) deallocate (d_p               )
            IF (allocated(obs_err          )) deallocate (obs_err           )
            IF (allocated(trans            )) deallocate (trans             )
            IF (allocated(wice_soi_ens     )) deallocate (wice_soi_ens      )
            IF (allocated(wice_soi_ens_da  )) deallocate (wice_soi_ens_da   )
            IF (allocated(wliq_soi_ens     )) deallocate (wliq_soi_ens      )
            IF (allocated(wliq_soi_ens_da  )) deallocate (wliq_soi_ens_da   )
         ENDDO
      ENDIF

   END SUBROUTINE run_DA_SYNOP

!-----------------------------------------------------------------------------

   SUBROUTINE end_DA_SYNOP()

!-----------------------------------------------------------------------------
      IMPLICIT NONE

!-----------------------------------------------------------------------------
      IF (allocated(synop_time             )) deallocate (synop_time             )
      IF (allocated(dt_b                   )) deallocate (dt_b                   )
      IF (allocated(dt_e                   )) deallocate (dt_e                   )
      IF (allocated(synop_lat              )) deallocate (synop_lat              )
      IF (allocated(synop_lon              )) deallocate (synop_lon              )
      IF (allocated(synop_id               )) deallocate (synop_id               )
      IF (allocated(synop_tref             )) deallocate (synop_tref             )
      IF (allocated(synop_qref             )) deallocate (synop_qref             )
      IF (allocated(qref_ens_o             )) deallocate (qref_ens_o             )
      IF (allocated(tref_ens_o             )) deallocate (tref_ens_o             )
      IF (allocated(obs_err                )) deallocate (obs_err                )
      IF (allocated(synop_idx              )) deallocate (synop_idx              )
      IF (allocated(index_p                )) deallocate (index_p                )
      IF (allocated(synop_lat_p            )) deallocate (synop_lat_p            )
      IF (allocated(synop_lon_p            )) deallocate (synop_lon_p            )
      IF (allocated(synop_qref_p           )) deallocate (synop_qref_p           )
      IF (allocated(synop_tref_p           )) deallocate (synop_tref_p           )
      IF (allocated(synop_p                )) deallocate (synop_p                )
      IF (allocated(qref_ens_p             )) deallocate (qref_ens_p             )
      IF (allocated(tref_ens_p             )) deallocate (tref_ens_p             )
      IF (allocated(pred_synop_ens_p       )) deallocate (pred_synop_ens_p       )
      IF (allocated(d_p                    )) deallocate (d_p                    )
      IF (allocated(trans                  )) deallocate (trans                  )
      IF (allocated(wice_soi_ens           )) deallocate (wice_soi_ens           )
      IF (allocated(wice_soi_ens_da        )) deallocate (wice_soi_ens_da        )
      IF (allocated(wliq_soi_ens           )) deallocate (wliq_soi_ens           )
      IF (allocated(wliq_soi_ens_da        )) deallocate (wliq_soi_ens_da        )
   END SUBROUTINE end_DA_SYNOP

!-----------------------------------------------------------------------------
END MODULE MOD_DA_SYNOP
#endif