!  SVN:$Id: ice_transport_driver.F90 925 2015-03-04 00:34:27Z eclare $
!=======================================================================
!
! Drivers for remapping and upwind ice transport
!
! authors: Elizabeth C. Hunke and William H. Lipscomb, LANL 
!
! 2004: Revised by William Lipscomb from ice_transport_mpdata.
!       Stripped out mpdata, retained upwind, and added block structure.
! 2006: Incorporated remap transport driver and renamed from
!       ice_transport_upwind.  
! 2011: ECH moved edgearea arrays into ice_transport_remap.F90

      module ice_transport_driver

      use ice_kinds_mod
      use ice_constants, only: bignum, rad_to_deg, secday
      use ice_communicate, only: my_task, master_task
      use ice_fileunits, only: nu_diag, nu_diag_l

      implicit none
      private
      public :: init_transport, transport_remap, transport_upwind, lagr_tracking
      save

      character (len=char_len), public ::     &
         advection   ! type of advection scheme used
                     ! 'upwind' => 1st order donor cell scheme
                     ! 'remap' => remapping scheme

      logical, parameter :: & ! if true, prescribe area flux across each edge  
         l_fixed_area = .false.

! NOTE: For remapping, hice and hsno are considered tracers.
!       ntrace is not equal to ntrcr!

      integer (kind=int_kind) ::                      &
         ntrace              ! number of tracers in use
                          
      integer (kind=int_kind), dimension(:), allocatable ::             &
         tracer_type       ,&! = 1, 2, or 3 (depends on 0, 1 or 2 other tracers)
         depend              ! tracer dependencies (see below)

      logical (kind=log_kind), dimension (:), allocatable ::             &
         has_dependents      ! true if a tracer has dependent tracers

      integer (kind=int_kind), parameter ::                      &
         integral_order = 3   ! polynomial order of quadrature integrals
                              ! linear=1, quadratic=2, cubic=3

      logical (kind=log_kind), parameter ::     &
         l_dp_midpt = .true.  ! if true, find departure points using
                              ! corrected midpoint velocity


!==================
! Lagrangian
!==================

     !PARAMETERS
      integer (kind=int_kind), parameter :: LAGR_BUFFER_SIZE_PARAM = 100000
      integer (kind=int_kind), parameter :: LAGR_BNDY_SIZE_PARAM   = 400


      type, public :: lagr_point
          real (kind=dbl_kind) ::      &
              lat_init, lon_init, time_init

          integer (kind=int_kind) ::   &
              type, stat, life, max_life

          real (kind=dbl_kind), dimension(2) ::       &
              latlon,     gpos,     bpos,     lpos,   &
              latlon_old, gpos_old, bpos_old

          real (kind=dbl_kind), dimension(0:2,0:2) :: &
              hte, htn, dxt, dyt, dxu, dyu

          real (kind=dbl_kind), dimension(0:1,0:1) :: &
              dsx, dsy

          integer (kind=int_kind) ::  &
              bid

          logical (kind=log_kind) :: migrate

      end type lagr_point
            

      real (kind=dbl_kind), parameter :: lagr_aice_thres = 0.05d0

 
      integer (kind=int_kind), parameter :: lagr_buffer_size = LAGR_BUFFER_SIZE_PARAM

      integer (kind=int_kind), parameter :: lagr_bndy_size = LAGR_BNDY_SIZE_PARAM
      integer (kind=int_kind), parameter :: lagr_bndy_info_len = 5

      integer (kind=int_kind), parameter ::  &
          LAGR_PHYSICAL_BUOY      =  1,      &
          LAGR_VIRTUAL_BUOY       =  2,      &
          LAGR_MIGRATING_BUOY     = -1,      &
          LAGR_MIGRATING_BUOY_ERR = -2

      integer (kind=int_kind), parameter ::  &
          LAGR_ACTIVE       =   1,           &
          LAGR_CLAIMED      =   0,           &
          LAGR_INACTIVE     = - 1,           &
          LAGR_ONLAND       = -11,           &
          LAGR_DEACTIVATED  = -12,           & 
          LAGR_MELTED       = -13,           &
          LAGR_RIDGED       = -14,           &
          LAGR_MIGRATED     = -15,           &
          LAGR_MIGRATED_ERR = -16,           &
          LAGR_XXX          = -99

      integer (kind=int_kind), parameter ::  &
          LAGR_ACTIVATION_LOC   = 0

      real(kind=dbl_kind), parameter ::  &
          LAGR_ACTIVATION_DNSTY = 1.0d0   ! Every T-cell

     !Activating every day, deactivating every 10 days
      real (kind=dbl_kind), parameter ::                            &
          LAGR_ACTIVATION_INTERVAL_DAILY = secday,                  &
          LAGR_ACTIVATION_INTERVAL_3DAY  = secday*3,                &
          LAGR_ACTIVATION_INTERVAL_5DAY  = secday*5,                &
          LAGR_ACTIVATION_INTERVAL_10DAY = secday*10,               &
          LAGR_ACTIVATION_INTERVAL_30DAY = secday*30,               &
          LAGR_ACTIVATION_INTERVAL_1YEAR = secday*365,              &
          LAGR_VIRTUAL_BUOY_MAX_LIFE_DURATION_10DAY  = secday*10,   &
          LAGR_VIRTUAL_BUOY_MAX_LIFE_DURATION_30DAY  = secday*30,   &
          LAGR_VIRTUAL_BUOY_MAX_LIFE_DURATION_1YEAR  = secday*365,  &
          LAGR_VIRTUAL_BUOY_MAX_LIFE_DURATION_5YEAR  = secday*1825, &
          LAGR_VIRTUAL_BUOY_MAX_LIFE_DURATION_10YEAR = secday*3650, &
          LAGR_VIRTUAL_BUOY_MAX_LIFE_DURATION_UNLIM  = bignum,      &
          LAGR_ACTIVATION_INTERVAL = LAGR_ACTIVATION_INTERVAL_DAILY,&! Every day
          LAGR_VIRTUAL_BUOY_MAX_LIFE_DURATION = LAGR_ACTIVATION_INTERVAL_30DAY


      integer (kind=int_kind) :: LAGR_ACTIVATION_FREQ 

      integer (kind=int_kind) :: LAGR_VIRTUAL_BUOY_MAX_LIFE


      integer (kind=int_kind), parameter ::  &
          LAGR_REMAPPING = 1,                &
          LAGR_TRACMASS  = 2,                &
          lagr_scheme = lagr_remapping

      real (kind=dbl_kind), parameter ::     &
          LAGR_REPORT_INTERVAL = 3600.0*6    ! Every 6-hour

      integer (kind=int_kind) ::  &
          LAGR_REPORT_FREQ

      logical (kind=log_kind), parameter ::  &
          lagr_debug   =  .true.,            &
          lagr_logging =  .true.


      character (len=20), parameter ::                          &
          LAGR_ACTIVATION_STR            = 'LAGR-ACTIVATE   @+  ',  &
          LAGR_DEACTIVATION_STR          = 'LAGR-DEACTIVATE @X  ',  &
          LAGR_DEACTIVATION_MELT_STR     = 'LAGR-MELTED     @$  ',  &
          LAGR_DEACTIVATION_LAND_STR     = 'LAGR-LANDED     @#  ',  &
          LAGR_DEACTIVATION_RIDG_STR     = 'LAGR-RIDGED     @%  ',  &
          LAGR_DEACTIVATION_MIGR_STR     = 'LAGR-MIGRATED   >>  ',  &
          LAGR_DEACTIVATION_MIGR_STR_ERR = 'LAGR-MIGRATED   >X  ',  &
          LAGR_DEACTIVATION_UNKN_STR     = 'LAGR-DEACTIVATE ??  ',  &
          LAGR_MIGRATION_STR             = 'LAGR-MIGRATION  <<  ',  &
          LAGR_MIGRATION_STR_ERR         = 'LAGR-MIGRATION  X<  ',  &
          LAGR_REPORT_STR                = 'LAGR-REPORT     @@  ',  &
          LAGR_ABNORMAL_STR              = 'LAGR-ABNORMAL   !!  '
              


      real (kind=dbl_kind) :: lagr_dt

      integer (kind=int_kind) :: lagr_step

      integer (kind=int_kind) :: lagr_active_count
      logical (kind=log_kind), dimension(:), allocatable :: lagr_slot_active 

      type(lagr_point), dimension(:), allocatable :: lagr_points

      type(lagr_point) :: lagr_point_temp

      real (kind=dbl_kind), dimension(:,:,:,:,:), allocatable :: lagr_bndy
      integer (kind=int_kind), dimension(:,:,:), allocatable :: lagr_bndy_count

      real (kind=dbl_kind), parameter :: lagr_bndy_defaultval = bignum


                          
!=======================================================================

      contains

!=======================================================================
!
! This subroutine is a wrapper for init_remap, which initializes the
! remapping transport scheme.  If the model is run with upwind
! transport, no initializations are necessary.
!
! authors William H. Lipscomb, LANL

      subroutine init_transport

      use ice_state, only: ntrcr, trcr_depend, nt_Tsfc, nt_qice, nt_qsno, &
          nt_sice, nt_fbri, nt_iage, nt_FY, nt_alvl, nt_vlvl, &
          nt_apnd, nt_hpnd, nt_ipnd, nt_bgc_n_sk
      use ice_exit, only: abort_ice
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_advect
      use ice_transport_remap, only: init_remap
      use ice_calendar, only: dt_dyn

      use ice_blocks, only: nx_block, ny_block
      use ice_domain_size, only: max_blocks

      integer (kind=int_kind) ::       &
         k, nt, nt1     ! tracer indices

      integer (kind=int_kind) :: li, alloc_error

      call ice_timer_start(timer_advect)  ! advection 

      ntrace = 2 + ntrcr ! hice,hsno,trcr

      if (allocated(tracer_type)) deallocate(tracer_type)
      if (allocated(depend)) deallocate(depend)
      if (allocated(has_dependents)) deallocate(has_dependents)

      allocate (tracer_type   (ntrace), &
                depend        (ntrace), &
                has_dependents(ntrace))

         ! define tracer dependency arrays
         ! see comments in remapping routine

          depend(1:2)         = 0 ! hice, hsno
          tracer_type(1:2)    = 1 ! no dependency
      
          k = 2

          do nt = 1, ntrcr
             depend(k+nt) = trcr_depend(nt) ! 0 for ice area tracers
                                            ! 1 for ice volume tracers
                                            ! 2 for snow volume tracers
             tracer_type(k+nt) = 2          ! depends on 1 other tracer
             if (trcr_depend(nt) == 0) then
                tracer_type(k+nt) = 1       ! depends on no other tracers
             elseif (trcr_depend(nt) > 2) then
                if (trcr_depend(trcr_depend(nt)-2) > 0) then
                   tracer_type(k+nt) = 3    ! depends on 2 other tracers
                endif
             endif
          enddo

          has_dependents = .false.
          do nt = 1, ntrace
             if (depend(nt) > 0) then
                nt1 = depend(nt)
                has_dependents(nt1) = .true.
                if (nt1 > nt) then
                   write(nu_diag,*)     &
                      'Tracer nt2 =',nt,' depends on tracer nt1 =',nt1
                   call abort_ice       &
                      ('ice: remap transport: Must have nt2 > nt1')
                endif
             endif
          enddo                 ! ntrace

          ! diagnostic output
          if (my_task == master_task) then
          write (nu_diag, *) 'tracer        index      depend        type has_dependents'
             nt = 1
                write(nu_diag,*) '   hi  ',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             nt = 2
                write(nu_diag,*) '   hs  ',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
          k=2
          do nt = k+1, k+ntrcr
             if (nt-k==nt_Tsfc) &
                write(nu_diag,*) 'nt_Tsfc',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             if (nt-k==nt_qice) &
                write(nu_diag,*) 'nt_qice',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             if (nt-k==nt_qsno) &
                write(nu_diag,*) 'nt_qsno',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             if (nt-k==nt_sice) &
                write(nu_diag,*) 'nt_sice',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             if (nt-k==nt_fbri) &
                write(nu_diag,*) 'nt_fbri',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             if (nt-k==nt_iage) &
                write(nu_diag,*) 'nt_iage',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             if (nt-k==nt_FY) &
                write(nu_diag,*) 'nt_FY  ',  nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             if (nt-k==nt_alvl) &
                write(nu_diag,*) 'nt_alvl',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             if (nt-k==nt_vlvl) &
                write(nu_diag,*) 'nt_vlvl',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             if (nt-k==nt_apnd) &
                write(nu_diag,*) 'nt_apnd',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             if (nt-k==nt_hpnd) &
                write(nu_diag,*) 'nt_hpnd',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             if (nt-k==nt_ipnd) &
                write(nu_diag,*) 'nt_ipnd',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
             if (nt-k==nt_bgc_N_sk) &
                write(nu_diag,*) 'nt_bgc_sk',nt,depend(nt),tracer_type(nt),&
                                              has_dependents(nt)
          enddo
          endif ! master_task

          if (trim(advection)=='remap') call init_remap    ! grid quantities

! Lagr initialization

          allocate(lagr_bndy_count(nx_block,ny_block,max_blocks),STAT=alloc_error)
          allocate(lagr_bndy(      nx_block,ny_block,                  &
                                   lagr_bndy_info_len,lagr_bndy_size,  &
                                                     max_blocks),STAT=alloc_error)
          allocate(lagr_points(lagr_buffer_size),STAT=alloc_error)
          allocate(lagr_slot_active(lagr_buffer_size),STAT=alloc_error)

          lagr_step = 0
          lagr_dt = dt_dyn    ! from ice_calender

          lagr_active_count = 0
          lagr_slot_active(:) = .false.
          do li = 1, lagr_buffer_size
             lagr_points(li)%stat = LAGR_INACTIVE
          end do  !li

          lagr_point_temp%stat = LAGR_INACTIVE

      call ice_timer_stop(timer_advect)  ! advection 

      end subroutine init_transport


!=======================================================================
!
! This subroutine solves the transport equations for one timestep
! using the conservative remapping scheme developed by John Dukowicz
! and John Baumgardner (DB) and modified for sea ice by William
! Lipscomb and Elizabeth Hunke.
!
! This scheme preserves monotonicity of ice area and tracers.  That is,
! it does not produce new extrema.  It is second-order accurate in space,
! except where gradients are limited to preserve monotonicity. 
!
! authors William H. Lipscomb, LANL

      subroutine transport_remap (dt)

      use ice_blocks, only: nx_block, ny_block
      use ice_boundary, only: ice_HaloUpdate
      use ice_constants, only: c0, &
          field_loc_center, field_loc_NEcorner, &
          field_type_scalar, field_type_vector
      use ice_global_reductions, only: global_sum, global_sum_prod
      use ice_domain, only: nblocks, distrb_info, blocks_ice, halo_info
      use ice_domain_size, only: ncat, max_blocks
      use ice_blocks, only: nx_block, ny_block, block, get_block, nghost
      use ice_state, only: aice0, aicen, vicen, vsnon, trcrn, ntrcr, &
          uvel, vvel, bound_state
      use ice_grid, only: tarea, HTE, HTN
      use ice_exit, only: abort_ice
      use ice_calendar, only: istep1
      use ice_timers, only: ice_timer_start, ice_timer_stop, &
          timer_advect, timer_bound
      use ice_transport_remap, only: horizontal_remap, make_masks, dpx, dpy, dsx, dsy

      real (kind=dbl_kind), intent(in) ::     &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) ::     &
         iblk           ,&! block index
         ilo,ihi,jlo,jhi,&! beginning and end of physical domain
         n              ,&! ice category index
         nt, nt1, nt2     ! tracer indices

      real (kind=dbl_kind),      &
         dimension (nx_block,ny_block,0:ncat,max_blocks) ::     &
         aim            ,&! mean ice category areas in each grid cell
         aimask           ! = 1. if ice is present, = 0. otherwise

      real (kind=dbl_kind),      &
         dimension (nx_block,ny_block,ntrace,ncat,max_blocks) ::     &
         trm            ,&! mean tracer values in each grid cell
         trmask           ! = 1. if tracer is present, = 0. otherwise

      logical (kind=log_kind) ::     &
         l_stop           ! if true, abort the model

      integer (kind=int_kind) ::     &
         istop, jstop     ! indices of grid cell where model aborts 

      integer (kind=int_kind), dimension(0:ncat,max_blocks) ::     &
         icellsnc         ! number of cells with ice

      integer (kind=int_kind),      &
         dimension(nx_block*ny_block,0:ncat,max_blocks) ::     &
         indxinc, indxjnc   ! compressed i/j indices

      type (block) :: &
         this_block           ! block information for current block
      
      ! variables related to optional bug checks

      logical (kind=log_kind), parameter ::     &
         l_conservation_check = .false. ,&! if true, check conservation
         l_monotonicity_check = .false.   ! if true, check monotonicity

      real (kind=dbl_kind), dimension(0:ncat) ::     &
         asum_init      ,&! initial global ice area
         asum_final       ! final global ice area

      real (kind=dbl_kind), dimension(ntrace,ncat) ::     &
         atsum_init     ,&! initial global ice area*tracer
         atsum_final      ! final global ice area*tracer

      real (kind=dbl_kind), dimension (:,:,:,:,:), allocatable ::     &
         tmin         ,&! local min tracer
         tmax           ! local max tracer

      integer (kind=int_kind) :: alloc_error

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1

      call ice_timer_start(timer_advect)  ! advection 

!---!-------------------------------------------------------------------
!---! Prepare for remapping.
!---! Initialize, update ghost cells, fill tracer arrays.
!---!-------------------------------------------------------------------

      l_stop = .false.
      istop = 0
      jstop = 0

    !-------------------------------------------------------------------
    ! Compute open water area in each grid cell.
    ! Note: An aggregate_area call is needed only if the open
    !       water area has changed since the previous call.
    !       Here we assume that aice0 is up to date.
    !-------------------------------------------------------------------

!      !$OMP PARALLEL DO PRIVATE(iblk)
!      do iblk = 1, nblocks
!         call aggregate_area (nx_block, ny_block,
!                              iblk,     &
!                              aicen(:,:,:,iblk),     &
!                              aice (:,:,  iblk),     &
!                              aice0(:,:,  iblk)) 
!      enddo
!      !$OMP END PARALLEL DO

    !-------------------------------------------------------------------
    ! Ghost cell updates for state variables.
    ! Commented out because ghost cells are updated after cleanup_itd.
    !-------------------------------------------------------------------
!      call ice_timer_start(timer_bound)

!      call ice_HaloUpdate (aice0,            halo_info,     &
!                           field_loc_center, field_type_scalar)

!      call bound_state (aicen, trcrn,     &
!                        vicen, vsnon)

!      call ice_timer_stop(timer_bound)

    !-------------------------------------------------------------------
    ! Ghost cell updates for ice velocity.
    ! Commented out because ghost cell velocities are computed
    !  in ice_dyn_evp.
    !-------------------------------------------------------------------

!      call ice_timer_start(timer_bound)
!      call ice_HaloUpdate (uvel,               halo_info,     &
!                           field_loc_NEcorner, field_type_vector)
!      call ice_HaloUpdate (vvel,               halo_info,     &
!                           field_loc_NEcorner, field_type_vector)
!      call ice_timer_stop(timer_bound)


      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks

    !-------------------------------------------------------------------
    ! Fill arrays with fields to be remapped.
    !-------------------------------------------------------------------

         call state_to_tracers(nx_block,          ny_block,             &
                               ntrcr,             ntrace,               &
                               aice0(:,:,  iblk), aicen(:,:,:,iblk),    &
                               trcrn(:,:,1:ntrcr,:,iblk),               &
                               vicen(:,:,:,iblk), vsnon(:,:,  :,iblk),  &
                               aim  (:,:,:,iblk), trm  (:,:,:,:,iblk))

      enddo
      !$OMP END PARALLEL DO

!---!-------------------------------------------------------------------
!---! Optional conservation and monotonicity checks.
!---!-------------------------------------------------------------------

      if (l_conservation_check) then

    !-------------------------------------------------------------------
    ! Compute initial values of globally conserved quantities.
    !-------------------------------------------------------------------

         do n = 0, ncat
            asum_init(n) = global_sum(aim(:,:,n,:),     distrb_info,       &
                                      field_loc_center, tarea)
         enddo

         do n = 1, ncat
            do nt = 1, ntrace
               if (tracer_type(nt)==1) then ! does not depend on another tracer
                  atsum_init(nt,n) =      &
                      global_sum_prod(trm(:,:,nt,n,:), aim(:,:,n,:),       &
                                      distrb_info,     field_loc_center,   &
                                      tarea)
               elseif (tracer_type(nt)==2) then ! depends on another tracer
                  nt1 = depend(nt)
                  work1(:,:,:) = trm(:,:,nt,n,:)*trm(:,:,nt1,n,:)
                  atsum_init(nt,n) =     &
                      global_sum_prod(work1(:,:,:), aim(:,:,n,:),          &
                                      distrb_info,  field_loc_center,      &
                                      tarea)
               elseif (tracer_type(nt)==3) then ! depends on two tracers
                  nt1 = depend(nt)
                  nt2 = depend(nt1)
                  work1(:,:,:) = trm(:,:,nt,n,:)*trm(:,:,nt1,n,:)          &
                                                *trm(:,:,nt2,n,:)
                  atsum_init(nt,n) =     &
                      global_sum_prod(work1(:,:,:), aim(:,:,n,:),          &
                                      distrb_info,  field_loc_center,      &
                                      tarea)
               endif            ! tracer_type
            enddo               ! nt
         enddo                  ! n

      endif                     ! l_conservation_check
      
      if (l_monotonicity_check) then

         allocate(tmin(nx_block,ny_block,ntrace,ncat,max_blocks),     &
                  tmax(nx_block,ny_block,ntrace,ncat,max_blocks),     &
                  STAT=alloc_error)

         if (alloc_error /= 0)      &
              call abort_ice ('ice: allocation error')

         tmin(:,:,:,:,:) = c0
         tmax(:,:,:,:,:) = c0

         !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block,n)
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)         
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

    !------------------------------------------------------------------- 
    ! Compute masks.
    ! Masks are used to prevent tracer values in cells without ice
    !  from being used in the monotonicity check.
    !------------------------------------------------------------------- 

            call make_masks (nx_block,          ny_block,              &
                             ilo, ihi,          jlo, jhi,              &
                             nghost,            ntrace,                &
                             has_dependents,                           &
                             icellsnc(:,iblk),                         &
                             indxinc(:,:,iblk), indxjnc(:,:,iblk),     &
                             aim(:,:,:,iblk),   aimask(:,:,:,iblk),    &
                             trm(:,:,:,:,iblk), trmask(:,:,:,:,iblk))

    !-------------------------------------------------------------------
    ! Compute local max and min of tracer fields.
    !-------------------------------------------------------------------

            do n = 1, ncat
               call local_max_min                                      &  
                            (nx_block,           ny_block,             &
                             ilo, ihi,           jlo, jhi,             &
                             trm (:,:,:,n,iblk),                       &
                             tmin(:,:,:,n,iblk), tmax  (:,:,:,n,iblk), &
                             aimask(:,:,n,iblk), trmask(:,:,:,n,iblk))
            enddo
         enddo
         !$OMP END PARALLEL DO

         call ice_timer_start(timer_bound)
         call ice_HaloUpdate (tmin,             halo_info,     &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (tmax,             halo_info,     &
                              field_loc_center, field_type_scalar)
         call ice_timer_stop(timer_bound)

         !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block,n)
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)         
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do n = 1, ncat
               call quasilocal_max_min (nx_block, ny_block,     &
                                        ilo, ihi, jlo, jhi,     &
                                        tmin(:,:,:,n,iblk),      &
                                        tmax(:,:,:,n,iblk))
            enddo
         enddo
         !$OMP END PARALLEL DO

      endif                     ! l_monotonicity_check

    !-------------------------------------------------------------------
    ! Main remapping routine: Step ice area and tracers forward in time.
    !-------------------------------------------------------------------
   
         call horizontal_remap (dt,                ntrace,             &
                                uvel      (:,:,:), vvel      (:,:,:),  &
                                dpx       (:,:,:), dpy       (:,:,:),  &
                                dsx       (:,:,:), dsy       (:,:,:),  &
                                aim     (:,:,:,:), trm   (:,:,:,:,:),  &
                                l_fixed_area,                          &
                                tracer_type,       depend,             &
                                has_dependents,    integral_order,     &
                                l_dp_midpt)
         
    !-------------------------------------------------------------------
    ! Given new fields, recompute state variables.
    !-------------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks

         call tracers_to_state (nx_block,          ny_block,            &
                                ntrcr,             ntrace,              &
                                aim  (:,:,:,iblk), trm  (:,:,:,:,iblk), &
                                aice0(:,:,  iblk), aicen(:,:,:,iblk),   &
                                trcrn(:,:,1:ntrcr,:,iblk),              &
                                vicen(:,:,:,iblk), vsnon(:,:,  :,iblk))

      enddo                     ! iblk
      !$OMP END PARALLEL DO

    !-------------------------------------------------------------------
    ! Ghost cell updates for state variables.
    !-------------------------------------------------------------------

      call ice_timer_start(timer_bound)

      call bound_state (aicen, trcrn,     &
                        vicen, vsnon)

      call ice_timer_stop(timer_bound)

!---!-------------------------------------------------------------------
!---! Optional conservation and monotonicity checks
!---!-------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! Compute final values of globally conserved quantities.
    ! Check global conservation of area and area*tracers.  (Optional)
    !-------------------------------------------------------------------

      if (l_conservation_check) then

         do n = 0, ncat
            asum_final(n) = global_sum(aim(:,:,n,:),     distrb_info,      &
                                       field_loc_center, tarea)
         enddo

         do n = 1, ncat
            do nt = 1, ntrace
               if (tracer_type(nt)==1) then ! does not depend on another tracer
                  atsum_final(nt,n) =      &
                      global_sum_prod(trm(:,:,nt,n,:), aim(:,:,n,:),       &
                                      distrb_info,     field_loc_center,   &
                                      tarea)
               elseif (tracer_type(nt)==2) then ! depends on another tracer
                  nt1 = depend(nt)
                  work1(:,:,:) = trm(:,:,nt,n,:)*trm(:,:,nt1,n,:)
                  atsum_final(nt,n) =     &
                      global_sum_prod(work1(:,:,:), aim(:,:,n,:),          &
                                      distrb_info,  field_loc_center,      &
                                      tarea)
               elseif (tracer_type(nt)==3) then ! depends on two tracers
                  nt1 = depend(nt)
                  nt2 = depend(nt1)
                  work1(:,:,:) = trm(:,:,nt,n,:)*trm(:,:,nt1,n,:)          &
                                                *trm(:,:,nt2,n,:)
                  atsum_final(nt,n) =     &
                      global_sum_prod(work1(:,:,:), aim(:,:,n,:),          &
                                      distrb_info,  field_loc_center,      &
                                      tarea)
               endif            ! tracer_type
            enddo               ! nt
         enddo                  ! n

         if (my_task == master_task) then
            call global_conservation (l_stop,     &
                                      asum_init(0), asum_final(0))

            if (l_stop) then
               write (nu_diag,*) 'istep1, my_task, iblk =',     &
                                  istep1, my_task, iblk
               write (nu_diag,*) 'transport: conservation error, cat 0'
               call abort_ice('ice remap transport: conservation error')
            endif

            do n = 1, ncat               
               call global_conservation                                 &
                                     (l_stop,                           &
                                      asum_init(n),    asum_final(n),   &
                                      atsum_init(:,n), atsum_final(:,n))

               if (l_stop) then
                  write (nu_diag,*) 'istep1, my_task, iblk, cat =',     &
                                     istep1, my_task, iblk, n
                  write (nu_diag,*) 'transport: conservation error, cat ',n
                  call abort_ice     &
                       ('ice remap transport: conservation error')
               endif
            enddo               ! n

         endif                  ! my_task = master_task

      endif                     ! l_conservation_check

    !-------------------------------------------------------------------
    ! Check tracer monotonicity.  (Optional)
    !-------------------------------------------------------------------

      if (l_monotonicity_check) then
         !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block,n,l_stop,istop,jstop)
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)         
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            l_stop = .false.
            istop = 0
            jstop = 0

            do n = 1, ncat
               call check_monotonicity      &
                               (nx_block,           ny_block,     &
                                ilo, ihi, jlo, jhi,     &
                                tmin(:,:,:,n,iblk), tmax(:,:,:,n,iblk),  &
                                aim (:,:,  n,iblk), trm (:,:,:,n,iblk),  &
                                l_stop,     &
                                istop,              jstop)

               if (l_stop) then
                  write (nu_diag,*) 'istep1, my_task, iblk, cat =',     &
                                     istep1, my_task, iblk, n
                  call abort_ice('ice remap transport: monotonicity error')
               endif
            enddo               ! n

         enddo                  ! iblk
         !$OMP END PARALLEL DO

         deallocate(tmin, tmax, STAT=alloc_error)
         if (alloc_error /= 0) call abort_ice ('deallocation error')

      endif                     ! l_monotonicity_check

      call ice_timer_stop(timer_advect)  ! advection 
           
      end subroutine transport_remap

!=======================================================================
!
! Computes the transport equations for one timestep using upwind. Sets
! several fields into a work array and passes it to upwind routine.

      subroutine transport_upwind (dt)

      use ice_boundary, only: ice_HaloUpdate
      use ice_blocks, only: nx_block, ny_block, block, get_block, nx_block, ny_block
      use ice_constants, only: p5, &
          field_loc_Nface, field_loc_Eface, field_type_vector
      use ice_domain, only: blocks_ice, halo_info, nblocks
      use ice_domain_size, only: ncat, max_blocks
      use ice_state, only: aice0, aicen, vicen, vsnon, trcrn, ntrcr, &
          uvel, vvel, trcr_depend, bound_state
      use ice_grid, only: HTE, HTN, tarea
      use ice_timers, only: ice_timer_start, ice_timer_stop, &
          timer_bound, timer_advect

      real (kind=dbl_kind), intent(in) ::     &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) ::     &
         narr               ! max number of state variable arrays

      integer (kind=int_kind) ::     &
         i, j, iblk       ,&! horizontal indices
         ilo,ihi,jlo,jhi    ! beginning and end of physical domain

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblocks) ::     &
         uee, vnn           ! cell edge velocities

      real (kind=dbl_kind),     &
         dimension (:,:,:,:), allocatable :: &
         works              ! work array

      type (block) ::     &
         this_block           ! block information for current block

      call ice_timer_start(timer_advect)  ! advection 

      narr = 1 + ncat*(3+ntrcr) ! max number of state variable arrays

      allocate (works(nx_block,ny_block,narr,max_blocks))

    !-------------------------------------------------------------------
    ! Get ghost cell values of state variables.
    ! (Assume velocities are already known for ghost cells, also.)
    !-------------------------------------------------------------------
!      call bound_state (aicen, trcrn,     &
!                        vicen, vsnon)

    !-------------------------------------------------------------------
    ! Average corner velocities to edges.
    !-------------------------------------------------------------------
      
      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            uee(i,j,iblk) = p5*(uvel(i,j,iblk) + uvel(i,j-1,iblk))
            vnn(i,j,iblk) = p5*(vvel(i,j,iblk) + vvel(i-1,j,iblk))
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (uee,             halo_info,     &
                           field_loc_Eface, field_type_vector)
      call ice_HaloUpdate (vnn,             halo_info,     &
                           field_loc_Nface, field_type_vector)
      call ice_timer_stop(timer_bound)

      !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi


      !-----------------------------------------------------------------
      ! fill work arrays with fields to be advected
      !-----------------------------------------------------------------

         call state_to_work (nx_block,             ny_block,             &
                             ntrcr,                                      &
                             narr,                 trcr_depend,          &
                             aicen (:,:,  :,iblk), trcrn (:,:,1:ntrcr,:,iblk), &
                             vicen (:,:,  :,iblk), vsnon (:,:,  :,iblk), &
                             aice0 (:,:,    iblk), works (:,:,  :,iblk))

      !-----------------------------------------------------------------
      ! advect
      !-----------------------------------------------------------------

         call upwind_field (nx_block,       ny_block,               &
                            ilo, ihi,       jlo, jhi,               &
                            dt,                                     &
                            narr,           works(:,:,:,iblk),      &
                            uee(:,:,iblk),  vnn    (:,:,iblk),      &
                            HTE(:,:,iblk),  HTN    (:,:,iblk),      &
                            tarea(:,:,iblk))

      !-----------------------------------------------------------------
      ! convert work arrays back to state variables
      !-----------------------------------------------------------------

         call work_to_state (nx_block,            ny_block,             &
                             ntrcr,                                     &
                             narr,                trcr_depend,          &
                             aicen(:,:,  :,iblk), trcrn (:,:,1:ntrcr,:,iblk), &
                             vicen(:,:,  :,iblk), vsnon (:,:,  :,iblk), &
                             aice0(:,:,    iblk), works (:,:,  :,iblk)) 

      enddo                     ! iblk
      !$OMP END PARALLEL DO
 
      deallocate (works)

    !-------------------------------------------------------------------
    ! Ghost cell updates for state variables.
    !-------------------------------------------------------------------

      call ice_timer_start(timer_bound)

      call bound_state (aicen, trcrn,     &
                        vicen, vsnon)

      call ice_timer_stop(timer_bound)

      call ice_timer_stop(timer_advect)  ! advection 

      end subroutine transport_upwind

!=======================================================================
! The next few subroutines (through check_monotonicity) are called
! by transport_remap.
!=======================================================================
!
! Fill ice area and tracer arrays.
! Assume that the advected tracers are hicen, hsnon, trcrn, 
!  qicen(1:nilyr), and qsnon(1:nslyr).
! This subroutine must be modified if a different set of tracers
!   is to be transported.  The rule for ordering tracers
!   is that a dependent tracer (such as qice) must have a larger
!   tracer index than the tracer it depends on (i.e., hice).
!
! author William H. Lipscomb, LANL

      subroutine state_to_tracers (nx_block, ny_block,   &
                                   ntrcr,    ntrace,     &
                                   aice0,    aicen,      &
                                   trcrn,                &
                                   vicen,    vsnon,      &
                                   aim,      trm)

      use ice_constants, only: c0, c1, rhos, Lfresh, puny
      use ice_domain_size, only: ncat, nslyr
      use ice_state, only: nt_qsno

      integer (kind=int_kind), intent(in) ::     &
           nx_block, ny_block, & ! block dimensions
           ntrcr             , & ! number of tracers in use
           ntrace                ! number of tracers in use incl. hi, hs

      real (kind=dbl_kind), dimension (nx_block,ny_block),     &
           intent(in) ::     &
           aice0     ! fractional open water area

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),     &
           intent(in) ::     &
           aicen   ,&! fractional ice area
           vicen   ,&! volume per unit area of ice          (m)
           vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat),     &
           intent(in) ::     &
           trcrn     ! ice area tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block,0:ncat),     &
            intent(out)::     &
           aim       ! mean ice area in each grid cell

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrace,ncat),  &
           intent(out) ::     &
           trm       ! mean tracer values in each grid cell

      ! local variables

      integer (kind=int_kind) ::     &
           i, j, n      ,&! standard indices
           it, kt       ,&! tracer indices
           ij             ! combined i/j index

      real (kind=dbl_kind) ::     &
           w1             ! work variable

      integer (kind=int_kind), dimension(nx_block*ny_block,0:ncat) ::  &
           indxi        ,&! compressed i/j indices
           indxj

      integer (kind=int_kind), dimension(0:ncat) ::     &
           icells         ! number of cells with ice

      aim(:,:,0) = aice0(:,:)

      do n = 1, ncat

         trm(:,:,:,n) = c0

    !-------------------------------------------------------------------
    ! Find grid cells where ice is present and fill area array.
    !-------------------------------------------------------------------

         icells(n) = 0
         do j = 1, ny_block
         do i = 1, nx_block
            aim(i,j,n) = aicen(i,j,n)
            if (aim(i,j,n) > puny) then
               icells(n) = icells(n) + 1
               ij = icells(n)
               indxi(ij,n) = i
               indxj(ij,n) = j
            endif               ! aim > puny
         enddo
         enddo
      
    !-------------------------------------------------------------------
    ! Fill tracer array
    ! Note: If aice > 0, then hice > 0, but we can have hsno = 0.
    ! Alse note: We transport qice*nilyr rather than qice, so as to
    !  avoid extra operations here and in tracers_to_state.
    !-------------------------------------------------------------------

         do ij = 1, icells(n)
            i = indxi(ij,n)
            j = indxj(ij,n)
            w1 = c1 / aim(i,j,n)
            trm(i,j,1,n) = vicen(i,j,n) * w1 ! hice
            trm(i,j,2,n) = vsnon(i,j,n) * w1 ! hsno
         enddo
         kt = 2

         do it = 1, ntrcr
            if (it >= nt_qsno .and. it < nt_qsno+nslyr) then
               do ij = 1, icells(n)
                  i = indxi(ij,n)
                  j = indxj(ij,n)
                  trm(i,j,kt+it,n) = trcrn(i,j,it,n) + rhos*Lfresh ! snow enthalpy
               enddo
            else
               do ij = 1, icells(n)
                  i = indxi(ij,n)
                  j = indxj(ij,n)
                  trm(i,j,kt+it,n) = trcrn(i,j,it,n) ! other tracers
               enddo
            endif
         enddo
      enddo                     ! ncat
 
      end subroutine state_to_tracers

!=======================================================================
!
! Convert area and tracer arrays back to state variables.
!
! author William H. Lipscomb, LANL

      subroutine tracers_to_state (nx_block, ny_block,   &
                                   ntrcr,    ntrace,     &
                                   aim,      trm,        &
                                   aice0,    aicen,      &
                                   trcrn,                &
                                   vicen,    vsnon)

      use ice_constants, only: c0, rhos, Lfresh
      use ice_domain_size, only: ncat, nslyr
      use ice_state, only: nt_qsno

      integer (kind=int_kind), intent(in) ::     &
           nx_block, ny_block, & ! block dimensions
           ntrcr             , & ! number of tracers in use
           ntrace                ! number of tracers in use incl. hi, hs

      real (kind=dbl_kind), dimension (nx_block,ny_block,0:ncat),     &
           intent(in) ::     &
           aim       ! fractional ice area

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrace,ncat),  &
           intent(in) ::     &
           trm       ! mean tracer values in each grid cell

      real (kind=dbl_kind), dimension (nx_block,ny_block),     &
           intent(inout) ::     &
           aice0     ! fractional ice area

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),     &
           intent(inout) ::     &
           aicen   ,&! fractional ice area
           vicen   ,&! volume per unit area of ice          (m)
           vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat),  &
           intent(inout) ::     &
           trcrn     ! tracers

      ! local variables

      integer (kind=int_kind) ::     &
           i, j, k, n      ,&! standard indices
           it, kt          ,&! tracer indices
           icells          ,&! number of cells with ice
           ij

      integer (kind=int_kind), dimension (nx_block*ny_block) ::     &
           indxi, indxj      ! compressed indices

      aice0(:,:) = aim(:,:,0)

      do n = 1, ncat

      icells = 0
      do j = 1, ny_block
      do i = 1, nx_block
         if (aim(i,j,n) > c0) then
            icells = icells + 1
            indxi(icells) = i
            indxj(icells) = j
         endif
      enddo
      enddo

    !-------------------------------------------------------------------
    ! Compute state variables.
    !-------------------------------------------------------------------

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            aicen(i,j,n) = aim(i,j,n)
            vicen(i,j,n) = aim(i,j,n)*trm(i,j,1,n) ! aice*hice
            vsnon(i,j,n) = aim(i,j,n)*trm(i,j,2,n) ! aice*hsno
         enddo                  ! ij
         kt = 2

         do it = 1, ntrcr
            if (it >= nt_qsno .and. it < nt_qsno+nslyr) then
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  trcrn(i,j,it,n) = trm(i,j,kt+it,n) - rhos*Lfresh ! snow enthalpy
               enddo
               else
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  trcrn(i,j,it,n) = trm(i,j,kt+it,n)  ! other tracers
               enddo
            endif
         enddo
      enddo                     ! ncat

      end subroutine tracers_to_state

!=======================================================================
!
! Check whether values of conserved quantities have changed.
! An error probably means that ghost cells are treated incorrectly.
!
! author William H. Lipscomb, LANL

      subroutine global_conservation (l_stop,                     &
                                      asum_init,  asum_final,     &
                                      atsum_init, atsum_final)

      use ice_constants, only: puny

      real (kind=dbl_kind), intent(in) ::     &
         asum_init   ,&! initial global ice area
         asum_final    ! final global ice area

      real (kind=dbl_kind), dimension(ntrace), intent(in), optional :: &
         atsum_init  ,&! initial global ice area*tracer
         atsum_final   ! final global ice area*tracer

      logical (kind=log_kind), intent(inout) ::     &
         l_stop    ! if true, abort on return

      ! local variables

      integer (kind=int_kind) ::     &
           nt            ! tracer index

      real (kind=dbl_kind) ::     &
           diff          ! difference between initial and final values


      if (asum_init > puny) then
         diff = asum_final - asum_init
         if (abs(diff/asum_init) > puny) then
            l_stop = .true.
            write (nu_diag,*)
            write (nu_diag,*) 'Ice area conserv error'
            write (nu_diag,*) 'Initial global area =', asum_init
            write (nu_diag,*) 'Final global area =', asum_final
            write (nu_diag,*) 'Fractional error =', abs(diff)/asum_init
            write (nu_diag,*) 'asum_final-asum_init =', diff
         endif
      endif

      if (present(atsum_init)) then
       do nt = 1, ntrace
         if (abs(atsum_init(nt)) > puny) then
            diff = atsum_final(nt) - atsum_init(nt)
            if (abs(diff/atsum_init(nt)) > puny) then
               l_stop = .true.
               write (nu_diag,*)
               write (nu_diag,*) 'area*tracer conserv error'
               write (nu_diag,*) 'tracer index =', nt
               write (nu_diag,*) 'Initial global area*tracer =',   &
                                  atsum_init(nt)
               write (nu_diag,*) 'Final global area*tracer =',     &
                                  atsum_final(nt)
               write (nu_diag,*) 'Fractional error =',             &
                                  abs(diff)/atsum_init(nt)
               write (nu_diag,*) 'atsum_final-atsum_init =', diff
            endif
         endif
       enddo
      endif                     ! present(atsum_init)

      end subroutine global_conservation

!=======================================================================
!
! At each grid point, compute the local max and min of a scalar
! field phi: i.e., the max and min values in the nine-cell region
! consisting of the home cell and its eight neighbors.
! 
! To extend to the neighbors of the neighbors (25 cells in all),
! follow this call with a call to quasilocal_max_min.
!
! author William H. Lipscomb, LANL

      subroutine local_max_min (nx_block, ny_block,     &
                                ilo, ihi, jlo, jhi,     &
                                trm,                    &
                                tmin,     tmax,         &
                                aimask,   trmask)

      use ice_constants, only: c1

      integer (kind=int_kind), intent(in) ::     &
           nx_block, ny_block,&! block dimensions
           ilo,ihi,jlo,jhi     ! beginning and end of physical domain

      real (kind=dbl_kind), intent(in),        &
           dimension(nx_block,ny_block) ::     &
           aimask         ! ice area mask

      real (kind=dbl_kind), intent(in),               &
           dimension (nx_block,ny_block,ntrace) ::    &
           trm          ,&! tracer fields
           trmask         ! tracer mask

      real (kind=dbl_kind), intent(out),              &
           dimension (nx_block,ny_block,ntrace) ::    &
           tmin         ,&! local min tracer
           tmax           ! local max tracer

      ! local variables

      integer (kind=int_kind) ::     &
           i, j         ,&! horizontal indices
           nt, nt1        ! tracer indices

      real (kind=dbl_kind), dimension(nx_block,ny_block) ::     &
           phimask        ! aimask or trmask, as appropriate

      real (kind=dbl_kind) ::     &
           phi_nw, phi_n, phi_ne ,&! field values in 8 neighbor cells
           phi_w, phi_e          ,&
           phi_sw, phi_s, phi_se

      do nt = 1, ntrace

         if (tracer_type(nt)==1) then  ! does not depend on another tracer

            do j = 1, ny_block
            do i = 1, nx_block
               phimask(i,j) = aimask(i,j)
            enddo
            enddo

         else   ! depends on another tracer

            nt1 = depend(nt)
            do j = 1, ny_block
            do i = 1, nx_block
               phimask(i,j) = trmask(i,j,nt1)
            enddo
            enddo

         endif

!-----------------------------------------------------------------------
!  Store values of trm in the 8 neighbor cells.
!  If aimask = 1, use the true value; otherwise use the home cell value
!  so that non-physical values of phi do not contribute to the gradient.
!-----------------------------------------------------------------------

         do j = jlo, jhi
            do i = ilo, ihi

               phi_nw = phimask(i-1,j+1) * trm(i-1,j+1,nt)     &
                  + (c1-phimask(i-1,j+1))* trm(i,  j,  nt)
               phi_n  = phimask(i,  j+1) * trm(i,  j+1,nt)     &
                  + (c1-phimask(i,  j+1))* trm(i,  j,  nt)
               phi_ne = phimask(i+1,j+1) * trm(i+1,j+1,nt)     &
                  + (c1-phimask(i+1,j+1))* trm(i,  j,  nt)
               phi_w  = phimask(i-1,j)   * trm(i-1,j,  nt)     &
                  + (c1-phimask(i-1,j))  * trm(i,  j,  nt)
               phi_e  = phimask(i+1,j)   * trm(i+1,j,  nt)     &
                  + (c1-phimask(i+1,j))  * trm(i,  j,  nt)
               phi_sw = phimask(i-1,j-1) * trm(i-1,j-1,nt)     &
                  + (c1-phimask(i-1,j-1))* trm(i,  j,  nt)
               phi_s  = phimask(i,  j-1) * trm(i,  j-1,nt)     &
                  + (c1-phimask(i,  j-1))* trm(i,  j,  nt)
               phi_se = phimask(i+1,j-1) * trm(i+1,j-1,nt)     &
                  + (c1-phimask(i+1,j-1))* trm(i,  j,  nt)

!-----------------------------------------------------------------------
!     Compute the minimum and maximum among the nine local cells.
!-----------------------------------------------------------------------

               tmax(i,j,nt) = max (phi_nw, phi_n,  phi_ne, phi_w,     &
                      trm(i,j,nt), phi_e,  phi_sw, phi_s,  phi_se)

               tmin(i,j,nt) = min (phi_nw, phi_n,  phi_ne, phi_w,     &
                      trm(i,j,nt), phi_e,  phi_sw, phi_s,  phi_se)

            enddo               ! i
         enddo                  ! j

      enddo                     ! nt

      end subroutine local_max_min

!=======================================================================
!
! Extend the local max and min by one grid cell in each direction.
! Incremental remapping is monotone for the "quasilocal" max and min,
! but in rare cases may violate monotonicity for the local max and min.
!
! author William H. Lipscomb, LANL

      subroutine quasilocal_max_min (nx_block, ny_block,     &
                                     ilo, ihi, jlo, jhi,     &
                                     tmin,     tmax)

      integer (kind=int_kind), intent(in) ::     &
         nx_block, ny_block,&! block dimensions
         ilo,ihi,jlo,jhi     ! beginning and end of physical domain

      real (kind=dbl_kind), intent(inout),     &
           dimension (nx_block,ny_block,ntrace) ::     &
           tmin         ,&! local min tracer
           tmax           ! local max tracer

      ! local variables

      integer (kind=int_kind) ::     &
           i, j          ,&! horizontal indices
           nt              ! tracer index

      do nt = 1, ntrace

         do j = jlo, jhi
         do i = ilo, ihi

            tmax(i,j,nt) =     &
              max (tmax(i-1,j+1,nt), tmax(i,j+1,nt), tmax(i+1,j+1,nt),     &
                   tmax(i-1,j,  nt), tmax(i,j,  nt), tmax(i+1,j,  nt),     &
                   tmax(i-1,j-1,nt), tmax(i,j-1,nt), tmax(i+1,j-1,nt))

            tmin(i,j,nt) =     &
              min (tmin(i-1,j+1,nt), tmin(i,j+1,nt), tmin(i+1,j+1,nt),     &
                   tmin(i-1,j,  nt), tmin(i,j,  nt), tmin(i+1,j,  nt),     &
                   tmin(i-1,j-1,nt), tmin(i,j-1,nt), tmin(i+1,j-1,nt))

         enddo                  ! i
         enddo                  ! j

      enddo

      end subroutine quasilocal_max_min

!======================================================================
!
! At each grid point, make sure that the new tracer values
! fall between the local max and min values before transport.
!
! author William H. Lipscomb, LANL

      subroutine check_monotonicity (nx_block, ny_block,     &
                                     ilo, ihi, jlo, jhi,     &
                                     tmin,     tmax,         &
                                     aim,      trm,          &
                                     l_stop,                 &
                                     istop,    jstop)

      use ice_constants, only: c1, puny

      integer (kind=int_kind), intent(in) ::     &
           nx_block, ny_block,&! block dimensions
           ilo,ihi,jlo,jhi     ! beginning and end of physical domain

      real (kind=dbl_kind), intent(in),         &
           dimension (nx_block,ny_block) ::     &
           aim            ! new ice area

      real (kind=dbl_kind), intent(in),                &
           dimension (nx_block,ny_block,ntrace) ::     &
           trm            ! new tracers

      real (kind=dbl_kind), intent(in),                &
           dimension (nx_block,ny_block,ntrace) ::     &
           tmin         ,&! local min tracer
           tmax           ! local max tracer

      logical (kind=log_kind), intent(inout) ::     &
         l_stop    ! if true, abort on return

      integer (kind=int_kind), intent(inout) ::     &
         istop, jstop     ! indices of grid cell where model aborts 

      ! local variables

      integer (kind=int_kind) ::     &
           i, j           ,&! horizontal indices
           nt, nt1, nt2     ! tracer indices

      real (kind=dbl_kind) ::     &
           w1, w2         ! work variables

      logical (kind=log_kind), dimension (nx_block, ny_block) ::   &
           l_check        ! if true, check monotonicity

      do nt = 1, ntrace

    !-------------------------------------------------------------------
    ! Load logical array to identify tracers that need checking.
    !-------------------------------------------------------------------

         if (tracer_type(nt)==1) then ! does not depend on another tracer

            do j = jlo, jhi
            do i = ilo, ihi
               if (aim(i,j) > puny) then 
                  l_check(i,j) = .true.
               else
                  l_check(i,j) = .false.
               endif
            enddo
            enddo

         elseif (tracer_type(nt)==2) then ! depends on another tracer

            nt1 = depend(nt)
            do j = jlo, jhi
            do i = ilo, ihi
               if (abs(trm(i,j,nt1)) > puny) then
                  l_check(i,j) = .true.
               else
                  l_check(i,j) = .false.
               endif
            enddo
            enddo

         elseif (tracer_type(nt)==3) then ! depends on two tracers

            nt1 = depend(nt)
            nt2 = depend(nt1)
            do j = jlo, jhi
            do i = ilo, ihi
               if (abs(trm(i,j,nt1)) > puny .and.     &
                   abs(trm(i,j,nt2)) > puny) then
                  l_check(i,j) = .true.
               else
                  l_check(i,j) = .false.
               endif
            enddo
            enddo
         endif

    !-------------------------------------------------------------------
    ! Make sure new values lie between tmin and tmax
    !-------------------------------------------------------------------

         do j = jlo, jhi
         do i = ilo, ihi

            if (l_check(i,j)) then
               ! w1 and w2 allow for roundoff error when abs(trm) is big
               w1 = max(c1, abs(tmin(i,j,nt)))
               w2 = max(c1, abs(tmax(i,j,nt)))
               if (trm(i,j,nt) < tmin(i,j,nt)-w1*puny) then
                  l_stop = .true.
                  istop = i
                  jstop = j
                  write (nu_diag,*) ' '
                  write (nu_diag,*) 'new tracer < tmin'
                  write (nu_diag,*) 'i, j, nt =', i, j, nt
                  write (nu_diag,*) 'new tracer =', trm (i,j,nt)
                  write (nu_diag,*) 'tmin ='      , tmin(i,j,nt)
                  write (nu_diag,*) 'ice area ='  , aim(i,j)
               elseif (trm(i,j,nt) > tmax(i,j,nt)+w2*puny) then
                  l_stop = .true.
                  istop = i
                  jstop = j
                  write (nu_diag,*) ' '
                  write (nu_diag,*) 'new tracer > tmax'
                  write (nu_diag,*) 'i, j, nt =', i, j, nt
                  write (nu_diag,*) 'new tracer =', trm (i,j,nt)
                  write (nu_diag,*) 'tmax ='      , tmax(i,j,nt)
                  write (nu_diag,*) 'ice area ='  , aim(i,j)
               endif
            endif

         enddo                  ! i
         enddo                  ! j

      enddo                     ! nt

      end subroutine check_monotonicity

!=======================================================================
! The remaining subroutines are called by transport_upwind.
!=======================================================================
!
! Fill work array with state variables in preparation for upwind transport

      subroutine state_to_work (nx_block, ny_block,        &
                                ntrcr,                     &
                                narr,     trcr_depend,     &
                                aicen,    trcrn,           &
                                vicen,    vsnon,           &
                                aice0,    works)

      use ice_domain_size, only: ncat
      use ice_state, only: nt_alvl, nt_apnd, nt_fbri, &
                           tr_pond_cesm, tr_pond_lvl, tr_pond_topo

      integer (kind=int_kind), intent(in) ::     &
         nx_block, ny_block, & ! block dimensions
         ntrcr             , & ! number of tracers in use
         narr        ! number of 2D state variable arrays in works array

      integer (kind=int_kind), dimension (ntrcr), intent(in) ::     &
         trcr_depend ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),     &
         intent(in) ::     &
         aicen   ,&! concentration of ice
         vicen   ,&! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat),     &
         intent(in) ::     &
         trcrn     ! ice tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block),         &
         intent(in) ::        &
         aice0     ! concentration of open water

      real (kind=dbl_kind), dimension(nx_block,ny_block,narr),     &
         intent (out) ::      &
         works     ! work array

      ! local variables

      integer (kind=int_kind) ::      &
         i, j, n, it    ,&! counting indices
         narrays          ! counter for number of state variable arrays

      !-----------------------------------------------------------------
      ! This array is used for performance (balance memory/cache vs
      ! number of bound calls);  a different number of arrays may perform
      ! better depending on the machine used, number of processors, etc.
      ! --tested on SGI R2000, using 4 pes for the ice model under MPI
      !-----------------------------------------------------------------

      do j = 1, ny_block
      do i = 1, nx_block
         works(i,j,1) = aice0(i,j)
      enddo
      enddo
      narrays = 1

      do n=1, ncat

         do j = 1, ny_block
         do i = 1, nx_block
            works(i,j,narrays+1) = aicen(i,j,n)
            works(i,j,narrays+2) = vicen(i,j,n)
            works(i,j,narrays+3) = vsnon(i,j,n)
         enddo                  ! i
         enddo                  ! j
         narrays = narrays + 3

         do it = 1, ntrcr
            if (trcr_depend(it) == 0) then
               do j = 1, ny_block
               do i = 1, nx_block
                  works(i,j,narrays+it) = aicen(i,j,n)*trcrn(i,j,it,n)
               enddo
               enddo
            elseif (trcr_depend(it) == 1) then
               do j = 1, ny_block
               do i = 1, nx_block
                  works(i,j,narrays+it) = vicen(i,j,n)*trcrn(i,j,it,n)
               enddo
               enddo
            elseif (trcr_depend(it) == 2) then
               do j = 1, ny_block
               do i = 1, nx_block
                  works(i,j,narrays+it) = vsnon(i,j,n)*trcrn(i,j,it,n)
               enddo
               enddo
            elseif (trcr_depend(it) == 2+nt_alvl) then
               do j = 1, ny_block
               do i = 1, nx_block
                  works(i,j,narrays+it) = aicen(i,j,n) &
                                        * trcrn(i,j,nt_alvl,n) &
                                        * trcrn(i,j,it,n)
               enddo
               enddo
            elseif (trcr_depend(it) == 2+nt_apnd .and. &
                    tr_pond_cesm .or. tr_pond_topo) then
               do j = 1, ny_block
               do i = 1, nx_block
                  works(i,j,narrays+it) = aicen(i,j,n) &
                                        * trcrn(i,j,nt_apnd,n) &
                                        * trcrn(i,j,it,n)
               enddo
               enddo
            elseif (trcr_depend(it) == 2+nt_apnd .and. &
                    tr_pond_lvl) then
               do j = 1, ny_block
               do i = 1, nx_block
                  works(i,j,narrays+it) = aicen(i,j,n) &
                                        * trcrn(i,j,nt_alvl,n) &
                                        * trcrn(i,j,nt_apnd,n) &
                                        * trcrn(i,j,it,n)
               enddo
               enddo
            elseif (trcr_depend(it) == 2+nt_fbri) then
               do j = 1, ny_block
               do i = 1, nx_block
                  works(i,j,narrays+it) = vicen(i,j,n) &
                                        * trcrn(i,j,nt_fbri,n) &
                                        * trcrn(i,j,it,n)
               enddo
               enddo
            endif
         enddo
         narrays = narrays + ntrcr

      enddo                     ! n

      if (narr /= narrays) write(nu_diag,*)      &
           "Wrong number of arrays in transport bound call"

      end subroutine state_to_work

!=======================================================================
!
! Convert work array back to state variables

      subroutine work_to_state (nx_block, ny_block,        &
                                ntrcr,                     &
                                narr,     trcr_depend,     &
                                aicen,    trcrn,           &
                                vicen,    vsnon,           &
                                aice0,    works)

      use ice_domain_size, only: ncat
      use ice_blocks, only: 
      use ice_itd, only: compute_tracers

      integer (kind=int_kind), intent (in) ::                       &
         nx_block, ny_block, & ! block dimensions
         ntrcr             , & ! number of tracers in use
         narr        ! number of 2D state variable arrays in works array

      integer (kind=int_kind), dimension (ntrcr), intent(in) ::     &
         trcr_depend ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon

      real (kind=dbl_kind), intent (in) ::                          &
         works (nx_block,ny_block,narr)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat),     &
         intent(out) ::     &
         aicen   ,&! concentration of ice
         vicen   ,&! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat), &
         intent(out) ::     &
         trcrn     ! ice tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block),          &
         intent(out) ::     &
         aice0     ! concentration of open water

      ! local variables

      integer (kind=int_kind) ::      &
         i, j, n        ,&! counting indices
         narrays        ,&! counter for number of state variable arrays
         icells           ! number of ocean/ice cells

      integer (kind=int_kind), dimension (nx_block*ny_block) ::        &
        indxi, indxj

      real (kind=dbl_kind), dimension (nx_block*ny_block,narr) ::      &
         work 

      ! for call to compute_tracers
      icells = 0
      do j = 1, ny_block
      do i = 1, nx_block
         icells = icells + 1
         indxi(icells) = i
         indxj(icells) = j
         work (icells,:) = works(i,j,:)
      enddo
      enddo

      do j=1,ny_block
      do i=1,nx_block
         aice0(i,j) = works(i,j,1)
      enddo
      enddo
      narrays = 1               ! aice0 is first array

      do n=1,ncat

         do j=1,ny_block
         do i=1,nx_block
            aicen(i,j,n) = works(i,j,narrays+1)
            vicen(i,j,n) = works(i,j,narrays+2)
            vsnon(i,j,n) = works(i,j,narrays+3)
         enddo
         enddo
         narrays = narrays + 3

         call compute_tracers (nx_block,     ny_block,               &
                               icells,       indxi,   indxj,         &
                               ntrcr,        trcr_depend,            &
                               work (:,narrays+1:narrays+ntrcr),     &
                               aicen(:,:,n),                         &
                               vicen(:,:,n), vsnon(:,:,n),           &
                               trcrn(:,:,:,n))

         narrays = narrays + ntrcr

      enddo                     ! ncat

      end subroutine work_to_state

!=======================================================================
!
! upwind transport algorithm

      subroutine upwind_field (nx_block, ny_block,   &
                               ilo, ihi, jlo, jhi,   &
                               dt,                   &
                               narrays,  phi,        &
                               uee,      vnn,        &
                               HTE,      HTN,        &
                               tarea)

      use ice_constants, only: p5

      integer (kind=int_kind), intent (in) ::     &
         nx_block, ny_block ,&! block dimensions
         ilo,ihi,jlo,jhi    ,&! beginning and end of physical domain
         narrays              ! number of 2D arrays to be transported

      real (kind=dbl_kind), intent(in) ::         &
         dt                   ! time step

      real (kind=dbl_kind), dimension(nx_block,ny_block,narrays), &
         intent(inout) ::                                         &
         phi                  ! scalar field

      real (kind=dbl_kind), dimension(nx_block,ny_block),         &
         intent(in)::     &
         uee, vnn             ! cell edge velocities

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: &
         HTE                ,&! length of east cell edge 
         HTN                ,&! length of north cell edge
         tarea                ! grid cell area

      ! local variables

      integer (kind=int_kind) :: &
         i, j, n              ! standard indices

      real (kind=dbl_kind) :: &
         upwind, y1, y2, a, h ! function

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         worka, workb

    !-------------------------------------------------------------------
    ! Define upwind function
    !-------------------------------------------------------------------

      upwind(y1,y2,a,h) = p5*dt*h*((a+abs(a))*y1+(a-abs(a))*y2)

    !-------------------------------------------------------------------
    ! upwind transport
    !-------------------------------------------------------------------

      do n = 1, narrays

         do j = 1, jhi
         do i = 1, ihi
            worka(i,j)=     &
               upwind(phi(i,j,n),phi(i+1,j,n),uee(i,j),HTE(i,j))
            workb(i,j)=     &
               upwind(phi(i,j,n),phi(i,j+1,n),vnn(i,j),HTN(i,j))
         enddo
         enddo

         do j = jlo, jhi
         do i = ilo, ihi
            phi(i,j,n) = phi(i,j,n) - ( worka(i,j)-worka(i-1,j)      &
                                      + workb(i,j)-workb(i,j-1) )    &
                                      / tarea(i,j)
         enddo
         enddo

      enddo                     ! narrays

      end subroutine upwind_field


!===========================================================
!
! Lagrangian tracking
!
! a. basic indexing: local, block, global
!
!      2   x---------x---------x---------x
!          |         |         |         |
!          |         |         |         |
!          |         |    +----|----+    |
!          |         |    |    |    |    |
!          |         |    |    |    |    |
!      1   x---------x---------x---------x
!          |         |    |    |    |    |
!          |         |    |    |    |    |
!          |         |    +----|----+    |
!          |         |         |         |
!          |         |         |         |
!      0   x---------x---------x---------x
!          |         |         |         |
!          |         |         |         |
!          |         |         |         |
!          |         |         |         |
!          |         |         |         |
!     -1   x---------x---------x---------x
!                                
!         -1         0         1         2
! 
!        lpos-2    lpos-1     lpos     lpos+1
!
!
!
!===========================================================


!========================================!==============================
!
! Entry point, called by dynamics after transport_remapping
!
!=========================================

      subroutine lagr_tracking( dt )

      use ice_blocks, only: block, get_block, nx_block, ny_block, nghost
      use ice_calendar, only: time
      use ice_domain, only: blocks_ice
      use ice_constants, only: field_loc_center, field_type_scalar, secday
      use ice_domain, only: nblocks, halo_info
      use ice_boundary, only: ice_HaloUpdate

      real (kind=dbl_kind), intent(in) :: dt

!
      type (block) :: this_block

      integer (kind=int_kind) :: i, j, iblk, ilgr, cidx, lidx, lagr_m_cnt
      
      real (kind=dbl_kind), dimension(2) :: cgpos, clatlon_init

      real (kind=dbl_kind) :: ctime, clon, clat

      logical (kind=log_kind), parameter :: DBG_INTERIM_STATUS = .false.
      logical (kind=log_kind), parameter :: DBG_LOCAL = .false.

!
      lagr_dt = dt

      LAGR_ACTIVATION_FREQ = nint(LAGR_ACTIVATION_INTERVAL / lagr_dt)
      LAGR_VIRTUAL_BUOY_MAX_LIFE = nint(LAGR_VIRTUAL_BUOY_MAX_LIFE_DURATION / lagr_dt)

      LAGR_REPORT_FREQ = nint(LAGR_REPORT_INTERVAL / lagr_dt)
      

    ! Deactivation (thermodynamics, life limit, etc.)
      do ilgr = 1, lagr_buffer_size
         if (is_lagr_point_active(lagr_points(ilgr))) then
            if (get_lagr_aice(lagr_points(ilgr)) < lagr_aice_thres) then
               call lagr_deactivate(lagr_points(ilgr),LAGR_MELTED)

               lagr_active_count = lagr_active_count -1
               lagr_slot_active(ilgr) = .false.
            elseif (lagr_points(ilgr)%life > lagr_points(ilgr)%max_life) then
               call lagr_deactivate(lagr_points(ilgr),LAGR_DEACTIVATED)

               lagr_active_count = lagr_active_count -1
               lagr_slot_active(ilgr) = .false.
            end if
         end if
      end do  !ilgr 

      lagr_step = lagr_step + 1

    ! Report on status
      if (mod(lagr_step,LAGR_REPORT_FREQ) == 1) then
      do ilgr = 1, lagr_buffer_size
         if (is_lagr_point_active(lagr_points(ilgr))) then
            call lagr_report_status(lagr_points(ilgr),nu_diag_l)
         end if
      end do
      end if

    ! Activation of Virtual Buoys
      if (mod(lagr_step,LAGR_ACTIVATION_FREQ) == 1) then 
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            call lagr_deploy_virtual_buoys(this_block,time-lagr_dt)
         end do
      end if

    ! Activation of Physical Buoys
    ! TODO

    ! Tracking
      lagr_bndy_count(    :,:,:) = 0
      lagr_bndy      (:,:,:,:,:) = lagr_bndy_defaultval

      do ilgr = 1, lagr_buffer_size
         if (.not. is_lagr_point_active(lagr_points(ilgr))) then
            cycle
         end if

        !Getting older
         call lagr_increase_life(lagr_points(ilgr))

        !Do tracking, with updated bpos, gpos, lpos
        !Optional update to local grid vars & latlon if NOT migrating
         call lagr_do_tracking(lagr_points(ilgr))

         if (lagr_debug .and. DBG_INTERIM_STATUS) then
            call lagr_report_status(lagr_points(ilgr),nu_diag_l)
         end if

        !Check for migration, filling boundary exchange fields
         if (lagr_points(ilgr)%migrate) then

           !Source locations in i & j
            i = ceiling(lagr_points(ilgr)%bpos_old(1))
            j = ceiling(lagr_points(ilgr)%bpos_old(2))
            iblk = lagr_points(ilgr)%bid

            if (lagr_bndy_count(i,j,iblk) < lagr_bndy_size) then
               lagr_bndy_count(i,j,iblk) = lagr_bndy_count(i,j,iblk) + 1

               cidx = lagr_bndy_count(i,j,iblk)

               lagr_bndy(i,j,1,cidx,iblk) = lagr_points(ilgr)%lon_init
               lagr_bndy(i,j,2,cidx,iblk) = lagr_points(ilgr)%lat_init

               if (lagr_points(ilgr)%type == LAGR_PHYSICAL_BUOY) then
                  lagr_bndy(i,j,3,cidx,iblk) =  lagr_points(ilgr)%time_init
               elseif (lagr_points(ilgr)%type == LAGR_VIRTUAL_BUOY) then
                  lagr_bndy(i,j,3,cidx,iblk) = -lagr_points(ilgr)%time_init -secday
               end if

               lagr_bndy(i,j,4,cidx,iblk) = lagr_points(ilgr)%gpos(1)
               lagr_bndy(i,j,5,cidx,iblk) = lagr_points(ilgr)%gpos(2)

              !Reclaim 
               call lagr_deactivate(lagr_points(ilgr),LAGR_MIGRATED)
            else
              !Reclaim 
               call lagr_deactivate(lagr_points(ilgr),LAGR_MIGRATED_ERR)
            end if

            lagr_active_count = lagr_active_count -1
            lagr_slot_active(ilgr) = .false.
         end if
      end do  !ilgr

    ! Exchanges
      if (DBG_LOCAL) then
      call ice_HaloUpdate( lagr_bndy_count,  halo_info, &
                           field_loc_center, field_type_scalar )
      end if
      call ice_HaloUpdate( lagr_bndy,        halo_info, &
                           field_loc_center, field_type_scalar )

    ! Maintenance
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         do j = this_block%jlo - nghost, this_block%jhi + nghost
         do i = this_block%ilo - nghost, this_block%ihi + nghost
           !Excluding inner region
            if ((i >= this_block%ilo) .and.  &
                (i <= this_block%ihi) .and.  &
                (j >= this_block%jlo) .and.  &
                (j <= this_block%jhi)) then
               cycle
            end if

            if (DBG_LOCAL) then
            lagr_m_cnt = lagr_bndy_count(i,j,iblk)
            else
            lagr_m_cnt = lagr_bndy_size
            end if

            do ilgr = 1, lagr_m_cnt
               if (lagr_debug .and. DBG_LOCAL) then
                  write(nu_diag_l,*) 'The ', lagr_m_cnt, '-th migrating Lagr. point'
               end if

               clatlon_init(1) = lagr_bndy(i,j,1,ilgr,iblk)
               clatlon_init(2) = lagr_bndy(i,j,2,ilgr,iblk)

               cgpos(1) = lagr_bndy(i,j,4,ilgr,iblk)
               cgpos(2) = lagr_bndy(i,j,5,ilgr,iblk)

               if (lagr_debug .and. DBG_LOCAL) then
                  write(nu_diag_l,*) '  LATLON_INIT:', clatlon_init(1)*rad_to_deg, &
                     clatlon_init(2)*rad_to_deg, &
                     ' ,GPOS:', cgpos(1), cgpos(2)
               end if

              !A valid migrating Lagr point
               if (is_valid_latlon_gpos(clatlon_init,cgpos)) then
                  if (lagr_debug .and. DBG_LOCAL) then
                     write(nu_diag_l,*) 'A migrating Lagr. point ',cgpos(1),cgpos(2), &
                        ', w/ init. latlon: ', clatlon_init(1), clatlon_init(2)
                  end if

                  if (is_lagr_gpos_on_block(cgpos,this_block)) then
                     call lagr_claim_slot(lidx)
                     if (lidx==0) then
                       !write(nu_diag_l,*) 'ERROR! No Lagrangian Slot Available! Losing a PT!'
                        ctime = lagr_bndy(i,j,3,ilgr,iblk)
                        call lagr_activate(                                   &
                                lagr_point_temp,                              &
                                ctime,  LAGR_MIGRATING_BUOY_ERR, this_block,  &
                                cgpos,  latlon_info=clatlon_init )
                        cycle
                     end if

                     ctime = lagr_bndy(i,j,3,ilgr,iblk)
                     call lagr_activate(                               &
                             lagr_points(lidx),                        &
                             ctime,  LAGR_MIGRATING_BUOY, this_block,  &
                             cgpos,  latlon_info=clatlon_init )
                  end if
               end if
            end do

         end do
         end do
      end do

    ! Check on land
      do ilgr = 1, lagr_buffer_size
         if (.not. is_lagr_point_active(lagr_points(ilgr))) then
            cycle
         end if

         if (is_lagr_point_on_land_strict(lagr_points(ilgr))) then
            call lagr_deactivate(lagr_points(ilgr),LAGR_ONLAND)

            lagr_active_count = lagr_active_count -1
            lagr_slot_active(ilgr) = .false.
         end if

      end do

      end subroutine lagr_tracking


!======================================== 
!
! Utilities - Queries
!
!======================================== DONE
!
! Whether the point is active OR not
!
!======================================== 

      function is_lagr_point_active(pt) result(active)

      type (lagr_point), intent(in) :: pt
      logical (kind=log_kind) :: active

!
      active = (pt%stat == LAGR_ACTIVE)

      end function is_lagr_point_active


!======================================== DONE
!
! Whether the Lagr point is on land
!
!========================================

      function is_lagr_point_on_land(pt) result(onland)

      use ice_grid, only: tmask, umask, tlat, tlon

      type (lagr_point), intent(in) :: pt

      logical (kind=log_kind) :: onland

!
      integer (kind=int_kind) :: i, j

!
      i = ceiling(pt%bpos(1))
      j = ceiling(pt%bpos(2))

      if (.not. tmask(i,j,pt%bid)) then
         onland = .true.
      else
         onland = .false.
      end if

      end function is_lagr_point_on_land    


!======================================== DONE
!
! Whether the Lagr point is on land, in strict sense
!
!========================================

      function is_lagr_point_on_land_strict(pt) result(onland)

      use ice_grid, only: tmask, umask, tlat, tlon

      type (lagr_point), intent(in) :: pt

      logical (kind=log_kind) :: onland

!
      integer (kind=int_kind) :: i, j

!
      i = ceiling(pt%bpos(1))
      j = ceiling(pt%bpos(2))

      if (.not. tmask(i,j,pt%bid)) then
         onland = .true.
      elseif ((.not. umask(i  ,j  ,pt%bid)) .or.  &
              (.not. umask(i-1,j  ,pt%bid)) .or.  &
              (.not. umask(i  ,j-1,pt%bid)) .or.  &
              (.not. umask(i-1,j-1,pt%bid)) ) then
         onland = .true.
      else
         onland = .false.
      end if

      end function is_lagr_point_on_land_strict


!======================================== DONE
!
! Whether the new position of the Lagr. point is within current block OR not
!    Checking bpos
!
!======================================== 

      function is_lagr_point_still_in_block(pt) result(inblock)

      use ice_blocks, only: block, get_block
      use ice_domain, only: blocks_ice

      type(lagr_point), intent(in) :: pt
      logical (kind=log_kind) :: inblock

!
      type(block) :: this_block
      integer (kind=int_kind) :: iblk

!
      if (.not. is_lagr_point_active(pt)) then
          inblock = .false.
          return
      end if

      iblk = pt%bid
      this_block = get_block(blocks_ice(iblk),iblk)

      if ((pt%bpos(1) <= this_block%ilo-1) .or.  &
          (pt%bpos(1) >  this_block%ihi  ) .or.  &
          (pt%bpos(2) <= this_block%jlo-1) .or.  &
          (pt%bpos(2) >  this_block%jhi  )) then
          inblock = .false.
      else
          inblock = .true.
      end if
         
      end function is_lagr_point_still_in_block


!======================================== DONE
!
! Whether the global position on current block
!
!======================================== 

      function is_lagr_gpos_on_block( gpos, cblock ) result(onblock)

      use ice_blocks, only: block
 
      type(block), intent(in) :: cblock
      real(kind=dbl_kind), dimension(2), intent(in) :: gpos

      logical(kind=log_kind) :: onblock

!
      integer (kind=int_kind) :: ilo_glob, ihi_glob, jlo_glob, jhi_glob

      logical (kind=log_kind), parameter :: DBG_LOCAL = .false.

!
      ilo_glob = cblock%i_glob(cblock%ilo)
      ihi_glob = cblock%i_glob(cblock%ihi)
      jlo_glob = cblock%j_glob(cblock%jlo)
      jhi_glob = cblock%j_glob(cblock%jhi)

      if ((gpos(1) >  ilo_glob-1) .and.    &
          (gpos(1) <= ihi_glob  ) .and.    &
          (gpos(2) >  jlo_glob-1) .and.    &
          (gpos(2) <= jhi_glob  ) ) then
         onblock = .true.
      else
         onblock = .false.
      end if

      if (lagr_debug .and. DBG_LOCAL) then
         write(nu_diag_l,*) 'I/J_LO/HI_GLO ', ilo_glob, ihi_glob, jlo_glob, jhi_glob, &
            ', GPOS ', gpos(1), gpos(2)
      end if

      end function is_lagr_gpos_on_block


!======================================== 
!
! Is valid LatLon and GPOS
!
!======================================== 

      function is_valid_latlon_gpos( latlon, gpos ) result (isvalid)

      use ice_constants, only: pi
      use ice_domain_size, only: nx_global, ny_global

      real (kind=dbl_kind), dimension(2), intent(in) :: latlon, gpos

      logical (kind=log_kind) :: isvalid

!
      if ((latlon(1) < -pi      ) .or. (latlon(1) > pi      ) .or.  &
          (latlon(2) < -pi/2.0d0) .or. (latlon(2) > pi/2.0d0)) then
         isvalid = .false.
      elseif ((gpos(1) < -1.0d0) .or. (gpos(1) > nx_global+1) .or.  &
              (gpos(2) < -1.0d0) .or. (gpos(2) > ny_global+1)) then
         isvalid = .false.
      else
         isvalid = .true.
      end if

      end function is_valid_latlon_gpos


!======================================== DONE
!
! Sea Ice concentration at current cell
!
!======================================== 

      function get_lagr_aice( pt ) result(caice)

      use ice_state, only: aice, vice, vsno

      type(lagr_point), intent(in) :: pt

      real(kind=dbl_kind) :: caice

!
      integer (kind=int_kind) :: bidxi, bidxj

!
      bidxi = ceiling(pt%bpos(1))
      bidxj = ceiling(pt%bpos(2))

      caice = aice(bidxi,bidxj,pt%bid)

      end function get_lagr_aice


!======================================== 
!
! Report Lagr. point status
!
!======================================== 

      subroutine lagr_report_status (pt,hdl)

      type(lagr_point), intent(in) :: pt

      integer (kind=int_kind), intent(in) :: hdl

!
      call lagr_info_dump(pt,LAGR_REPORT_STR,hdl)

     !TODO: report of sea ice status

      end subroutine lagr_report_status


!======================================== DONE
!
! Sea Ice status of the current cell
!
!======================================== 

      subroutine get_lagr_ice_status( pt, cai, chi, chs, ctsfc )

      use ice_state, only: aice, vice, vsno

      type(lagr_point), intent(in) :: pt

      real(kind=dbl_kind), intent(out) :: cai, chi, chs, ctsfc

!
      integer (kind=int_kind) :: bidxi, bidxj

!
      bidxi = ceiling(pt%bpos(1))
      bidxj = ceiling(pt%bpos(2))

      cai   = aice(bidxi,bidxj,pt%bid)
      chi   = vice(bidxi,bidxj,pt%bid)
      chs   = vsno(bidxi,bidxj,pt%bid)

      ! TODO: Tsfc & more parameters
      ctsfc = 0.0

      end subroutine get_lagr_ice_status


!========================================!===========
!
! Utilities - Logical-Physical Mappings
!
!======================================== DONE
!
! Mapping logical location to geolocation
!    Operates on the structure of Lagrangian point
!
!======================================== 
    
    subroutine logical_to_physical( pt )

    use ice_grid, only: ULON, ULAT
    
    type (lagr_point), intent(inout) :: pt

!
    integer(int_kind) :: bposi, bposj, bidxi, bidxj
    real(dbl_kind), dimension(2,0:1,0:1) :: ulatlon
    real(dbl_kind), dimension(2) :: lpos, latlon
    
!
    bidxi = ceiling(pt%bpos(1))
    bidxj = ceiling(pt%bpos(2))

    ulatlon(1,0,0) = ULON(bidxi-1,bidxj-1,pt%bid)
    ulatlon(1,1,0) = ULON(bidxi  ,bidxj-1,pt%bid) 
    ulatlon(1,1,1) = ULON(bidxi  ,bidxj  ,pt%bid)
    ulatlon(1,0,1) = ULON(bidxi-1,bidxj  ,pt%bid)

    ulatlon(2,0,0) = ULAT(bidxi-1,bidxj-1,pt%bid)
    ulatlon(2,1,0) = ULAT(bidxi  ,bidxj-1,pt%bid) 
    ulatlon(2,1,1) = ULAT(bidxi  ,bidxj  ,pt%bid)
    ulatlon(2,0,1) = ULAT(bidxi-1,bidxj  ,pt%bid)
    
    call logical_to_physical_core( ulatlon, pt%lpos, latlon )
    
    pt%latlon(:) = latlon(:)
    
    end subroutine logical_to_physical
    
 
!======================================== DONE
!
! Core subroutine for mapping logical location to geolocation
!    LATLON should be LONG-LAT in terms of layout
!
!
!      IV          III
!     (0,1)-------(1,1)
!       |           |
!       |           |
!       |           |
!       |           |
!       |           |
!     (0,0)-------(1,0)
!       I          II
!
!
!======================================== 

    subroutine logical_to_physical_core( ulatlon, lpos, latlon )
    
    real(dbl_kind), intent(in) :: ulatlon(2,0:1,0:1), lpos(2)

    real(dbl_kind), intent(out) :: latlon(2)

!    
    real(dbl_kind) :: vtcs(3,4), dummy(3), p1(3), p2(3), p(3)
    integer(int_kind) :: i
    
!   
    vtcs(:,1) = latlon2pos(ulatlon(:,0,0))
    vtcs(:,2) = latlon2pos(ulatlon(:,1,0))
    vtcs(:,3) = latlon2pos(ulatlon(:,1,1))
    vtcs(:,4) = latlon2pos(ulatlon(:,0,1))

    if( (lpos(1) < 0.0) .or.   &
        (lpos(1) > 1.0) .or.   &
        (lpos(2) < 0.0) .or.   &
        (lpos(2) > 1.0) ) then
       write(*,*) 'Unbounded Lagrangian point!'

       latlon(:) = 0.0
       return 
    end if  

    dummy(:) = bignum
    
    if     (lpos(1) <= 0.5 .and. lpos(2) >= lpos(1) .and. lpos(1) + lpos(2) <= 1.0) then
        p1 = interpolate( vtcs(:,1), vtcs(:,2),     dummy, vtcs(:,4), 3, lpos )
        p2 = interpolate( vtcs(:,1),     dummy, vtcs(:,3), vtcs(:,4), 2, lpos )
    elseif (lpos(2) <= 0.5 .and. lpos(1) >  lpos(2) .and. lpos(1) + lpos(2) <= 1.0) then
        p1 = interpolate( vtcs(:,1), vtcs(:,2),     dummy, vtcs(:,4), 3, lpos )
        p2 = interpolate( vtcs(:,1), vtcs(:,2), vtcs(:,3),     dummy, 4, lpos )
    elseif (lpos(1) >  0.5 .and. lpos(1) >  lpos(2) .and. lpos(1) + lpos(2) >  1.0) then
        p1 = interpolate( vtcs(:,1), vtcs(:,2), vtcs(:,3),     dummy, 4, lpos )
        p2 = interpolate(     dummy, vtcs(:,2), vtcs(:,3), vtcs(:,4), 1, lpos )
    elseif (lpos(2) >  0.5 .and. lpos(2) >= lpos(1) .and. lpos(1) + lpos(2) >  1.0) then
        p1 = interpolate( vtcs(:,1),     dummy, vtcs(:,3), vtcs(:,4), 2, lpos )
        p2 = interpolate(     dummy, vtcs(:,2), vtcs(:,3), vtcs(:,4), 1, lpos )
    end if
    
    do i = 1, 3
        p(i) = (p1(i)+p2(i))/2
    end do
    
    latlon = pos2latlon(p)
    
    end subroutine logical_to_physical_core

    
!======================================== DONE
!
! Core function for triangular interpolation
!
!======================================== 

    function interpolate( p1, p2, p3, p4, miss, uv ) result(p)
    
    real(dbl_kind), dimension(3), intent(inout) :: p1, p2, p3, p4

    integer(int_kind), intent(in) :: miss

    real(dbl_kind), dimension(2), intent(in) :: uv
    
    real(dbl_kind), dimension(3) :: p

!
    integer(int_kind) :: i
    real(dbl_kind), dimension(3) :: t1, t2
  
!
    if (miss == 1) then
        do i = 1, 3
            p1(i) = p2(i) + p4(i) - p3(i)
        end do
    elseif (miss == 2) then
        do i = 1, 3
            p2(i) = p1(i) + p3(i) - p4(i)
        end do
    elseif (miss == 3) then
        do i = 1, 3
            p3(i) = p2(i) + p4(i) - p1(i)
        end do
    elseif (miss == 4) then
        do i = 1, 3
            p4(i) = p1(i) + p3(i) - p2(i)
        end do
    end if
        
    do i = 1, 3
        t1(i) = p1(i)*(1-uv(1))+p2(i)*uv(1)
        t2(i) = p4(i)*(1-uv(1))+p3(i)*uv(1)
    end do
    
    do i = 1, 3
        p(i) = t1(i)*(1-uv(2))+t2(i)*uv(2)
    end do
    
    end function interpolate
    
    
!======================================== DONE
!
! Core function for transforming LatLon to 3D location 
!
!======================================== 
    
    function latlon2pos( latlon ) result(pos3)
    
    real(dbl_kind), dimension(2), intent(in) :: latlon
    real(dbl_kind), dimension(3) :: pos3
    
!
   !pos3(1) = dcos(latlon(1))*dcos(latlon(2))
   !pos3(2) = dcos(latlon(1))*dsin(latlon(2))
   !pos3(3) = dsin(latlon(1))
    pos3(1) = dcos(latlon(2))*dcos(latlon(1))
    pos3(2) = dcos(latlon(2))*dsin(latlon(1))
    pos3(3) = dsin(latlon(2))
    
    end function latlon2pos
    
    
!======================================== DONE
!
! Core function for transforming 3D location to LatLon
!
!======================================== 
    
    function pos2latlon( pos3 ) result(latlon)
    
    use ice_constants, only: pi
    real(dbl_kind), dimension(3), intent(in) :: pos3
    real(dbl_kind), dimension(2) :: latlon
    
!
    real(dbl_kind) :: r, xr
    
!
    r = dsqrt( pos3(1)*pos3(1) + pos3(2)*pos3(2) )

    if ( r > 0.0d0 ) then
        xr = dacos( pos3(1) / r )
        if (pos3(2) >= 0.0) then
            latlon(1) =  xr
        else 
           !latlon(1) =  2*pi - xr
            latlon(1) = -xr
        end if
    
        r = dsqrt( pos3(1)*pos3(1) + pos3(2)*pos3(2) + pos3(3)*pos3(3) )
        latlon(2) = dasin( pos3(3)/r )
    else
        if ( pos3(3) >= 0.0d0 ) then
            latlon(2) =  pi/2.0d0
        else
            latlon(2) = -pi/2.0d0
        end if

        latlon(1) = 0.0d0
    end if
    
    end function pos2latlon
    

!====================================================
!
! Utilities - State Transitions
!
!======================================== 
!
! Activate a Lagrangian point with location within the block (bpos)
!     Type info and block info (i.e., index) is provided
!
!========================================

      subroutine lagr_activate(pt, ctime, ptype, cblock, pos_info, latlon_info)

      use ice_blocks, only: block, get_block, nghost
      use ice_domain, only: blocks_ice
      use ice_constants, only: secday
      use ice_calendar, only: time

      type (lagr_point), intent(inout) :: pt

      real(kind=dbl_kind), intent(in) :: ctime

      integer(kind=int_kind), intent(in) :: ptype

      type (block), intent(in) :: cblock

      real(kind=dbl_kind), dimension(2), intent(in) :: pos_info

      real(kind=dbl_kind), dimension(2), intent(in), optional :: latlon_info

!
      integer (kind=int_kind) :: ptype_local, ilo_glob, jlo_glob
      real(kind=dbl_kind) :: ctime_local, life_local
      real(kind=dbl_kind), dimension(2) :: bpos, gpos

!
      if (is_lagr_point_active(pt)) then
         write(nu_diag_l,*) 'ERROR! Activiating an active point'
         return
      end if

    ! Setting status
      pt%stat = LAGR_ACTIVE

    ! Decoding PTYPE
      if ((ptype == LAGR_MIGRATING_BUOY) .or. (ptype == LAGR_MIGRATING_BUOY_ERR)) then
         if (ctime >= 0) then
            ptype_local = LAGR_PHYSICAL_BUOY
         else
            ptype_local = LAGR_VIRTUAL_BUOY
         end if
      else
         ptype_local = ptype
      end if

    ! Decoding CTIME
      if ((ptype == LAGR_MIGRATING_BUOY) .or. (ptype == LAGR_MIGRATING_BUOY_ERR)) then
         if (ctime >=0) then
            ctime_local =  ctime
         else
            ctime_local = -ctime -secday
         end if

         life_local = nint((time - ctime_local)/lagr_dt)
      else
        !Special treatment to newly deployed points
         ctime_local = ctime

         life_local = nint((time - lagr_dt - ctime_local)/lagr_dt)
      end if

    ! Decoding POS_INFO
      if (ptype == LAGR_VIRTUAL_BUOY) then
         bpos(:) = pos_info(:)
      elseif (ptype == LAGR_PHYSICAL_BUOY) then
         ! TODO: no support for physical buoys yet
      elseif ((ptype == LAGR_MIGRATING_BUOY) .or. (ptype == LAGR_MIGRATING_BUOY_ERR)) then
         gpos(:) = pos_info(:)

         ilo_glob = cblock%i_glob(cblock%ilo)
        !ihi_glob = cblock%i_glob(cblock%ihi)
         jlo_glob = cblock%j_glob(cblock%jlo)
        !jhi_glob = cblock%j_glob(cblock%jhi)

         bpos(1) = gpos(1) - (ilo_glob-1) + nghost
         bpos(2) = gpos(2) - (jlo_glob-1) + nghost
      end if

     !Sanity check
      if ((bpos(1) <= cblock%ilo-1) .or.  &
          (bpos(1) >  cblock%ihi  ) .or.  &
          (bpos(2) <= cblock%jlo-1) .or.  &
          (bpos(2) >  cblock%jhi  )) then
         write(nu_diag_l,*) 'ERROR! Activating a Lagrangian Point NOT in Block!'
         pt%stat = LAGR_INACTIVE   ! resetting from CLAIMED
         return
      end if

      pt%bid = cblock%local_id

      pt%bpos(:) = bpos(:)
      call lagr_update_local_grid_vars(pt)

      pt%lpos(1) = pt%bpos(1) - floor(pt%bpos(1))
      pt%lpos(2) = pt%bpos(2) - floor(pt%bpos(2))

      call lagr_compute_gpos(pt)
      call lagr_compute_latlon(pt)

      if (present(latlon_info)) then
         pt%lon_init = latlon_info(1)
         pt%lat_init = latlon_info(2)
      else
         pt%lon_init = pt%latlon(1)
         pt%lat_init = pt%latlon(2)
      end if

      pt%time_init = ctime_local
      pt%type      = ptype_local
      pt%life      = life_local

      if (pt%type == LAGR_VIRTUAL_BUOY) then
         pt%max_life = LAGR_VIRTUAL_BUOY_MAX_LIFE
      end if


     !Optional output for logging?
      if (ptype == LAGR_MIGRATING_BUOY) then
         call lagr_info_dump(pt,LAGR_MIGRATION_STR,     nu_diag_l)
      elseif (ptype == LAGR_MIGRATING_BUOY_ERR) then
         call lagr_info_dump(pt,LAGR_MIGRATION_STR_ERR, nu_diag_l)
         pt%stat = LAGR_INACTIVE
      else
         call lagr_info_dump(pt,LAGR_ACTIVATION_STR,nu_diag_l)
      endif

      end subroutine lagr_activate


!======================================== DONE
! 
! Activation of virtual buoys
!
!   Buoys are distributed evenly within the cell
!   This subroutine is called every ACTIVATION_FREQ steps
!
!========================================

      subroutine lagr_deploy_virtual_buoys (cblock,ctime)

      use ice_domain, only: nblocks
      use ice_blocks, only: block
      use ice_grid, only: tmask, umask, tlat, tlon
      use ice_state, only: aice

      type(block), intent(in) :: cblock

      real(kind=dbl_kind), intent(in) :: ctime

!
      integer (kind=int_kind) :: i, j, li, lj, iblk, lidx, cnt, ig, jg

      real(kind=dbl_kind), dimension(2) :: cpos

      real(kind=dbl_kind), parameter :: STEP  = 1.0d0/(LAGR_ACTIVATION_DNSTY+1)
      integer (kind=int_kind), parameter :: STEPL = nint(1.0d0/LAGR_ACTIVATION_DNSTY)

      real(kind=dbl_kind), parameter :: LAT_LIMIT_DEPLOY = 0.0d0

      logical (kind=log_kind), parameter :: DEPLOY_DBG  = .true.
      logical (kind=log_kind), parameter :: DEPLOY_DBG2 = .false.

!
      iblk = cblock%local_id     

      if (lagr_debug .and. DEPLOY_DBG2) then
         write(nu_diag_l,*) 'Activating virtual buoys on block ', iblk, nblocks
      end if


      cnt = 0

      do j = cblock%jlo, cblock%jhi
      do i = cblock%ilo, cblock%ihi

         if (LAGR_ACTIVATION_DNSTY >= 1.0) then
            if (.not. tmask(i,j,iblk)) then
               cycle     ! No buoy on land
            elseif ((.not. umask(i  ,j  ,iblk)) .or.  &
                    (.not. umask(i-1,j  ,iblk)) .or.  &
                    (.not. umask(i  ,j-1,iblk)) .or.  &
                    (.not. umask(i-1,j-1,iblk)) ) then
               cycle     ! No buoy on stagnant locations  ??? 
            elseif (aice(i,j,iblk) < lagr_aice_thres) then
               cycle     ! No buoy in water
            end if

            if (DEPLOY_DBG) then
               ! Killing deployment in south to 85N
               if (tlat(i,j,iblk) < LAT_LIMIT_DEPLOY/rad_to_deg) then
                  cycle
               end if
            end if

           !Actual activations
            cpos(2) = j - 1 + STEP
            do lj = 1, nint(LAGR_ACTIVATION_DNSTY)

               cpos(1) = i - 1 + STEP
               do li = 1, nint(LAGR_ACTIVATION_DNSTY)

                  call lagr_claim_slot(lidx)
                  if (lidx == 0) then
                     cycle  ! Can also return
                  end if

                  call lagr_activate( lagr_points(lidx),          &
                                      ctime,  LAGR_VIRTUAL_BUOY,  &
                                      cblock, cpos )

                  cnt = cnt + 1

                  cpos(1) = cpos(1) + STEP
               end do  !li

               cpos(2) = cpos(2) + STEP
            end do  !lj

         else
           !Latitude control
            if (tlat(i,j,iblk) < LAT_LIMIT_DEPLOY/rad_to_deg) then
               cycle
            end if

           !Location control
            ig = cblock%i_glob(i)
            jg = cblock%j_glob(j)

            if ((mod(ig+STEPL/2,STEPL) .ne. 0) .or.  &
                (mod(jg+STEPL/2,STEPL) .ne. 0)) then
               cycle
            end if

           !Mask & AICE control
            if (.not. tmask(i,j,iblk)) then
               cycle     ! No buoy on land
            elseif ((.not. umask(i  ,j  ,iblk)) .or.  &
                    (.not. umask(i-1,j  ,iblk)) .or.  &
                    (.not. umask(i  ,j-1,iblk)) .or.  &
                    (.not. umask(i-1,j-1,iblk)) ) then
               cycle     ! No buoy on stagnant locations  ??? 
            elseif (aice(i,j,iblk) < lagr_aice_thres) then
               cycle     ! No buoy in water
            end if
           
           !Actual activations
            cpos(1) = i - 0.5d0
            cpos(2) = j - 0.5d0

            call lagr_claim_slot(lidx)
            if (lidx == 0) then
               cycle  ! Can also return
            end if

            call lagr_activate( lagr_points(lidx),          &
                                ctime,  LAGR_VIRTUAL_BUOY,  &
                                cblock, cpos )

            cnt = cnt + 1

         end if 
      end do
      end do

      if (lagr_debug .and. cnt > 0) then
         write(nu_diag_l,*) 'Activated ', cnt, ' virtual buoys for block:', iblk
      end if
     
      end subroutine lagr_deploy_virtual_buoys


!======================================== 

! Deactivate a Lagrangian point with location within the block (bpos)
!     Type info and block info (i.e., index) is provided
!
!========================================

      subroutine lagr_deactivate(pt, reason)

      type(lagr_point), intent(inout) :: pt

      integer (kind=int_kind), intent(in) :: reason

!
      pt%stat = reason

      if (reason == LAGR_DEACTIVATED) then
         call lagr_info_dump(pt,LAGR_DEACTIVATION_STR,         nu_diag_l)
      elseif (reason == LAGR_MELTED)  then
         call lagr_info_dump(pt,LAGR_DEACTIVATION_MELT_STR,    nu_diag_l)
      elseif (reason == LAGR_ONLAND)  then
         call lagr_info_dump(pt,LAGR_DEACTIVATION_LAND_STR,    nu_diag_l)
      elseif (reason == LAGR_RIDGED)  then
         call lagr_info_dump(pt,LAGR_DEACTIVATION_RIDG_STR,    nu_diag_l)
      elseif (reason == LAGR_MIGRATED)  then
         call lagr_info_dump(pt,LAGR_DEACTIVATION_MIGR_STR,    nu_diag_l)
      elseif (reason == LAGR_MIGRATED_ERR)  then
         call lagr_info_dump(pt,LAGR_DEACTIVATION_MIGR_STR_ERR,nu_diag_l)
      elseif (reason < 0)  then
         call lagr_info_dump(pt,LAGR_DEACTIVATION_UNKN_STR,    nu_diag_l)
      else
         call lagr_info_dump(pt,LAGR_ABNORMAL_STR,             nu_diag_l)
         pt%stat = LAGR_XXX
      endif
      
      end subroutine lagr_deactivate


!======================================== DONE
!
! Claiming a Lagrangian point slot
!    Returning index of the slot, setting its status to CLAIMED
!    Returning 0 on failure (no slot available)
!
!========================================

      subroutine lagr_claim_slot(idx)

      integer (kind=int_kind), intent(out) :: idx

!
      integer (kind=int_kind) :: i

!
      if (lagr_active_count == lagr_buffer_size) then
         idx = 0
         return
      end if

      do i = 1, lagr_buffer_size
         if (.not. lagr_slot_active(i)) then
            idx = i  

            lagr_slot_active(i) = .true.
            lagr_active_count = lagr_active_count + 1

            lagr_points(idx)%stat = LAGR_CLAIMED

            return
         end if
      end do

      end subroutine lagr_claim_slot


!======================================== DONE
!
! Increase life of a Lagrangian point
!
!========================================

      subroutine lagr_increase_life( pt )

      type (lagr_point), intent(inout) :: pt

!
      pt%life = pt%life + 1

      end subroutine lagr_increase_life


!======================================== DONE
!
! Dump Lagr. information to string
!
!========================================
      
      subroutine lagr_info_str( pt, ostr )
      
      use ice_grid, only: tlat, tlon

      type (lagr_point), intent(in) :: pt
      
      character (len=char_len_long), intent(out) :: ostr
      
!
      integer (kind=int_kind) :: i,j

      character (len=char_len_long), parameter ::  &
         ostr_fmt  = '(a6,a,a3, F11.5,a,F11.5,a2,E15.8,a2, i9, a4,F11.5,a,F11.5,a2)',  &
         ostr_fmt2 = '(a6,a,a3, F11.5,a,F11.5,a2,E15.8,a2, i9, a4,F11.5,a,F11.5,a5, a4,F9.3,a,F9.3,a2, a4,F7.3,a,F7.3,a2, a4,F6.3,a,F6.3,a2, a5,F8.3,a,F8.3,a2, a4,F8.3,a,F8.3,a)'
                     ! Header   lon/lat/time_init              lon/lat              gpos               bpos               lpos               T-lon/lat          HTN/E(1,1)
         
      character (len=1) :: type_str
      
!
      if (pt%type == LAGR_PHYSICAL_BUOY) then
         type_str = 'P'
      elseif (pt%type == LAGR_VIRTUAL_BUOY) then
         type_str = 'V'
      else
         type_str = 'U'
      end if
      
      if (lagr_debug) then
         i = ceiling(pt%bpos(1))
         j = ceiling(pt%bpos(2))

         write(ostr,ostr_fmt2) &
            'LAGR-<', type_str, '>:[',                        &
            pt%lon_init*rad_to_deg, ',', pt%lat_init*rad_to_deg, ',', pt%time_init, '] ',  &
            pt%life,   &
            ' GL[', pt%latlon(1)*rad_to_deg, ',', pt%latlon(2)*rad_to_deg, ' ] ~~',  &
            ' GP{', pt%gpos(1),              ',', pt%gpos(2),              ' }',     &
            ' BP(', pt%bpos(1),              ',', pt%bpos(2),              ' )',     &
            ' LP<', pt%lpos(1),              ',', pt%lpos(2),              ' >',     &
           ' TGL[', tlon(i,j,pt%bid)*rad_to_deg,                                    &
                                             ',', tlat(i,j,pt%bid)*rad_to_deg,      &
                                                                           ' ]',     &
            ' SC`', pt%dxt(1,1),             ',', pt%dyt(1,1),             '`'
      else
         write(ostr,ostr_fmt ) 'LAGR-<', type_str, '>:[',  &
            pt%lon_init, ',', pt%lat_init, ',', pt%time_init, '] ',  &
            pt%life,   &
            ' GL[', pt%latlon(1), ',', pt%latlon(2), ']'
      end if
      
      end subroutine lagr_info_str


!======================================== DONE
!
! Dump Lagr. information 
!
!========================================

      subroutine lagr_info_dump( pt, prefix, hdl )

      type(lagr_point), intent(in) :: pt

      character (len=20), intent(in) :: prefix

      integer (kind=int_kind), intent(in) :: hdl

!
      character (len=char_len_long) :: ostr
      character (len=char_len+char_len_long) :: lstr

!
      call lagr_info_str(pt,ostr)

      write(hdl,*) prefix, ostr

      end subroutine lagr_info_dump
 
      
!========================================
!
! Utilities - Filling & Updating
!
!======================================== 
!
! Subroutine to compute global position with location within block (i.e., bpos)
!    NOTE: periodic boundary & tripolar boundary
!
!========================================

      subroutine lagr_compute_gpos( pt )

      use ice_blocks, only: block, get_block, nghost, nx_block, ny_block
      use ice_domain_size, only: nx_global, ny_global
      use ice_domain, only: blocks_ice
      use ice_grid, only: grid_type

      type (lagr_point), intent(inout) :: pt

!
      integer (kind=int_kind) :: i, j, i_glob, j_glob, ilo_glob, jlo_glob
      type(block) :: this_block
      real (kind=dbl_kind), dimension(2) :: gpos

!
      this_block = get_block(blocks_ice(pt%bid),pt%bid)
      
      ilo_glob = this_block%i_glob(this_block%ilo)
      jlo_glob = this_block%j_glob(this_block%jlo)

      gpos(1) = pt%bpos(1) - nghost + (ilo_glob-1) 
      gpos(2) = pt%bpos(2) - nghost + (jlo_glob-1) 

      pt%gpos(:) = gpos(:)


      !East-West
      if ((trim(grid_type) == 'displaced_pole') .or.  &
          (trim(grid_type) == 'tripole')) then
         if (pt%gpos(1) <= 0.0d0) then
            pt%gpos(1) = pt%gpos(1) + nx_global
         elseif (pt%gpos(1) > nx_global) then
            pt%gpos(1) = pt%gpos(1) - nx_global
         end if
      else
        !No need to be grumpy, the Lagr. point will die
        !write(nu_diag_l,*) 'Warning! Regional Model!'
      end if
      
      !South

      !North
      if (trim(grid_type) == 'tripole') then
         if (pt%gpos(2) > ny_global) then
            gpos(1) =   nx_global - gpos(1)
            gpos(2) = 2*ny_global - gpos(2)

            pt%gpos(:) = gpos(:)
         end if
      end if

      end subroutine lagr_compute_gpos


!======================================== DONE
!
! Subroutine called when geolocation needs to be computed
!    Position within block is required
!    A simple wrapper to call to logical_to_physical
!
!========================================

      subroutine lagr_compute_latlon (pt)

      type (lagr_point), intent(inout) :: pt

!
      call logical_to_physical(pt)

      end subroutine lagr_compute_latlon


!======================================== DONE

      subroutine lagr_update_local_grid_vars( pt )

      use ice_grid, only: htn, hte, dxt, dyt, dxu, dyu

      type (lagr_point), intent(inout) :: pt

!
      integer (kind=int_kind) :: il,ih,jl,jh, iblk

      logical (kind=log_kind), parameter :: DBG_LOCAL = .false.

!
      if (.not. is_lagr_point_active(pt)) then
          return
      end if

      il = floor(pt%bpos(1))
      jl = floor(pt%bpos(2))
      ih = ceiling(pt%bpos(1))
      jh = ceiling(pt%bpos(2))

      if (il == ih) then
         il = ih - 1
      end if
      if (jl == jh) then
         jl = jh - 1
      end if


      iblk = pt%bid

      pt%htn(0:2,0:2) = htn(il:ih+1,jl:jh+1,iblk)
      pt%hte(0:2,0:2) = hte(il:ih+1,jl:jh+1,iblk)

      pt%dxt(0:2,0:2) = dxt(il:ih+1,jl:jh+1,iblk)
      pt%dyt(0:2,0:2) = dyt(il:ih+1,jl:jh+1,iblk)

      pt%dxu(0:2,0:2) = dxu(il:ih+1,jl:jh+1,iblk)
      pt%dyu(0:2,0:2) = dyu(il:ih+1,jl:jh+1,iblk)

      if (lagr_debug .and. DBG_LOCAL) then
         write(nu_diag_l,*) 'DXT   ', il, ih, jl, jh, iblk, dxt(ih,jh,iblk)
         write(nu_diag_l,*) 'DXT@PT', il, ih, jl, jh, iblk, dxt(ih,jh,iblk)
      end if

      end subroutine lagr_update_local_grid_vars


!======================================== DONE
!
! Updating offset of the downstream location of corner (U) points
!    This info is available through transport remapping
!
!======================================== 

      subroutine lagr_update_downstream(pt)

      use ice_transport_remap, only: dsx, dsy

      type (lagr_point), intent(inout) :: pt

!
      integer (kind=int_kind) :: il,ih,jl,jh, iblk

      logical (kind=log_kind), parameter :: FORCE_IMMOBILE_DBG = .false.
      logical (kind=log_kind), parameter :: FORCE_SHIFT_DBG = .false.

!
      if (.not. is_lagr_point_active(pt)) then
          return
      end if

      il = floor(pt%bpos(1))
      jl = floor(pt%bpos(2))
      ih = ceiling(pt%bpos(1))
      jh = ceiling(pt%bpos(2))

      iblk = pt%bid

      pt%dsx(:,:) = dsx(il:ih,jl:jh,iblk)
      pt%dsy(:,:) = dsy(il:ih,jl:jh,iblk)


      if (FORCE_IMMOBILE_DBG) then
         pt%dsx(:,:) = 0.0d0
         pt%dsy(:,:) = 0.0d0
      end if
      if (FORCE_SHIFT_DBG) then
         pt%dsx(:,:) =  2000.0d0
         pt%dsy(:,:) =  2000.0d0
      end if

      end subroutine lagr_update_downstream


!========================================
!
! Utilities - Lagrangian Tracking
!
!======================================== DONE
!
! Core subroutine to carry out Lagrangian tracking for a specific point
!
!======================================== 

      subroutine lagr_do_tracking( pt )

      type(lagr_point), intent(inout) :: pt

!
      real(kind=dbl_kind), dimension(2) :: new_lpos
      logical (kind=log_kind) :: t_changed

!
    ! Get corner displacements
      call lagr_update_downstream(pt)

      if (lagr_scheme == LAGR_REMAPPING) then
         call lagr_core_bmap( pt%dxt, pt%dyt, pt%dxu, pt%dyu,  &
                              pt%dsx, pt%dsy, pt%lpos,         &
                              new_lpos )

      elseif (lagr_scheme == LAGR_TRACMASS) then
         ! TODO: no support for TRACMASS algorithm yet
      end if

      if ((new_lpos(1) > 1.0d0) .or. (new_lpos(1) < 0.0d0) .or.  &
          (new_lpos(2) > 1.0d0) .or. (new_lpos(1) < 0.0d0)) then
         t_changed = .true.
      else
         t_changed = .false.
      end if

    ! Backup & update position within block
      pt%bpos_old(:) = pt%bpos(:)
      pt%bpos(1) = pt%bpos(1) + (new_lpos(1) - pt%lpos(1))
      pt%bpos(2) = pt%bpos(2) + (new_lpos(2) - pt%lpos(2))

      ! local position
      pt%lpos(1) = pt%bpos(1) - floor(pt%bpos(1))
      pt%lpos(2) = pt%bpos(2) - floor(pt%bpos(2))

      ! global position
      call lagr_compute_gpos(pt)

      if (.not. is_lagr_point_still_in_block(pt)) then
         pt%migrate = .true.
      else
         pt%migrate = .false.

         if (t_changed) then
            call lagr_update_local_grid_vars(pt)
         end if
         call lagr_compute_latlon(pt)
      end if

      end subroutine lagr_do_tracking


!======================================== DONE
!
! Core subroutine for Lagrangian tracking
!    Stateless, no structure info is passed 
!
!======================================== 

      subroutine lagr_core_bmap (dxt, dyt,  &
                                 dxu, dyu,  &
                                 dx,  dy,   &
                                 old_pos,   &
                                 new_pos)
!
      real (kind=dbl_kind), dimension(0:2,0:2), intent(in) :: dxt, dyt, dxu, dyu

      real (kind=dbl_kind), dimension(0:1,0:1), intent(in) :: dx, dy

      real (kind=dbl_kind), dimension(2), intent(in) :: old_pos

      real (kind=dbl_kind), dimension(2), intent(out) :: new_pos

!
      real (kind=dbl_kind), dimension(0:1,0:1) :: ndx, ndy
      real (kind=dbl_kind), dimension(4) :: alpha, beta
      integer (kind=int_kind) :: i,j

      logical (kind=log_kind), parameter :: DBG_LOCAL = .false.

!
      if (lagr_debug .and. DBG_LOCAL) then
         write(nu_diag_l,*) 'DXT:', dxt(0,0), dxt(1,0), dxt(1,1), dxt(1,0)
         write(nu_diag_l,*) 'DYT:', dyt(0,0), dyt(1,0), dyt(1,1), dyt(1,0)
         write(nu_diag_l,*) 'DXU:', dxu(0,0), dxu(1,0), dxu(1,1), dxu(1,0)
         write(nu_diag_l,*) 'DYU:', dyu(0,0), dyu(1,0), dyu(1,1), dyu(1,0)
         write(nu_diag_l,*) 'DX :', dx (0,0), dx (1,0), dx (1,1), dx (1,0)
         write(nu_diag_l,*) 'DY :', dy (0,0), dy (1,0), dy (1,1), dy (1,0)
         write(nu_diag_l,*) 'POS:', old_pos(1), old_pos(2)
      end if

      do j = 0, 1
      do i = 0, 1
          ndx(i,j) = i + dx(i,j) / dxu(i,j)
          ndy(i,j) = j + dy(i,j) / dyu(i,j)
      end do
      end do

      alpha(1) = ndx(0,0)
      alpha(2) = ndx(1,0) - ndx(0,0)
      alpha(3) = ndx(0,1) - ndx(0,0)
      alpha(4) = ndx(0,0) + ndx(1,1) - ndx(0,1) - ndx(1,0)

      beta(1)  = ndy(0,0)
      beta(2)  = ndy(1,0) - ndy(0,0)
      beta(3)  = ndy(0,1) - ndy(0,0)
      beta(4)  = ndy(0,0) + ndy(1,1) - ndy(0,1) - ndy(1,0)
 
      new_pos(1) = alpha(1) + alpha(2)* old_pos(1)  &
                 + alpha(3) * old_pos(2) + alpha(4) * old_pos(1) * old_pos(2)
      new_pos(2) = beta(1)  + beta(2) * old_pos(1)  &
                 + beta(3)  * old_pos(2) + beta(4)  * old_pos(1) * old_pos(2)

      end subroutine lagr_core_bmap


!***********************************************************************   

!     subroutine lagr_util_fill(pt)

!     use ice_blocks, only: block, get_block

!     type (lagr_point), intent(inout) :: pt

!
!     type(block) :: this_block

!     integer (kind=int_kind), dimension(2) :: bidx

!     integer (kind=int_kind) :: iblk
!
!     if (pt%status .lt. 1) then
!         return
!     end if

!     iblk = pt%bid

!     this_block = get_block(blocks_ice(iblk),iblk)

!     bidx(1) = ceiling(pt%bpos(1))
!     bidx(2) = ceiling(pt%bpos(2))

!     pt%uvel(0:2) = uvel(bidx(1)-1:bidx(1)+1,bidx(2)-1:bidx(2)+1,iblk)
!     pt%vvel(0:2) = vvel(bidx(1)-1:bidx(1)+1,bidx(2)-1:bidx(2)+1,iblk)

!     pt%htn(0:2)  =  htn(bidx(1)-1:bidx(1)+1,bidx(2)-1:bidx(2)+1,iblk)
!     pt%hte(0:2)  =  hte(bidx(1)-1:bidx(1)+1,bidx(2)-1:bidx(2)+1,iblk)

!     end subroutine lagr_util_fill



      end module ice_transport_driver


!=======================================================================
