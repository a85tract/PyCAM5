module restart_physics

  use shr_kind_mod,       only: r8 => shr_kind_r8
  use iso_c_binding,      only: c_double, c_int64_t, c_loc, c_ptr
  use spmd_utils,         only: masterproc
  use ppgrid,             only: pcols, pver, pverp, begchunk, endchunk
  use constituents,       only: pcnst
  use co2_cycle,          only: c_i, co2_transport
  use dyn_grid,           only: ptimelevels
  use radae,              only: abstot_3d, absnxt_3d, emstot_3d, initialize_radbuffer, ntoplw
  use comsrf,             only: sgh, sgh30, landm, trefmxav, trefmnav, & 
       fsnt, flns, fsns, fsds, flnt, initialize_comsrf
  use ioFileMod
  use cam_abortutils,     only: endrun
  use units,              only: getunit
  use camsrfexch,         only: cam_in_t, cam_out_t
  use cam_control_mod,    only: adiabatic, ideal_phys
  use cam_logfile,        only: iulog
  use pio,                only: file_desc_t, io_desc_t, var_desc_t, &
                                pio_double, pio_int, pio_noerr, &
                                pio_seterrorhandling, pio_internal_error, pio_bcast_error, &
                                pio_inq_dimid, pio_inq_varname, pio_inq_varid, &
                                pio_def_var, pio_def_dim, &
                                pio_put_att, pio_put_var, pio_get_var, &
                                pio_write_darray, pio_read_darray
  use cospsimulator_intr, only: docosp
  use radiation,          only: cosp_cnt_init, cosp_cnt

  implicit none
  private
  save
!
! Public interfaces
!
  public :: write_restart_physics    ! Write the physics restart info out
  public :: read_restart_physics     ! Read the physics restart info in
  public :: get_abs_restart_filepath ! Get the name of the restart filepath
  public :: init_restart_physics

!
! Private data
!
    integer :: nrg2 = -1         ! Abs/ems restart dataset unit number
    character(len=256) :: pname  ! Full abs-ems restart filepath
    logical :: restart_write_entered_logged = .false.

    type(var_desc_t) :: trefmxav_desc, trefmnav_desc, flwds_desc, landm_desc, sgh_desc, &
         sgh30_desc, solld_desc, co2prog_desc, co2diag_desc, sols_desc, soll_desc, &
         solsd_desc, fsnt_desc, flns_desc, emstot_desc, absnxt_desc(4), &
         pblh_desc,  tpert_desc, qpert_desc, flnt_desc, fsds_desc, fsns_desc

    type(var_desc_t) :: bcphidry_desc, bcphodry_desc, ocphidry_desc, ocphodry_desc, &
       dstdry1_desc, dstdry2_desc, dstdry3_desc, dstdry4_desc

    type(var_desc_t) :: cflx_desc(pcnst)

    type(var_desc_t), allocatable :: abstot_desc(:)

    type(var_desc_t) :: cospcnt_desc

    logical :: use_native_restart_pack_impl = .false.
    logical :: restart_pack_impl_selected = .false.
    logical :: restart_pack_entered_logged = .false.
    logical :: restart_init_entered_logged = .false.

    interface
       function init_restart_physics_codon(stage_c) result(stage_out) bind(c, name="init_restart_physics_codon")
         import c_int64_t
         integer(c_int64_t), value :: stage_c
         integer(c_int64_t) :: stage_out
       end function init_restart_physics_codon
       function write_restart_physics_codon(stage_c) result(stage_out) bind(c, name="write_restart_physics_codon")
         import c_int64_t
         integer(c_int64_t), value :: stage_c
         integer(c_int64_t) :: stage_out
       end function write_restart_physics_codon
    end interface

  CONTAINS
    subroutine restart_physics_log_init_entry()
      integer(c_int64_t) :: touch_c

      if (restart_init_entered_logged) return
      touch_c = init_restart_physics_codon(1101_c_int64_t)
      if (masterproc .and. touch_c == 1101_c_int64_t) then
         write(iulog,'(A)') 'init_restart_physics direct = codon; PIO/chem/subcol restart definition native island'
      end if
      restart_init_entered_logged = .true.
    end subroutine restart_physics_log_init_entry

    subroutine restart_physics_log_write_entry()
      integer(c_int64_t) :: touch_c

      if (restart_write_entered_logged) return
      touch_c = write_restart_physics_codon(1102_c_int64_t)
      if (masterproc .and. touch_c == 1102_c_int64_t) then
         write(iulog,'(A)') 'write_restart_physics direct = codon; pack helpers direct = codon; PIO/chem/subcol restart native island'
      end if
      restart_write_entered_logged = .true.
    end subroutine restart_physics_log_write_entry

    subroutine init_restart_physics ( File, pbuf2d, hdimids)
      
    use cam_pio_utils,       only: fillvalue
    
    use physics_buffer,              only: pbuf_init_restart, physics_buffer_desc
    use dyn_grid,            only: get_horiz_grid_dim_d
    use radiation,           only: radiation_do
    use chemistry,           only: chem_init_restart
    use prescribed_ozone,    only: init_prescribed_ozone_restart
    use prescribed_ghg,      only: init_prescribed_ghg_restart
    use prescribed_aero,     only: init_prescribed_aero_restart
    use prescribed_volcaero, only: init_prescribed_volcaero_restart
    use subcol_utils,        only: is_subcol_on
    use subcol,              only: subcol_init_restart

    type(file_desc_t), intent(inout) :: file
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    integer,intent(in) :: hdimids(:)

    integer :: hdimcnt, ierr, hdim1_d, hdim2_d, i, vsize, lwrdim
    integer :: dimids(4)
    integer :: pver_id, pverp_id, pcnst_id
    integer, pointer :: ldof(:)
    character(len=4) :: num
    integer(c_int64_t) :: init_touch_c

    init_touch_c = init_restart_physics_codon(1101_c_int64_t)
    call restart_physics_log_init_entry()

    hdimcnt=size(hdimids)
    dimids(1:hdimcnt) = hdimids
    call get_horiz_grid_dim_d(hdim1_d, hdim2_d)

    call pio_seterrorhandling(File, PIO_BCAST_ERROR)
    ierr = pio_inq_dimid(File, 'lev', pver_id)
    call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)
    if(ierr/=PIO_NOERR) then
       ierr = pio_def_dim(File, 'lev', pver, pver_id)
    end if
    call pio_seterrorhandling(File, PIO_BCAST_ERROR)
    ierr = pio_inq_dimid(File, 'ilev', pverp_id)
    call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)
    if(ierr/=PIO_NOERR) then
       ierr = pio_def_dim(File, 'ilev', pverp, pverp_id)
    end if

    call pbuf_init_restart(File, pbuf2d)

    if ( .not. adiabatic .and. .not. ideal_phys )then
       
       call chem_init_restart(File)

       call init_prescribed_ozone_restart(File)
       call init_prescribed_ghg_restart(File)
       call init_prescribed_aero_restart(File)
       call init_prescribed_volcaero_restart(File)

       call pio_seterrorhandling(File, PIO_BCAST_ERROR)
       ierr = pio_inq_dimid(File, 'pcnst', pcnst_id)
       call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)
       if(ierr/=PIO_NOERR) then
          ierr = pio_def_dim(File, 'pcnst', pcnst, pcnst_id)
       end if
    
       ierr = pio_def_var(File, 'FSNT',     pio_double, hdimids, fsnt_desc)
       ierr = pio_def_var(File, 'FSNS',     pio_double, hdimids, fsns_desc)
       ierr = pio_def_var(File, 'FSDS',     pio_double, hdimids, fsds_desc)
       ierr = pio_def_var(File, 'FLNT',     pio_double, hdimids, flnt_desc)
       ierr = pio_def_var(File, 'FLNS',     pio_double, hdimids, flns_desc)
       ierr = pio_def_var(File, 'LANDM',    pio_double, hdimids, landm_desc)
       ierr = pio_def_var(File, 'SGH',      pio_double, hdimids, sgh_desc)
       ierr = pio_def_var(File, 'SGH30',    pio_double, hdimids, sgh30_desc)
       ierr = pio_def_var(File, 'TREFMXAV', pio_double, hdimids, trefmxav_desc)
       ierr = pio_def_var(File, 'TREFMNAV', pio_double, hdimids, trefmnav_desc)
       
       ierr = pio_def_var(File, 'FLWDS', pio_double, hdimids, flwds_desc)
       ierr = pio_def_var(File, 'SOLS', pio_double, hdimids, sols_desc)
       ierr = pio_def_var(File, 'SOLL', pio_double, hdimids, soll_desc)
       ierr = pio_def_var(File, 'SOLSD', pio_double, hdimids, solsd_desc)
       ierr = pio_def_var(File, 'SOLLD', pio_double, hdimids, solld_desc)

       ierr = pio_def_var(File, 'BCPHIDRY', pio_double, hdimids, bcphidry_desc)
       ierr = pio_def_var(File, 'BCPHODRY', pio_double, hdimids, bcphodry_desc)
       ierr = pio_def_var(File, 'OCPHIDRY', pio_double, hdimids, ocphidry_desc)
       ierr = pio_def_var(File, 'OCPHODRY', pio_double, hdimids, ocphodry_desc)
       ierr = pio_def_var(File, 'DSTDRY1',  pio_double, hdimids, dstdry1_desc)
       ierr = pio_def_var(File, 'DSTDRY2',  pio_double, hdimids, dstdry2_desc)
       ierr = pio_def_var(File, 'DSTDRY3',  pio_double, hdimids, dstdry3_desc)
       ierr = pio_def_var(File, 'DSTDRY4',  pio_double, hdimids, dstdry4_desc)

       if(co2_transport()) then
          ierr = pio_def_var(File, 'CO2PROG', pio_double, hdimids, co2prog_desc)
          ierr = pio_def_var(File, 'CO2DIAG', pio_double, hdimids, co2diag_desc)
       end if

       ! cam_import variables -- write the constituent surface fluxes as individual 2D arrays
       ! rather than as a single variable with a pcnst dimension.  Note that the cflx components
       ! are only needed for those constituents that are not passed to the coupler.  The restart
       ! for constituents passed through the coupler are handled by the .rs. restart file.  But
       ! we don't currently have a mechanism to know whether the constituent is handled by the
       ! coupler or not, so we write all of cflx to the CAM restart file.
       do i = 1, pcnst
          write(num,'(i4.4)') i
          ierr = pio_def_var(File, 'CFLX'//num,  pio_double, hdimids, cflx_desc(i))
       end do

    end if


    if( radiation_do('aeres')  ) then
       vsize = (pverp-ntoplw+1)
       if(vsize/=pverp) then
          ierr = pio_def_dim(File, 'lwcols', vsize, dimids(hdimcnt+1))
       else
          dimids(hdimcnt+1) = pverp_id
       end if
!
! split this into vsize variables to avoid excessive memory usage in IO
!
       allocate(abstot_desc(ntoplw:pverp))
       do i=ntoplw,pverp
          write(pname,'(a,i3.3)') 'NAL_absorp',i
          ierr = pio_def_var(File, trim(pname), pio_double, dimids(1:hdimcnt+1), abstot_desc(i))
       end do
	
       dimids(hdimcnt+1) = pverp_id
       ierr = pio_def_var(File, 'Emissivity', pio_double, dimids(1:hdimcnt+1), emstot_desc)

       dimids(hdimcnt+1) = pver_id
       do i=1,4
          write(pname,'(a,i3.3)') 'NN_absorp',i
          ierr = pio_def_var(File, pname, pio_double, dimids(1:hdimcnt+1), absnxt_desc(i))
       end do


    end if
    if (docosp) then
      ierr = pio_def_var(File, 'cosp_cnt_init', pio_int, cospcnt_desc)
    end if

    if (is_subcol_on()) then
      call subcol_init_restart(file, hdimids)
    end if
      
  end subroutine init_restart_physics

  subroutine write_restart_physics (File, cam_in, cam_out, pbuf2d)

      !-----------------------------------------------------------------------
      use physics_buffer,             only: physics_buffer_desc, pbuf_write_restart
      
      use chemistry,          only: chem_write_restart
      use prescribed_ozone,   only: write_prescribed_ozone_restart
      use prescribed_ghg,     only: write_prescribed_ghg_restart
      use prescribed_aero,    only: write_prescribed_aero_restart
      use prescribed_volcaero,only: write_prescribed_volcaero_restart
      use radiation,          only: radiation_do
      use cam_pio_utils,      only: get_phys_decomp, fillvalue
      use spmd_utils,         only: iam
      use subcol_utils,       only: is_subcol_on
      use subcol,             only: subcol_write_restart
      !
      ! Input arguments
      !
      type(file_desc_t), intent(inout) :: File
      type(cam_in_t),    intent(in)    :: cam_in(begchunk:endchunk)
      type(cam_out_t),   intent(in)    :: cam_out(begchunk:endchunk)
      type(physics_buffer_desc), pointer        :: pbuf2d(:,:)
      !
      ! Local workspace
      !
      real(r8), target :: tmpfield(pcols*(endchunk-begchunk+1))
      integer :: i, ii, j, m       ! loop index
      integer :: n3tmp             ! timestep index
      character(len=256) fname  ! abs-ems restart filename
      integer :: ioerr             ! I/O status
      integer :: ncol          ! number of vertical columns
      integer :: ierr
      type(io_desc_t), pointer :: iodesc
      integer(c_int64_t) :: write_touch_c
      !-----------------------------------------------------------------------

      write_touch_c = write_restart_physics_codon(1102_c_int64_t)
      call restart_physics_log_write_entry()

      ! Physics buffer
      if (is_subcol_on()) then
         call subcol_write_restart(File)
      end if

      call pbuf_write_restart(File, pbuf2d)

      if ( .not. adiabatic .and. .not. ideal_phys )then

         ! data for chemistry
         call chem_write_restart(File)

         call write_prescribed_ozone_restart(File)
         call write_prescribed_ghg_restart(File)
         call write_prescribed_aero_restart(File)
         call write_prescribed_volcaero_restart(File)
 
         call restart_physics_select_pack_impl()

	 do i=begchunk,endchunk
            ncol = cam_out(i)%ncol
            call restart_physics_fill_tail_wrap(ncol, fillvalue, fsnt(:,i), fsns(:,i), &
                 fsds(:,i), flnt(:,i), flns(:,i), landm(:,i), sgh(:,i), sgh30(:,i), &
                 trefmxav(:,i), trefmnav(:,i))
         end do

! the transfer intrinsic function fails if we are writting a 0 sized array, but the call to pio_write_darray 
! needs to be made because it is collective. 

         call get_phys_decomp(iodesc, 1,1,1,pio_double)
            
            ! Comsrf module variables (can following coup_csm definitions be removed?)
         call pio_write_darray(File, fsnt_desc, iodesc, fsnt, ierr)
         call pio_write_darray(File, fsns_desc, iodesc, fsns, ierr)
         call pio_write_darray(File, fsds_desc, iodesc, fsds, ierr)
         call pio_write_darray(File, flnt_desc, iodesc, flnt, ierr)
         
         call pio_write_darray(File, flns_desc,  iodesc,  flns, ierr)
         call pio_write_darray(File, landm_desc, iodesc, landm, ierr)
         call pio_write_darray(File, sgh_desc,   iodesc,   sgh, ierr)
         call pio_write_darray(File, sgh30_desc, iodesc, sgh30, ierr)
         
         call pio_write_darray(File, trefmxav_desc, iodesc, trefmxav, ierr)
         call pio_write_darray(File, trefmnav_desc, iodesc, trefmnav, ierr)

         call restart_physics_tmpfield_fill_wrap(size(tmpfield), fillvalue, tmpfield)
         do i=begchunk,endchunk
            ncol = cam_out(i)%ncol
            call restart_physics_pack_chunk_field_wrap(ncol, i-begchunk+1, cam_out(i)%flwds, tmpfield)
         end do
         call pio_write_darray(File, flwds_desc, iodesc, tmpfield, ierr)

         do i=begchunk,endchunk
            ncol = cam_out(i)%ncol
            call restart_physics_pack_chunk_field_wrap(ncol, i-begchunk+1, cam_out(i)%sols, tmpfield)
         end do
         call pio_write_darray(File, sols_desc, iodesc, tmpfield, ierr)

         do i=begchunk,endchunk
            ncol = cam_out(i)%ncol
            call restart_physics_pack_chunk_field_wrap(ncol, i-begchunk+1, cam_out(i)%soll, tmpfield)
         end do
         call pio_write_darray(File, soll_desc, iodesc, tmpfield, ierr)

         do i=begchunk,endchunk
            ncol = cam_out(i)%ncol
            call restart_physics_pack_chunk_field_wrap(ncol, i-begchunk+1, cam_out(i)%solsd, tmpfield)
         end do
         call pio_write_darray(File, solsd_desc, iodesc, tmpfield, ierr)

         do i=begchunk,endchunk
            ncol = cam_out(i)%ncol
            call restart_physics_pack_chunk_field_wrap(ncol, i-begchunk+1, cam_out(i)%solld, tmpfield)
         end do
         call pio_write_darray(File, solld_desc, iodesc, tmpfield, ierr)

         do i=begchunk,endchunk
            ncol = cam_out(i)%ncol
            call restart_physics_pack_chunk_field_wrap(ncol, i-begchunk+1, cam_out(i)%bcphidry, tmpfield)
         end do
         call pio_write_darray(File, bcphidry_desc, iodesc, tmpfield, ierr)

         do i=begchunk,endchunk
            ncol = cam_out(i)%ncol
            call restart_physics_pack_chunk_field_wrap(ncol, i-begchunk+1, cam_out(i)%bcphodry, tmpfield)
         end do
         call pio_write_darray(File, bcphodry_desc, iodesc, tmpfield, ierr)

         do i=begchunk,endchunk
            ncol = cam_out(i)%ncol
            call restart_physics_pack_chunk_field_wrap(ncol, i-begchunk+1, cam_out(i)%ocphidry, tmpfield)
         end do
         call pio_write_darray(File, ocphidry_desc, iodesc, tmpfield, ierr)

         do i=begchunk,endchunk
            ncol = cam_out(i)%ncol
            call restart_physics_pack_chunk_field_wrap(ncol, i-begchunk+1, cam_out(i)%ocphodry, tmpfield)
         end do
         call pio_write_darray(File, ocphodry_desc, iodesc, tmpfield, ierr)

         do i=begchunk,endchunk
            ncol = cam_out(i)%ncol
            call restart_physics_pack_chunk_field_wrap(ncol, i-begchunk+1, cam_out(i)%dstdry1, tmpfield)
         end do
         call pio_write_darray(File, dstdry1_desc, iodesc, tmpfield, ierr)

         do i=begchunk,endchunk
            ncol = cam_out(i)%ncol
            call restart_physics_pack_chunk_field_wrap(ncol, i-begchunk+1, cam_out(i)%dstdry2, tmpfield)
         end do
         call pio_write_darray(File, dstdry2_desc, iodesc, tmpfield, ierr)

         do i=begchunk,endchunk
            ncol = cam_out(i)%ncol
            call restart_physics_pack_chunk_field_wrap(ncol, i-begchunk+1, cam_out(i)%dstdry3, tmpfield)
         end do
         call pio_write_darray(File, dstdry3_desc, iodesc, tmpfield, ierr)

         do i=begchunk,endchunk
            ncol = cam_out(i)%ncol
            call restart_physics_pack_chunk_field_wrap(ncol, i-begchunk+1, cam_out(i)%dstdry4, tmpfield)
         end do
         call pio_write_darray(File, dstdry4_desc, iodesc, tmpfield, ierr)

         if (co2_transport()) then
            do i=begchunk,endchunk
               ncol = cam_out(i)%ncol
               call restart_physics_pack_chunk_field_wrap(ncol, i-begchunk+1, cam_out(i)%co2prog, tmpfield)
            end do
            call pio_write_darray(File, co2prog_desc, iodesc, tmpfield, ierr)
            do i=begchunk,endchunk
               ncol = cam_out(i)%ncol
               call restart_physics_pack_chunk_field_wrap(ncol, i-begchunk+1, cam_out(i)%co2diag, tmpfield)
            end do
            call pio_write_darray(File, co2diag_desc, iodesc, tmpfield, ierr)
         end if

         ! cam_in components
         do m = 1, pcnst
            do i = begchunk, endchunk
               ncol = cam_in(i)%ncol
               call restart_physics_pack_chunk_field_wrap(ncol, i-begchunk+1, cam_in(i)%cflx(:,m), tmpfield)
            end do
            call pio_write_darray(File, cflx_desc(m), iodesc, tmpfield, ierr)
         end do

      end if
      !
      !-----------------------------------------------------------------------
      ! Write the abs/ems restart dataset if necessary    
      !-----------------------------------------------------------------------
      !

      if ( radiation_do('aeres')  ) then
         
	 do i=begchunk,endchunk
            ncol = cam_out(i)%ncol
            if(ncol<pcols) then
               abstot_3d(ncol+1:pcols,:,:,i) = fillvalue
               absnxt_3d(ncol+1:pcols,:,:,i) = fillvalue
               emstot_3d(ncol+1:pcols,:,i) = fillvalue
            end if
         end do
         call get_phys_decomp(iodesc, 1,(pverp-ntoplw+1),1,pio_double)
         do i=ntoplw,pverp
            call pio_write_darray(File, abstot_desc(i), iodesc, abstot_3d(:,:,i,:), ierr)
         end do
         if(ntoplw/=1) then
            call get_phys_decomp(iodesc, 1,pverp,1,pio_double)
         end if
         call pio_write_darray(File, emstot_desc, iodesc, emstot_3d, ierr)
         call get_phys_decomp(iodesc, 1,pver,1,pio_double)

         do i=1,4
            call pio_write_darray(File, absnxt_desc(i), iodesc, absnxt_3d(:,:,i,:), ierr)
         end do

         deallocate(abstot_desc)
      end if

      if (docosp) then
        ierr = pio_put_var(File, cospcnt_desc, (/cosp_cnt(begchunk)/))
      end if

      
    end subroutine write_restart_physics

!#######################################################################

    subroutine read_restart_physics(File, cam_in, cam_out, pbuf2d)

     !-----------------------------------------------------------------------
     use physics_buffer,            only: physics_buffer_desc, pbuf_read_restart
     
     use chemistry,          only: chem_read_restart
     use cam_pio_utils,      only: get_phys_decomp, fillvalue
     use radiation,          only: radiation_do
     use prescribed_ozone,   only: read_prescribed_ozone_restart
     use prescribed_ghg,     only: read_prescribed_ghg_restart
     use prescribed_aero,    only: read_prescribed_aero_restart
     use prescribed_volcaero,only: read_prescribed_volcaero_restart
     use subcol_utils,       only: is_subcol_on
     use subcol,             only: subcol_read_restart
     !
     ! Arguments
     !
     type(file_desc_t),   intent(inout) :: File
     type(cam_in_t),            pointer :: cam_in(:)
     type(cam_out_t),           pointer :: cam_out(:)
     type(physics_buffer_desc), pointer :: pbuf2d(:,:)
     !
     ! Local workspace
     !
     real(r8), allocatable :: tmpfield(:)
     integer :: i, c, ii, m       ! loop index
     integer :: n3tmp             ! timestep index
     character*80  locfn          ! Local filename
     integer :: ioerr             ! I/O status
     integer :: ncol           ! number of columns in a chunk
     type(io_desc_t), pointer :: iodesc
     type(var_desc_t) :: vardesc
     integer :: ierr, csize, vsize
     character(len=4) :: num
     !-----------------------------------------------------------------------

     ! Allocate memory in physics buffer, buffer, comsrf, and radbuffer modules.
     ! (This is done in subroutine initial_conds for an initial run.)
     call initialize_comsrf()
     call initialize_radbuffer()

     ! Physics buffer

     ! subcol_read_restart must be called before pbuf_read_restart
     if (is_subcol_on()) then
        call subcol_read_restart(File)
     end if

     call pbuf_read_restart(File, pbuf2d)

     csize=endchunk-begchunk+1
     
     if ( .not. adiabatic .and. .not. ideal_phys )then

        ! data for chemistry
        call chem_read_restart(File)

        call read_prescribed_ozone_restart(File)
        call read_prescribed_ghg_restart(File)
        call read_prescribed_aero_restart(File)
        call read_prescribed_volcaero_restart(File)

        allocate(tmpfield(pcols*csize))
        tmpfield(:)=fillvalue

        call get_phys_decomp(iodesc,1,1,1,pio_double)

        ierr = pio_inq_varid(File, 'FSNT', vardesc)
        call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
        fsnt(:,:) = reshape(tmpfield, (/pcols, csize/))

        ierr = pio_inq_varid(File, 'FSNS', vardesc)
        call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
        fsns(:,:) = reshape(tmpfield, (/pcols, csize/))

        ierr = pio_inq_varid(File, 'FSDS', vardesc)
        call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
        fsds(:,:) = reshape(tmpfield, (/pcols, csize/))

        ierr = pio_inq_varid(File, 'FLNT', vardesc)
        call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
        flnt(:,:) = reshape(tmpfield, (/pcols, csize/))

        ierr = pio_inq_varid(File, 'FLNS', vardesc)
        call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
        flns(:,:) = reshape(tmpfield, (/pcols, csize/))

        ierr = pio_inq_varid(File, 'LANDM', vardesc)
        call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
        landm(:,:) = reshape(tmpfield, (/pcols, csize/))

        ierr = pio_inq_varid(File, 'SGH', vardesc)
        call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
        sgh(:,:) = reshape(tmpfield, (/pcols, csize/))

        ierr = pio_inq_varid(File, 'SGH30', vardesc)
        call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
        sgh30(:,:) = reshape(tmpfield, (/pcols, csize/))

        ierr = pio_inq_varid(File, 'TREFMXAV', vardesc)
        call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
        trefmxav(:,:) = reshape(tmpfield, (/pcols, csize/))

        ierr = pio_inq_varid(File, 'TREFMNAV', vardesc)
        call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
        trefmnav(:,:) = reshape(tmpfield, (/pcols, csize/))

        ierr = pio_inq_varid(File, 'FLWDS', vardesc)
        call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
        ii=0
        do c=begchunk,endchunk
           do i=1,pcols
              ii=ii+1
              cam_out(c)%flwds(i) = tmpfield(ii)
           end do
        end do

        ierr = pio_inq_varid(File, 'SOLS', vardesc)
        call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
        ii=0
        do c=begchunk,endchunk
           do i=1,pcols
              ii=ii+1
              cam_out(c)%sols(i) = tmpfield(ii)
           end do
        end do

        ierr = pio_inq_varid(File, 'SOLL', vardesc)
        call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
        ii=0
        do c=begchunk,endchunk
           do i=1,pcols
              ii=ii+1
              cam_out(c)%soll(i) = tmpfield(ii)
           end do
        end do

        ierr = pio_inq_varid(File, 'SOLSD', vardesc)
        call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
        ii=0
        do c=begchunk,endchunk
           do i=1,pcols
              ii=ii+1
              cam_out(c)%solsd(i) = tmpfield(ii)
           end do
        end do

        ierr = pio_inq_varid(File, 'SOLLD', vardesc)
        call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
        ii=0
        do c=begchunk,endchunk
           do i=1,pcols
              ii=ii+1
              cam_out(c)%solld(i) = tmpfield(ii)
           end do
        end do

        ierr = pio_inq_varid(File, 'BCPHIDRY', vardesc)
        call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
        ii=0
        do c=begchunk,endchunk
           do i=1,pcols
              ii=ii+1
              cam_out(c)%bcphidry(i) = tmpfield(ii)
           end do
        end do

        ierr = pio_inq_varid(File, 'BCPHODRY', vardesc)
        call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
        ii=0
        do c=begchunk,endchunk
           do i=1,pcols
              ii=ii+1
              cam_out(c)%bcphodry(i) = tmpfield(ii)
           end do
        end do

        ierr = pio_inq_varid(File, 'OCPHIDRY', vardesc)
        call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
        ii=0
        do c=begchunk,endchunk
           do i=1,pcols
              ii=ii+1
              cam_out(c)%ocphidry(i) = tmpfield(ii)
           end do
        end do

        ierr = pio_inq_varid(File, 'OCPHODRY', vardesc)
        call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
        ii=0
        do c=begchunk,endchunk
           do i=1,pcols
              ii=ii+1
              cam_out(c)%ocphodry(i) = tmpfield(ii)
           end do
        end do

        ierr = pio_inq_varid(File, 'DSTDRY1', vardesc)
        call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
        ii=0
        do c=begchunk,endchunk
           do i=1,pcols
              ii=ii+1
              cam_out(c)%dstdry1(i) = tmpfield(ii)
           end do
        end do

        ierr = pio_inq_varid(File, 'DSTDRY2', vardesc)
        call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
        ii=0
        do c=begchunk,endchunk
           do i=1,pcols
              ii=ii+1
              cam_out(c)%dstdry2(i) = tmpfield(ii)
           end do
        end do

        ierr = pio_inq_varid(File, 'DSTDRY3', vardesc)
        call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
        ii=0
        do c=begchunk,endchunk
           do i=1,pcols
              ii=ii+1
              cam_out(c)%dstdry3(i) = tmpfield(ii)
           end do
        end do

        ierr = pio_inq_varid(File, 'DSTDRY4', vardesc)
        call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
        ii=0
        do c=begchunk,endchunk
           do i=1,pcols
              ii=ii+1
              cam_out(c)%dstdry4(i) = tmpfield(ii)
           end do
        end do

        if (co2_transport()) then
           ierr = pio_inq_varid(File, 'CO2PROG', vardesc)
           call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
           ii=0
           do c=begchunk,endchunk
              do i=1,pcols
                 ii=ii+1
                 cam_out(c)%co2prog(i) = tmpfield(ii)
              end do
           end do

           ierr = pio_inq_varid(File, 'CO2DIAG', vardesc)
           call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
           ii=0
           do c=begchunk,endchunk
              do i=1,pcols
                 ii=ii+1
                 cam_out(c)%co2diag(i) = tmpfield(ii)
              end do
           end do
        end if

        ! Reading the CFLX* components from the restart is optional for
        ! backwards compatibility.  These fields were not needed for an
        ! exact restart until the UNICON scheme was added.  More generally,
        ! these components are only needed if they are not handled by the
        ! coupling layer restart (the ".rs." file), and if the values are
        ! used in the tphysbc physics before the tphysac code has a chance
        ! to update the values that are coming from boundary datasets.
        do m = 1, pcnst

           write(num,'(i4.4)') m

           call pio_seterrorhandling(File, PIO_BCAST_ERROR)
           ierr = pio_inq_varid(File, 'CFLX'//num, vardesc)
           call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)

           if (ierr == PIO_NOERR) then ! CFLX variable found on restart file
              call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
              ii = 0
              do c= begchunk, endchunk
                 do i = 1, pcols
                    ii = ii + 1
                    cam_in(c)%cflx(i,m) = tmpfield(ii)
                 end do
              end do
           end if

        end do

        deallocate(tmpfield)	

     end if

     !
     !-----------------------------------------------------------------------
     ! Read the abs/ems restart dataset if necessary    
     !-----------------------------------------------------------------------
     !
     if ( radiation_do('aeres')  ) then
        call pio_seterrorhandling( File, PIO_BCAST_ERROR)
        ierr = pio_inq_varid(File, 'Emissivity', vardesc)
        call pio_seterrorhandling( File, PIO_INTERNAL_ERROR)
        if(ierr/=PIO_NOERR) then
           if(masterproc) write(iulog,*) 'Warning: Emissivity variable not found on restart file.'
           return
        end if

        call get_phys_decomp(iodesc,1,pverp,1,pio_double)
        allocate(tmpfield(pcols*pverp*csize))         
        tmpfield(:)=fillvalue
        call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
        emstot_3d = reshape(tmpfield, (/pcols, pverp, csize/))
        
        vsize = pverp-ntoplw+1
        if(vsize/=pverp) then
           deallocate(tmpfield)
           call get_phys_decomp(iodesc,1,(pverp-ntoplw+1),1,pio_double)
           allocate(tmpfield(pcols*vsize*csize))         
        end if
        tmpfield(:)=fillvalue
        
        do i=ntoplw,pverp
           write(pname,'(a,i3.3)') 'NAL_absorp',i
           ierr = pio_inq_varid(File, trim(pname), vardesc)
           call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
           abstot_3d(:,:,i,:) = reshape(tmpfield, (/pcols, vsize, csize/))
        end do

        deallocate(tmpfield)
        call get_phys_decomp(iodesc,1,pver,1,pio_double)

        allocate(tmpfield(pcols*pver*csize))         
        tmpfield(:)=fillvalue

        do i=1,4
           write(pname,'(a,i3.3)') 'NN_absorp',i
           ierr = pio_inq_varid(File, trim(pname), vardesc)
           call pio_read_darray(File, vardesc, iodesc, tmpfield, ierr)
           absnxt_3d(:,:,i,:) = reshape(tmpfield, (/pcols, pver, csize/))
        end do
        deallocate(tmpfield)
     end if

     if (docosp) then
           call pio_seterrorhandling( File, PIO_BCAST_ERROR)
           ierr = pio_inq_varid(File, 'cosp_cnt_init', vardesc)
           call pio_seterrorhandling( File, PIO_INTERNAL_ERROR)
           if(ierr/=PIO_NOERR) then
             cosp_cnt_init=0
           else
             ierr = pio_get_var(File, vardesc, cosp_cnt_init)
           end if
     end if

   end subroutine read_restart_physics


   character(len=256) function get_abs_restart_filepath ( )
     !	
     ! Return the full filepath to the abs-ems restart file
     !	
     get_abs_restart_filepath = pname
   end function get_abs_restart_filepath

   subroutine restart_physics_select_pack_impl()
     character(len=32) :: impl_name
     integer :: n, status

     if (restart_pack_impl_selected) return

     call cam_codon_get_impl('RESTART_PHYSICS_PACK_IMPL', impl_name, n, status)
     if (status == 0 .and. n > 0) then
        use_native_restart_pack_impl = trim(adjustl(impl_name(:n))) == 'native'
     else
        use_native_restart_pack_impl = .false.
     end if

     if (masterproc) then
        if (use_native_restart_pack_impl) then
           write(iulog,*) 'restart_physics_pack implementation = native'
        else
           write(iulog,*) 'restart_physics_pack implementation = codon'
        end if
     end if

     restart_pack_impl_selected = .true.
   end subroutine restart_physics_select_pack_impl

   subroutine restart_physics_log_pack_entry()
     if (masterproc .and. .not. restart_pack_entered_logged) then
        write(iulog,*) 'restart_physics_pack entered (tail fill/tmpfield pack helpers = codon)'
        restart_pack_entered_logged = .true.
     end if
   end subroutine restart_physics_log_pack_entry

   subroutine restart_physics_fill_tail_wrap(ncol_local, fillvalue_local, fsnt_local, fsns_local, &
        fsds_local, flnt_local, flns_local, landm_local, sgh_local, sgh30_local, &
        trefmxav_local, trefmnav_local)
     integer, intent(in) :: ncol_local
     real(r8), intent(in) :: fillvalue_local
     real(r8), target, intent(inout) :: fsnt_local(pcols), fsns_local(pcols), fsds_local(pcols)
     real(r8), target, intent(inout) :: flnt_local(pcols), flns_local(pcols), landm_local(pcols)
     real(r8), target, intent(inout) :: sgh_local(pcols), sgh30_local(pcols)
     real(r8), target, intent(inout) :: trefmxav_local(pcols), trefmnav_local(pcols)

     interface
        subroutine restart_physics_fill_tail_codon(ncol_c, pcols_c, fillvalue_c, &
             fsnt_p, fsns_p, fsds_p, flnt_p, flns_p, landm_p, sgh_p, sgh30_p, &
             trefmxav_p, trefmnav_p) bind(c, name="restart_physics_fill_tail_codon")
          import c_double, c_int64_t, c_ptr
          integer(c_int64_t), value :: ncol_c, pcols_c
          real(c_double), value :: fillvalue_c
          type(c_ptr), value :: fsnt_p, fsns_p, fsds_p, flnt_p, flns_p, landm_p
          type(c_ptr), value :: sgh_p, sgh30_p, trefmxav_p, trefmnav_p
        end subroutine restart_physics_fill_tail_codon
     end interface

     if (use_native_restart_pack_impl) then
        if (ncol_local < pcols) then
           fsnt_local(ncol_local+1:pcols) = fillvalue_local
           fsns_local(ncol_local+1:pcols) = fillvalue_local
           fsds_local(ncol_local+1:pcols) = fillvalue_local
           flnt_local(ncol_local+1:pcols) = fillvalue_local
           flns_local(ncol_local+1:pcols) = fillvalue_local
           landm_local(ncol_local+1:pcols) = fillvalue_local
           sgh_local(ncol_local+1:pcols) = fillvalue_local
           sgh30_local(ncol_local+1:pcols) = fillvalue_local
           trefmxav_local(ncol_local+1:pcols) = fillvalue_local
           trefmnav_local(ncol_local+1:pcols) = fillvalue_local
        end if
     else
        call restart_physics_log_pack_entry()
        call restart_physics_fill_tail_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), &
             real(fillvalue_local, c_double), c_loc(fsnt_local(1)), c_loc(fsns_local(1)), &
             c_loc(fsds_local(1)), c_loc(flnt_local(1)), c_loc(flns_local(1)), &
             c_loc(landm_local(1)), c_loc(sgh_local(1)), c_loc(sgh30_local(1)), &
             c_loc(trefmxav_local(1)), c_loc(trefmnav_local(1)))
     end if
   end subroutine restart_physics_fill_tail_wrap

   subroutine restart_physics_tmpfield_fill_wrap(total_len_local, fillvalue_local, tmpfield_local)
     integer, intent(in) :: total_len_local
     real(r8), intent(in) :: fillvalue_local
     real(r8), target, intent(inout) :: tmpfield_local(total_len_local)

     interface
        subroutine restart_physics_tmpfield_fill_codon(total_len_c, fillvalue_c, tmpfield_p) &
             bind(c, name="restart_physics_tmpfield_fill_codon")
          import c_double, c_int64_t, c_ptr
          integer(c_int64_t), value :: total_len_c
          real(c_double), value :: fillvalue_c
          type(c_ptr), value :: tmpfield_p
        end subroutine restart_physics_tmpfield_fill_codon
     end interface

     if (use_native_restart_pack_impl) then
        tmpfield_local(:) = fillvalue_local
     else
        call restart_physics_log_pack_entry()
        call restart_physics_tmpfield_fill_codon(int(total_len_local, c_int64_t), &
             real(fillvalue_local, c_double), c_loc(tmpfield_local(1)))
     end if
   end subroutine restart_physics_tmpfield_fill_wrap

   subroutine restart_physics_pack_chunk_field_wrap(ncol_local, chunk_pos_local, field_local, tmpfield_local)
     integer, intent(in) :: ncol_local, chunk_pos_local
     real(r8), target, intent(in) :: field_local(pcols)
     real(r8), target, intent(inout) :: tmpfield_local(pcols*(endchunk-begchunk+1))
     integer :: j, offset

     interface
        subroutine restart_physics_pack_chunk_field_codon(ncol_c, pcols_c, chunk_pos_c, field_p, tmpfield_p) &
             bind(c, name="restart_physics_pack_chunk_field_codon")
          import c_int64_t, c_ptr
          integer(c_int64_t), value :: ncol_c, pcols_c, chunk_pos_c
          type(c_ptr), value :: field_p, tmpfield_p
        end subroutine restart_physics_pack_chunk_field_codon
     end interface

     if (use_native_restart_pack_impl) then
        offset = (chunk_pos_local - 1) * pcols
        do j = 1, pcols
           if (j <= ncol_local) tmpfield_local(offset + j) = field_local(j)
        end do
     else
        call restart_physics_log_pack_entry()
        call restart_physics_pack_chunk_field_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), &
             int(chunk_pos_local, c_int64_t), c_loc(field_local(1)), c_loc(tmpfield_local(1)))
     end if
   end subroutine restart_physics_pack_chunk_field_wrap

 end module restart_physics
