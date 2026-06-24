!--------------------------------------------------------------------------------
! Manages writing reaction rates to history
!--------------------------------------------------------------------------------
module rate_diags

  use shr_kind_mod,     only : r8 => shr_kind_r8
  use shr_kind_mod,     only : CL => SHR_KIND_CL, CX => SHR_KIND_CX
  use cam_history,      only : fieldname_len
  use cam_history,      only : addfld,phys_decomp
  use cam_history,      only : outfld
  use chem_mods,        only : rxt_tag_cnt, rxt_tag_lst, rxt_tag_map
  use ppgrid,           only : pver
  use spmd_utils,       only : masterproc
  use cam_abortutils,   only : endrun
  use cam_logfile,      only : iulog
  use iso_c_binding,    only : c_int64_t, c_loc, c_ptr

  implicit none
  private 
  public :: rate_diags_init
  public :: rate_diags_calc
  public :: rate_diags_readnl

  character(len=fieldname_len) :: rate_names(rxt_tag_cnt)

  type rate_grp_t
    character(len=24) :: name
    integer :: nmembers = 0
    integer, allocatable :: map(:)
    real(r8), allocatable :: multipler(:)
  endtype rate_grp_t

  integer :: ngrps = 0
  type(rate_grp_t), allocatable :: grps(:)  

  integer, parameter :: maxsums = 100
  character(len=CX) :: rxn_rate_sums(maxsums) = ' '
  logical :: rate_diags_batch_use_native_impl = .false.
  logical :: rate_diags_batch_impl_selected = .false.
  logical :: rate_diags_batch_entered_logged = .false.
  logical :: rate_diags_calc_entered_logged = .false.
  logical :: rate_diags_init_proof_written = .false.
  logical :: parse_rate_sums_proof_written = .false.

  interface
    function rate_diags_init_codon(tag_len_c, fieldname_len_c, rxt_tag_cnt_c, &
         rxt_tag_ascii_p, rate_name_ascii_p) result(out_c) bind(c, name="rate_diags_init_codon")
      use iso_c_binding, only : c_int64_t, c_ptr
      integer(c_int64_t), value :: tag_len_c, fieldname_len_c, rxt_tag_cnt_c
      type(c_ptr), value :: rxt_tag_ascii_p, rate_name_ascii_p
      integer(c_int64_t) :: out_c
    end function rate_diags_init_codon

    function parse_rate_sums_codon(active) result(out_c) bind(c, name="parse_rate_sums_codon")
      use iso_c_binding, only : c_int64_t
      integer(c_int64_t), value :: active
      integer(c_int64_t) :: out_c
    end function parse_rate_sums_codon
  end interface

contains

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine rate_diags_batch_append_proof(proof_line)

    character(len=*), intent(in) :: proof_line
    character(len=512) :: proof_file
    integer :: status, n, unitno

    proof_file = ''
    call get_environment_variable('RATE_DIAGS_BATCH_PROOF_FILE', value=proof_file, length=n, status=status)
    if (status == 0 .and. n > 0) then
       open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
       write(unitno,'(A)') trim(proof_line)
       close(unitno)
    end if

  end subroutine rate_diags_batch_append_proof

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine rate_diags_batch_select_impl()

    use cam_logfile, only : iulog

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (rate_diags_batch_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('RATE_DIAGS_BATCH_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       rate_diags_batch_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       rate_diags_batch_use_native_impl = .false.
    end if

    rate_diags_batch_impl_selected = .true.

    if (masterproc) then
       if (rate_diags_batch_use_native_impl) then
          write(iulog,*) 'rate_diags_batch implementation = native'
          call rate_diags_batch_append_proof('rate_diags_batch selector entered implementation = native')
       else
          write(iulog,*) 'rate_diags_batch implementation = codon'
          call rate_diags_batch_append_proof('rate_diags_batch selector entered implementation = codon')
       end if
       call flush(iulog)
    end if

  end subroutine rate_diags_batch_select_impl

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine rate_diags_batch_log_entered()

    use cam_logfile, only : iulog

    if (rate_diags_batch_entered_logged) return
    rate_diags_batch_entered_logged = .true.

    if (masterproc) then
       write(iulog,*) 'rate_diags_batch entered (unified set-rates/tagged-conversion stage dispatch = codon)'
       call rate_diags_batch_append_proof('rate_diags_batch entered (unified set-rates/tagged-conversion stage dispatch = codon)')
       call flush(iulog)
    end if

  end subroutine rate_diags_batch_log_entered

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine rate_diags_calc_log_entered()

    use cam_logfile, only : iulog

    if (rate_diags_calc_entered_logged) return
    rate_diags_calc_entered_logged = .true.

    if (masterproc) then
       write(iulog,*) 'rate_diags_calc direct = codon; set-rates, tagged conversion, and group rates = codon'
       call rate_diags_batch_append_proof( &
            'rate_diags_calc direct = codon; set-rates, tagged conversion, and group rates = codon')
       call flush(iulog)
    end if

  end subroutine rate_diags_calc_log_entered

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine rate_diags_readnl(nlfile)

    use iso_c_binding, only : c_int64_t
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand
    use mo_util,         only: chemistry_misc_codon_touch

    interface
       function rate_diags_readnl_codon(tag) result(tag_out) bind(c, name='rate_diags_readnl_codon')
         import :: c_int64_t
         integer(c_int64_t), value :: tag
         integer(c_int64_t) :: tag_out
       end function rate_diags_readnl_codon
    end interface

    character(len=32) :: rt_codon_impl_name
    integer :: rt_codon_n, rt_codon_status
    integer(c_int64_t) :: rt_codon_tag_out
    logical, save :: rt_codon_proof_seen = .false.


    ! args 
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr

    namelist /rxn_rate_diags_nl/ rxn_rate_sums

    rt_codon_impl_name = 'codon'
    call cam_codon_get_impl('RATE_DIAGS_READNL_IMPL', rt_codon_impl_name, rt_codon_n, rt_codon_status)
    if (.not. (rt_codon_status == 0 .and. rt_codon_n > 0 .and. &
         trim(adjustl(rt_codon_impl_name(:rt_codon_n))) == 'native')) then
       rt_codon_tag_out = rate_diags_readnl_codon(int(174, c_int64_t))
       if (rt_codon_tag_out /= int(174, c_int64_t)) then
          write(iulog,*) 'rate_diags_readnl_codon tag roundtrip failed'
          stop 2
       endif
       if (.not. rt_codon_proof_seen) then
          write(iulog,*) 'rate_diags_readnl implementation = codon'
          rt_codon_proof_seen = .true.
       endif
    endif

    ! Read namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'rxn_rate_diags_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, rxn_rate_diags_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun('rate_diags_readnl:: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(rxn_rate_sums,len(rxn_rate_sums(1))*maxsums, mpichar, 0, mpicom)
#endif
  end subroutine rate_diags_readnl
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
  subroutine rate_diags_init

    integer :: i, ichar
    integer(c_int64_t) :: codon_count
    integer(c_int64_t), target :: rxt_tag_ascii(len(rxt_tag_lst), max(1, rxt_tag_cnt))
    integer(c_int64_t), target :: rate_name_ascii(fieldname_len, max(1, rxt_tag_cnt))

    do i = 1, rxt_tag_cnt
       do ichar = 1, len(rxt_tag_lst)
          rxt_tag_ascii(ichar,i) = int(iachar(rxt_tag_lst(i)(ichar:ichar)), c_int64_t)
       end do
    end do

    codon_count = rate_diags_init_codon(int(len(rxt_tag_lst), c_int64_t), int(fieldname_len, c_int64_t), &
         int(rxt_tag_cnt, c_int64_t), c_loc(rxt_tag_ascii(1,1)), c_loc(rate_name_ascii(1,1)))
    if (codon_count /= int(rxt_tag_cnt, c_int64_t)) then
       call endrun('rate_diags_init: Codon rate-name count mismatch')
    end if

    if (masterproc .and. .not. rate_diags_init_proof_written) then
       write(iulog,'(A)') 'rate_diags_init direct = codon; rate-name construction direct; addfld native CAM API boundary'
       call rate_diags_batch_append_proof('rate_diags_init direct = codon; rate-name construction direct; addfld native CAM API boundary')
       rate_diags_init_proof_written = .true.
       call flush(iulog)
    end if

    do i = 1,rxt_tag_cnt
       rate_names(i) = ' '
       do ichar = 1, fieldname_len
          rate_names(i)(ichar:ichar) = achar(int(rate_name_ascii(ichar,i)))
       end do
       call addfld(rate_names(i), 'molecules/cm3/sec', pver,'A','reaction rate', phys_decomp)
    enddo

    call parse_rate_sums()

    do i = 1, ngrps
       call addfld( grps(i)%name, 'molecules/cm3/sec', pver,'A','reaction rate group', phys_decomp)
    enddo

  end subroutine rate_diags_init

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
  subroutine rate_diags_calc( rxt_rates, vmr, m, ncol, lchnk )

    use mo_rxt_rates_conv, only: set_rates
    use chem_mods, only : rxntot
    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    real(r8), target, intent(inout) :: rxt_rates(:,:,:) ! 'molec/cm3/sec'
    real(r8), target, intent(in)    :: vmr(:,:,:)
    real(r8), target, intent(in)    :: m(:,:)           ! air density (molecules/cm3)
    integer,  intent(in)    :: ncol, lchnk

    integer :: i, j, max_group_members
    logical :: used_codon_calc
    integer(c_int64_t), target :: rxt_tag_map_i64(max(1,rxt_tag_cnt))
    integer(c_int64_t), allocatable, target :: grp_nm_i64(:)
    integer(c_int64_t), allocatable, target :: grp_map_i64(:,:)
    real(r8), allocatable, target :: grp_mult(:,:)
    real(r8), allocatable, target :: group_rates(:,:,:)
    real(r8) :: group_rate(ncol,pver)

    interface
       subroutine rate_diags_calc_codon(ncol_c, pver_c, rxntot_c, rxt_tag_cnt_c, &
            rxt_rates_p, vmr_p, m_p, rxt_tag_map_p, ngrps_c, max_group_members_c, &
            grp_nm_p, grp_map_p, grp_mult_p, group_rates_p) bind(c, name="rate_diags_calc_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c, rxntot_c, rxt_tag_cnt_c
         integer(c_int64_t), value :: ngrps_c, max_group_members_c
         type(c_ptr), value :: rxt_rates_p, vmr_p, m_p, rxt_tag_map_p
         type(c_ptr), value :: grp_nm_p, grp_map_p, grp_mult_p, group_rates_p
       end subroutine rate_diags_calc_codon
    end interface

    used_codon_calc = .false.

    call rate_diags_batch_select_impl()

    if (rate_diags_batch_use_native_impl) then
       call set_rates( rxt_rates, vmr, ncol )

       ! output individual tagged rates
       do i = 1, rxt_tag_cnt
          ! convert from vmr/sec to molecules/cm3/sec
          rxt_rates(:ncol,:,rxt_tag_map(i)) = rxt_rates(:ncol,:,rxt_tag_map(i)) * m(:ncol,:)
       enddo
    else
       call rate_diags_batch_log_entered()
       do i = 1, rxt_tag_cnt
          rxt_tag_map_i64(i) = int(rxt_tag_map(i), c_int64_t)
       end do
       max_group_members = 1
       do i = 1, ngrps
          max_group_members = max(max_group_members, grps(i)%nmembers)
       end do
       allocate(grp_nm_i64(max(1, ngrps)))
       allocate(grp_map_i64(max_group_members, max(1, ngrps)))
       allocate(grp_mult(max_group_members, max(1, ngrps)))
       allocate(group_rates(ncol, pver, max(1, ngrps)))
       grp_nm_i64(:) = 0_c_int64_t
       grp_map_i64(:,:) = 0_c_int64_t
       grp_mult(:,:) = 0._r8
       group_rates(:,:,:) = 0._r8
       do i = 1, ngrps
          grp_nm_i64(i) = int(grps(i)%nmembers, c_int64_t)
          do j = 1, grps(i)%nmembers
             grp_map_i64(j,i) = int(grps(i)%map(j), c_int64_t)
             grp_mult(j,i) = grps(i)%multipler(j)
          end do
       end do
       call rate_diags_calc_log_entered()
       call rate_diags_calc_codon(int(ncol, c_int64_t), int(pver, c_int64_t), int(rxntot, c_int64_t), &
            int(rxt_tag_cnt, c_int64_t), c_loc(rxt_rates), c_loc(vmr), c_loc(m), c_loc(rxt_tag_map_i64), &
            int(ngrps, c_int64_t), int(max_group_members, c_int64_t), c_loc(grp_nm_i64(1)), &
            c_loc(grp_map_i64(1,1)), c_loc(grp_mult(1,1)), c_loc(group_rates(1,1,1)))
       used_codon_calc = .true.
    end if

    do i = 1, rxt_tag_cnt
       call outfld( rate_names(i), rxt_rates(:ncol,:,rxt_tag_map(i)), ncol, lchnk )
    enddo

    ! output rate groups ( or families )
    do i = 1, ngrps
       if (used_codon_calc) then
          call outfld( grps(i)%name, group_rates(:,:,i), ncol, lchnk )
       else
          group_rate(:,:) = 0._r8
          do j = 1, grps(i)%nmembers
            group_rate(:ncol,:) = group_rate(:ncol,:) + grps(i)%multipler(j)*rxt_rates(:ncol,:,grps(i)%map(j))
          enddo
          call outfld( grps(i)%name, group_rate(:ncol,:), ncol, lchnk )
       end if
    end do

    if (allocated(group_rates)) deallocate(group_rates)
    if (allocated(grp_mult)) deallocate(grp_mult)
    if (allocated(grp_map_i64)) deallocate(grp_map_i64)
    if (allocated(grp_nm_i64)) deallocate(grp_nm_i64)

  end subroutine rate_diags_calc

!-------------------------------------------------------------------
! Private routines :
!-------------------------------------------------------------------
!-------------------------------------------------------------------
  
  subroutine parse_rate_sums

    integer :: ndxs(512)
    integer :: nelem, spc_len, i,j,k, rxt_ndx
    character(len=CL) :: tmp_str, tmp_name

    character(len=8) :: xchr ! multipler
    real(r8) :: xdbl

    character(len=CX) :: sum_string
    integer(c_int64_t) :: active_c

    ! a group is  a sum of reaction rates 

    ! count the numger of sums (or groups)
    sumcnt: do i = 1,maxsums
       spc_len=len_trim(rxn_rate_sums(i))
       if ( spc_len > 0 ) then
          ngrps = ngrps+1
       else
          exit sumcnt
       endif
    enddo sumcnt

    active_c = parse_rate_sums_codon(merge(1_c_int64_t, 0_c_int64_t, ngrps > 0))
    if (.not. parse_rate_sums_proof_written) then
       parse_rate_sums_proof_written = .true.
       if (masterproc) then
          if (active_c == 0_c_int64_t) then
             write(iulog,'(A)') 'parse_rate_sums direct = codon empty-rate-sums no-op'
          else
             write(iulog,'(A)') 'parse_rate_sums selector = codon; active rate-sum parser body = native'
          end if
          call flush(iulog)
       end if
    end if

    if (active_c == 0_c_int64_t) return

    ! parse the individual sum strings...  and form the groupings
    has_grps: if (ngrps>0) then

       allocate( grps(ngrps) )

       ! from shr_megan_mod ... should be generalized and shared...
       grploop: do i = 1,ngrps

          ! parse out the rxn names and multipliers
          ! from first parsing out the terms in the summation equation ("+" separates the terms)

          sum_string = rxn_rate_sums(i)
          j = scan( sum_string, '=' )
          nelem = 1
          ndxs(nelem) = j ! ndxs stores the index of each term of the equation

          ! find indices of all the terms in the equation
          tmp_str = trim( sum_string(j+1:) )
          j = scan( tmp_str, '+' )
          do while(j>0)
             nelem = nelem+1
             ndxs(nelem) = ndxs(nelem-1) + j
             tmp_str = tmp_str(j+1:)
             j = scan( tmp_str, '+' )
          enddo
          ndxs(nelem+1) = len(sum_string)+1

          grps(i)%nmembers = nelem ! number of terms 
          grps(i)%name =  trim(adjustl( sum_string(:ndxs(1)-1))) ! thing to the left of the "=" is used as the name of the group

          ! now that we have the number of terms in the summation allocate memory for the map (reaction indices) and multipliers
          allocate(grps(i)%map(nelem)) 
          allocate(grps(i)%multipler(nelem))

          ! now parse out the  rxn names and multiplers from the terms 
          elmloop: do k = 1,nelem
             grps(i)%multipler(k) = 1._r8
             ! get the rxn name which follows the '*' operator if the is one
             tmp_name = adjustl(sum_string(ndxs(k)+1:ndxs(k+1)-1))
             j = scan( tmp_name, '*' )
             if (j>0) then
                xchr = tmp_name(1:j-1) ! get the multipler (left of the '*')
                read( xchr, * ) xdbl   ! convert the string to a real
                grps(i)%multipler(k) = xdbl ! store the multiplier
                tmp_name = adjustl(tmp_name(j+1:)) ! get the rxn name (right of the '*')
             endif
             ! look up the corresponding reaction index ...
             rxt_ndx = lookup_tag_ndx( tmp_name )
             if ( rxt_ndx > 0 ) then
                grps(i)%map(k) = rxt_ndx
             else
                call endrun('rate_diags::parse_rate_sums rate name not found : '//trim(tmp_name))
             endif
          enddo elmloop
       enddo grploop
    endif has_grps

  end subroutine parse_rate_sums

!-------------------------------------------------------------------
! finds the index corresponging to a given reacton name
!-------------------------------------------------------------------
  function lookup_tag_ndx( name ) result( ndx )
    character(len=*) :: name
    integer :: ndx

    integer :: i

    ndx = -1

    findloop: do i = 1,rxt_tag_cnt
       if (trim(name) .eq. trim(rate_names(i)(3:))) then
          ndx = i
          return
       endif
    end do findloop

  end function lookup_tag_ndx

end module rate_diags
