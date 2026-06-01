
      module charge_neutrality

      use shr_kind_mod,      only : r8 => shr_kind_r8
      use iso_c_binding,     only : c_int64_t
      use cam_logfile,       only : iulog

      implicit none

      private
      public :: charge_balance
      public :: charge_fix     ! temporary, for fixing charge balance after vertical diffusion
                               ! without converting mass mixing ratios to volume
                               ! mean mass assumed to be mwdry
      logical :: charge_balance_use_native_impl = .false.
      logical :: charge_balance_impl_selected = .false.
      logical :: charge_fix_proof_written = .false.

      interface
         function charge_fix_active_codon(active) result(out_c) bind(c, name="charge_fix_active_codon")
           import :: c_int64_t
           integer(c_int64_t), value :: active
           integer(c_int64_t) :: out_c
         end function charge_fix_active_codon
      end interface

      contains

      subroutine charge_balance( ncol, conc )
!-----------------------------------------------------------------------      
!        ... force ion/electron balance
!-----------------------------------------------------------------------      

        use ppgrid,       only : pver
        use mo_chem_utls, only : get_spc_ndx
        use chem_mods,    only : gas_pcnst
        use iso_c_binding, only : c_int64_t, c_loc, c_ptr

        implicit none

!-----------------------------------------------------------------------      
!        ... dummy arguments
!-----------------------------------------------------------------------      
      integer,  intent(in)          :: ncol
      real(r8), target, intent(inout) :: conc(ncol,pver,gas_pcnst)         ! concentration

!-----------------------------------------------------------------------      
!        ... local variables
!-----------------------------------------------------------------------      
      integer  :: k, n
      integer  :: elec_ndx, np_ndx, n2p_ndx, op_ndx, o2p_ndx, nop_ndx
      real(r8), target :: wrk(ncol,pver)

      interface
         subroutine charge_balance_codon(ncol_c, pver_c, gas_pcnst_c, np_ndx_c, n2p_ndx_c, op_ndx_c, o2p_ndx_c, &
              nop_ndx_c, conc_p, wrk_p) bind(c, name="charge_balance_codon")
           use iso_c_binding, only : c_int64_t, c_ptr
           integer(c_int64_t), value :: ncol_c, pver_c, gas_pcnst_c
           integer(c_int64_t), value :: np_ndx_c, n2p_ndx_c, op_ndx_c, o2p_ndx_c, nop_ndx_c
           type(c_ptr), value :: conc_p, wrk_p
         end subroutine charge_balance_codon
      end interface

      call charge_balance_select_impl()

      elec_ndx = get_spc_ndx('e')
#ifdef CB_DEBUG
      write(iulog,*) ' '
      write(iulog,*) '------------------------------------------------------------------'
      write(iulog,*) 'charge_balance: e ndx,offset = ',elec_ndx,offset
      write(iulog,*) 'charge_balance: size of conc = ',size(conc,dim=1),' x ',size(conc,dim=2),' x ',size(conc,dim=3)
#endif
      if (charge_balance_use_native_impl) then
         if( elec_ndx > 0 ) then
	    wrk(:,:) = 0._r8
            n = get_spc_ndx('Np')
            if( n > 0 ) then
	       do k = 1,pver
	         wrk(:,k) = wrk(:,k) + conc(:ncol,k,n)
	       end do
            end if
            n = get_spc_ndx('N2p')
            if( n > 0 ) then
	       do k = 1,pver
	         wrk(:,k) = wrk(:,k) + conc(:ncol,k,n)
	       end do
            end if
            n = get_spc_ndx('Op')
            if( n > 0 ) then
	       do k = 1,pver
	         wrk(:,k) = wrk(:,k) + conc(:ncol,k,n)
	       end do
            end if
            n = get_spc_ndx('O2p')
            if( n > 0 ) then
	       do k = 1,pver
	         wrk(:,k) = wrk(:,k) + conc(:ncol,k,n)
	       end do
            end if
            n = get_spc_ndx('NOp')
            if( n > 0 ) then
	       do k = 1,pver
	         wrk(:,k) = wrk(:,k) + conc(:ncol,k,n)
	       end do
            end if
#ifdef CB_DEBUG
            write(iulog,*) 'charge_balance: electron concentration before balance'
            write(iulog,'(1p,5g15.7)') conc(1,:,elec_ndx)
            write(iulog,*) 'charge_balance: electron concentration after  balance'
            write(iulog,'(1p,5g15.7)') wrk(1,:)
            write(iulog,*) '------------------------------------------------------------------'
            write(iulog,*) ' '
#endif
            conc(:ncol,:,elec_ndx) = wrk(:ncol,:)
         end if
         return
      end if

      if( elec_ndx > 0 ) then
	 wrk(:,:) = 0._r8
         np_ndx = get_spc_ndx('Np')
         n2p_ndx = get_spc_ndx('N2p')
         op_ndx = get_spc_ndx('Op')
         o2p_ndx = get_spc_ndx('O2p')
         nop_ndx = get_spc_ndx('NOp')

         call charge_balance_codon( &
              int(ncol, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), int(np_ndx, c_int64_t), &
              int(n2p_ndx, c_int64_t), int(op_ndx, c_int64_t), int(o2p_ndx, c_int64_t), int(nop_ndx, c_int64_t), &
              c_loc(conc), c_loc(wrk) &
         )

         conc(:ncol,:,elec_ndx) = wrk(:ncol,:)
      end if

      end subroutine charge_balance

      subroutine charge_balance_select_impl()

      use spmd_utils, only : masterproc

      character(len=32) :: impl_name
      integer :: status, n, i, code

      if (charge_balance_impl_selected) return

      impl_name = 'codon'
      call get_environment_variable('CHARGE_BALANCE_IMPL', value=impl_name, length=n, status=status)

      if (status == 0 .and. n > 0) then
         do i = 1, n
            code = iachar(impl_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
         end do
         charge_balance_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
      else
         charge_balance_use_native_impl = .false.
      end if

      charge_balance_impl_selected = .true.

      if (masterproc) then
         if (charge_balance_use_native_impl) then
            write(iulog,*) 'charge_balance implementation = native'
         else
            write(iulog,*) 'charge_balance implementation = codon'
         end if
         call flush(iulog)
      end if

      end subroutine charge_balance_select_impl


      subroutine charge_fix(state, pbuf)
!-----------------------------------------------------------------------      
!        ... force ion/electron balance
!-----------------------------------------------------------------------      

      use ppgrid,              only : pcols, pver
      use constituents,        only : cnst_get_ind, cnst_mw
      use physconst,           only : mwdry                   ! molecular weight of dry air
      use physconst,           only : mbarv                       ! Constituent dependent mbar
      use phys_control,        only : waccmx_is
      use short_lived_species, only : slvd_index,slvd_pbf_ndx => pbf_idx ! Routines to access short lived species in pbuf
      use mo_chem_utls,        only : get_spc_ndx
      use chem_mods,           only : adv_mass
      use physics_buffer,      only : pbuf_get_field,physics_buffer_desc ! Needed to get variables from physics buffer
      use physics_types,       only : physics_state
      use spmd_utils,          only : masterproc

      implicit none

!-----------------------------------------------------------------------      
!        ... dummy arguments
!-----------------------------------------------------------------------      
      type(physics_state), intent(inout), target :: state
      type(physics_buffer_desc), pointer :: pbuf(:)    ! physics buffer

!-----------------------------------------------------------------------      
!        ... local variables
!-----------------------------------------------------------------------      
      integer  :: k, n, ns, nc
      integer  :: elec_ndx, elec_sndx
      integer  :: lchnk                 !Chunk number from state structure
      integer  :: ncol                  !Number of columns in this chunk from state structure

      real(r8) :: wrk(pcols,pver)
      real(r8) :: mbar(pcols,pver)  ! mean mass (=mwdry) used to fake out optimizer to get
                                    ! identical answers to old code
				   
      real(r8), dimension(:,:,:), pointer   :: q         ! model mass mixing ratios
      real(r8), dimension(:,:),   pointer   :: qs        ! Pointer to access fields in pbuf
      real(r8), dimension(:,:),   pointer   :: qse       ! Pointer to access electrons in pbuf
      integer(c_int64_t) :: active_c

!-----------------------------------------------------------------------
      lchnk = state%lchnk
      ncol  = state%ncol
      q => state%q

     !-----------------------------------------------------
     ! Get index to access electron mass mixing ratio
     !-----------------------------------------------------
      elec_sndx = -1
      call cnst_get_ind( 'e', elec_ndx, abort=.false. )
      if (elec_ndx < 0) elec_sndx = slvd_index( 'e' )
      active_c = charge_fix_active_codon(merge(1_c_int64_t, 0_c_int64_t, &
           elec_ndx > 0 .or. elec_sndx > 0))
      if (.not. charge_fix_proof_written) then
         charge_fix_proof_written = .true.
         if (masterproc) then
            if (active_c == 0_c_int64_t) then
               write(iulog,'(A)') 'charge_fix direct = codon no-electron-species no-op'
            else
               write(iulog,'(A)') 'charge_fix selector = codon; active electron/ion balance body = native'
            end if
            call flush(iulog)
         end if
      end if
      if (active_c == 0_c_int64_t) return

     !------------------------------------------------
     ! assume that mbar = mwdry except for WACCM-X
     !------------------------------------------------
      mbar(:ncol,:) = mwdry
      if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) mbar(:ncol,:) = mbarv(:ncol,:,lchnk)

      !--------------------------------------------------------------------
      ! If electrons are in state%q or pbuf, add up ions to get electrons
      !--------------------------------------------------------------------
      if( elec_ndx > 0 .or. elec_sndx > 0) then
	 wrk(:,:) = 0._r8
         call cnst_get_ind( 'Np', n, abort=.false. )
         if (n < 0) then 
            ns = slvd_index( 'Np' )
	    call pbuf_get_field(pbuf, slvd_pbf_ndx, qs, start=(/1,1,ns/), kount=(/pcols,pver,1/) )
            nc = get_spc_ndx('Np')
      	    do k = 1,pver
	      wrk(:ncol,k) = wrk(:ncol,k) + mbar(:ncol,k) * qs(:ncol,k) / adv_mass(nc)
	    end do	   
         else
	    do k = 1,pver
	      wrk(:ncol,k) = wrk(:ncol,k) + mbar(:ncol,k) * q(:ncol,k,n) / cnst_mw(n)
	    end do
         end if
         call cnst_get_ind( 'N2p', n, abort=.false. )
         if (n < 0) then 
            ns = slvd_index( 'N2p' )
	    call pbuf_get_field(pbuf, slvd_pbf_ndx, qs, start=(/1,1,ns/), kount=(/pcols,pver,1/) )
            nc = get_spc_ndx('N2p')
      	    do k = 1,pver
	      wrk(:ncol,k) = wrk(:ncol,k) + mbar(:ncol,k) * qs(:ncol,k) / adv_mass(nc)
	    end do	   
         else
	    do k = 1,pver
	      wrk(:ncol,k) = wrk(:ncol,k) + mbar(:ncol,k) * q(:ncol,k,n) / cnst_mw(n)
	    end do
         end if
         call cnst_get_ind( 'Op', n, abort=.false. )
         if (n < 0) then 
            ns = slvd_index( 'Op' )
	    call pbuf_get_field(pbuf, slvd_pbf_ndx, qs, start=(/1,1,ns/), kount=(/pcols,pver,1/) )
            nc = get_spc_ndx('Op')
      	    do k = 1,pver
	      wrk(:ncol,k) = wrk(:ncol,k) + mbar(:ncol,k) * qs(:ncol,k) / adv_mass(nc)
	    end do	   
         else
	    do k = 1,pver
	      wrk(:ncol,k) = wrk(:ncol,k) + mbar(:ncol,k) * q(:ncol,k,n) / cnst_mw(n)
	    end do
         end if
         call cnst_get_ind( 'O2p', n, abort=.false. )
         if (n < 0) then 
            ns = slvd_index( 'O2p' )
	    call pbuf_get_field(pbuf, slvd_pbf_ndx, qs, start=(/1,1,ns/), kount=(/pcols,pver,1/) )
            nc = get_spc_ndx('O2p')
      	    do k = 1,pver
	      wrk(:ncol,k) = wrk(:ncol,k) + mbar(:ncol,k) * qs(:ncol,k) / adv_mass(nc)
	    end do	   
         else
	    do k = 1,pver
	      wrk(:ncol,k) = wrk(:ncol,k) + mbar(:ncol,k) * q(:ncol,k,n) / cnst_mw(n)
	    end do
         end if
         call cnst_get_ind( 'NOp', n, abort=.false. )
         if (n < 0) then 
            ns = slvd_index( 'NOp' )
	    call pbuf_get_field(pbuf, slvd_pbf_ndx, qs, start=(/1,1,ns/), kount=(/pcols,pver,1/) )
            nc = get_spc_ndx('NOp')
      	    do k = 1,pver
	      wrk(:ncol,k) = wrk(:ncol,k) + mbar(:ncol,k) * qs(:ncol,k) / adv_mass(nc)
	    end do	   
         else
	    do k = 1,pver
	      wrk(:ncol,k) = wrk(:ncol,k) + mbar(:ncol,k) * q(:ncol,k,n) / cnst_mw(n)
	    end do
         end if
	 
         !--------------------------------------------------------------------------------------
         !  Total ions now in wrk array so determine electrons.  qse is a pointer to pbuf and
	 !  q is a pointer to state%q 
         !--------------------------------------------------------------------------------------
         if (elec_ndx < 0) then 
            ns = slvd_index( 'e' )
            call pbuf_get_field(pbuf, slvd_pbf_ndx, qse, start=(/1,1,ns/), kount=(/pcols,pver,1/) )
            nc = get_spc_ndx('e')
      	    do k = 1,pver
	      qse(:ncol,k) = adv_mass(nc) * wrk(:ncol,k) / mbar(:ncol,k)
	    end do	   
         else
           do k = 1,pver
	     q(:ncol,k,elec_ndx) = cnst_mw(elec_ndx) * wrk(:ncol,k) / mbar(:ncol,k)
	   end do
	 end if  
      end if

      end subroutine charge_fix

      end module charge_neutrality
