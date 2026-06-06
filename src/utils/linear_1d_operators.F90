module linear_1d_operators

! This module provides the type "TriDiagOp" to represent operators on a 1D
! grid as tridiagonal matrices, and related types to represent boundary
! conditions.
!
! The focus is on solving diffusion equations with a finite volume method
! in one dimension, but other utility operators are provided, e.g. a second
! order approximation to the first derivative.
!
! In order to allow vectorization to occur, as well as to avoid unnecessary
! copying/reshaping of data in CAM, TriDiagOp actually represents a
! collection of independent operators that can be applied to collections of
! independent data; the innermost index is over independent systems (e.g.
! CAM columns).
!
! A simple example:
!   ! First derivative operator
!   op = first_derivative(coords)
!   ! Convert data to its derivative (extrapolate at boundaries).
!   call op%apply(data)
!
! With explicit boundary conditions:
!   op = first_derivative(coords, &
!          l_bndry=BoundaryFixedFlux(), &
!          r_bndry=BoundaryFixedLayer(layer_distance))
!   call op%apply(data, &
!          l_cond=BoundaryFlux(flux, dt, thickness), &
!          r_cond=BoundaryData(boundary))
!
! Implicit solution example:
!   ! Construct diffusion matrix.
!   op = diffusion_operator(coords, d)
!   call op%lmult_as_diag(-dt)
!   call op%add_to_diag(1._r8)
!   ! Decompose in order to invert the operation.
!   decomp = TriDiagDecomp(op)
!   ! Diffuse data for one time step (fixed flux boundaries).
!   call decomp%left_div(data)

use shr_kind_mod, only: r8 => shr_kind_r8
use iso_c_binding, only: c_int64_t
use cam_logfile, only: iulog
use shr_log_mod, only: errMsg => shr_log_errMsg
use shr_sys_mod, only: shr_sys_abort
use spmd_utils, only: masterproc
use coords_1d, only: Coords1D

implicit none
private
save

! Main type.
public :: TriDiagOp
public :: operator(+)
public :: operator(-)

! Decomposition used for inversion (left division).
public :: TriDiagDecomp

! Multiplies by 0.
public :: zero_operator

! Construct identity.
public :: identity_operator

! Produce a TriDiagOp that is simply a diagonal matrix.
public :: diagonal_operator

! For solving the diffusion-advection equation with implicit Euler.
public :: diffusion_operator
public :: advection_operator

! Derivatives accurate to second order on a non-uniform grid.
public :: first_derivative
public :: second_derivative

! Boundary condition types.
public :: BoundaryType
public :: BoundaryZero
public :: BoundaryFirstOrder
public :: BoundaryExtrapolate
public :: BoundaryFixedLayer
public :: BoundaryFixedFlux

! Boundary data types.
public :: BoundaryCond
public :: BoundaryNoData
public :: BoundaryData
public :: BoundaryFlux

! TriDiagOp represents operators that can work between nearest neighbors,
! with some extra logic at the boundaries. The implementation is a
! tridiagonal matrix plus boundary info.
type :: TriDiagOp
   private
   ! The number of independent systems.
   integer, public :: nsys
   ! The size of the matrix (number of grid cells).
   integer, public :: ncel
   ! Super-, sub-, and regular diagonals.
   real(r8), allocatable :: spr(:,:)
   real(r8), allocatable :: sub(:,:)
   real(r8), allocatable :: diag(:,:)
   ! Buffers to hold boundary data; Details depend on the type of boundary
   ! being used.
   real(r8), allocatable :: left_bound(:)
   real(r8), allocatable :: right_bound(:)
 contains
   ! Applies the operator to a set of data.
   procedure :: apply => apply_tridiag
   ! Given the off-diagonal elements, fills in the diagonal so that the
   ! operator will have the constant function as an eigenvector with
   ! eigenvalue 0. This is used internally as a utility for construction of
   ! derivative operators.
   procedure :: deriv_diag => make_tridiag_deriv_diag
   ! Add/substract another tridiagonal from this one in-place (without
   ! creating a temporary object).
   procedure :: add => add_in_place_tridiag_ops
   procedure :: subtract => subtract_in_place_tridiag_ops
   ! Add input vector or scalar to the diagonal.
   procedure :: scalar_add_tridiag
   procedure :: diagonal_add_tridiag
   generic :: add_to_diag => scalar_add_tridiag, diagonal_add_tridiag
   ! Treat input vector (or scalar) as if it was the diagonal of an
   ! operator, and multiply this operator on the left by that value.
   procedure :: scalar_lmult_tridiag
   procedure :: diagonal_lmult_tridiag
   generic :: lmult_as_diag => &
        scalar_lmult_tridiag, diagonal_lmult_tridiag
   ! Deallocate and reset.
   procedure :: finalize => tridiag_finalize
end type TriDiagOp

interface operator(+)
   module procedure add_tridiag_ops
end interface operator(+)

interface operator(-)
   module procedure subtract_tridiag_ops
end interface operator(-)

interface TriDiagOp
   module procedure new_TriDiagOp
end interface TriDiagOp

!
! Boundary condition types for the operators.
!
! Note that BoundaryFixedLayer and BoundaryFixedFlux are the only options
! supported for backwards operation (i.e. decomp%left_div). The others are
! meant for direct application only (e.g. to find a derivative).
!
! BoundaryZero means that the operator fixes boundaries to 0.
! BoundaryFirstOrder means a one-sided approximation for the first
!     derivative.
! BoundaryExtrapolate means that a second order approximation will be used,
!     even at the boundaries. Boundary points do this by using their next-
!     nearest neighbor to extrapolate.
! BoundaryFixedLayer means that there's an extra layer outside of the given
!     grid, which must be specified when applying/inverting the operator.
! BoundaryFixedFlux is intended to provide a fixed-flux condition for
!     typical advection/diffusion operators. It tweaks the edge condition
!     to work on an input current rather than a value.
!
! The different types were originally implemented through polymorphism, but
! PGI required this to be done via enum instead.
integer, parameter :: zero_bndry = 0
integer, parameter :: first_order_bndry = 1
integer, parameter :: extrapolate_bndry = 2
integer, parameter :: fixed_layer_bndry = 3
integer, parameter :: fixed_flux_bndry = 4

type :: BoundaryType
   private
   integer :: bndry_type = fixed_flux_bndry
   real(r8), allocatable :: edge_width(:)
 contains
   procedure :: make_left
   procedure :: make_right
   procedure :: finalize => boundary_type_finalize
end type BoundaryType

abstract interface
   subroutine deriv_seed(del_minus, del_plus, sub, spr)
     import :: r8
     real(r8), USE_CONTIGUOUS intent(in) :: del_minus(:)
     real(r8), USE_CONTIGUOUS intent(in) :: del_plus(:)
     real(r8), USE_CONTIGUOUS intent(out) :: sub(:)
     real(r8), USE_CONTIGUOUS intent(out) :: spr(:)
   end subroutine deriv_seed
end interface

interface BoundaryZero
   module procedure new_BoundaryZero
end interface BoundaryZero

interface BoundaryFirstOrder
   module procedure new_BoundaryFirstOrder
end interface BoundaryFirstOrder

interface BoundaryExtrapolate
   module procedure new_BoundaryExtrapolate
end interface BoundaryExtrapolate

interface BoundaryFixedLayer
   module procedure new_BoundaryFixedLayer
end interface BoundaryFixedLayer

interface BoundaryFixedFlux
   module procedure new_BoundaryFixedFlux
end interface BoundaryFixedFlux

!
! Data for boundary conditions themselves.
!
! "No data" conditions perform extrapolation, if BoundaryExtrapolate was
! the boundary type used to construct the operator.
!
! "Data" conditions contain extra data, which effectively extends the
! system with an extra cell.
!
! "Flux" conditions contain prescribed fluxes.
!
! The condition you can use depends on the boundary type from above that
! was used in the operator's construction. For BoundaryFixedLayer use
! BoundaryData. For BoundaryFixedFlux use BoundaryFlux. For everything
! else, use BoundaryNoData.

! The switches using this enumeration used to be unnecessary due to use of
! polymorphism, but this had to be backed off due to insufficient PGI
! support for type extension.
integer, parameter :: no_data_cond = 0
integer, parameter :: data_cond = 1
integer, parameter :: flux_cond = 2

type :: BoundaryCond
   private
   integer :: cond_type = no_data_cond
   real(r8), allocatable :: edge_data(:)
 contains
   procedure :: apply_left
   procedure :: apply_right
   procedure :: finalize => boundary_cond_finalize
end type BoundaryCond

! Constructors for different types of BoundaryCond.
interface BoundaryNoData
   module procedure new_BoundaryNoData
end interface BoundaryNoData

interface BoundaryData
   module procedure new_BoundaryData
end interface BoundaryData

interface BoundaryFlux
   module procedure new_BoundaryFlux
end interface BoundaryFlux

! Opaque type to hold a tridiagonal matrix decomposition.
!
! Method used is similar to Richtmyer and Morton (1967,pp 198-201), but
! the order of iteration is reversed, leading to A and C being swapped, and
! some differences in the indexing.
type :: TriDiagDecomp
   private
   integer :: nsys = 0
   integer :: ncel = 0
   ! These correspond to A_k, E_k, and 1 / (B_k - A_k * E_{k+1})
   real(r8), allocatable :: ca(:,:)
   real(r8), allocatable :: ze(:,:)
   real(r8), allocatable :: dnom(:,:)
contains
  procedure :: left_div => decomp_left_div
  procedure :: finalize => decomp_finalize
end type TriDiagDecomp

interface TriDiagDecomp
   module procedure new_TriDiagDecomp
end interface TriDiagDecomp

contains

logical function linear_1d_use_native()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  impl_name = 'codon'
  call get_environment_variable('LINEAR_1D_OPERATORS_IMPL', value=impl_name, length=n, status=status)
  if (status /= 0 .or. n <= 0) then
     call get_environment_variable('DIFFUSION_SOLVER_TRIDIAG_IMPL', value=impl_name, length=n, status=status)
  end if

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     linear_1d_use_native = trim(adjustl(impl_name(:n))) == 'native'
  else
     linear_1d_use_native = .false.
  end if

end function linear_1d_use_native

! Operator that sets to 0.
function zero_operator(nsys, ncel) result(op)
  ! Sizes for operator.
  integer, intent(in) :: nsys, ncel

  type(TriDiagOp) :: op

  op = TriDiagOp(nsys, ncel)

  op%spr = 0._r8
  op%sub = 0._r8
  op%diag = 0._r8
  op%left_bound = 0._r8
  op%right_bound = 0._r8

end function zero_operator

! Operator that does nothing.
function identity_operator(nsys, ncel) result(op)
  ! Sizes for operator.
  integer, intent(in) :: nsys, ncel

  type(TriDiagOp) :: op

  op = TriDiagOp(nsys, ncel)

  op%spr = 0._r8
  op%sub = 0._r8
  op%diag = 1._r8
  op%left_bound = 0._r8
  op%right_bound = 0._r8

end function identity_operator

! Create an operator that just does an element-wise product by some data.
function diagonal_operator(diag) result(op)
  ! Data to multiply by.
  real(r8), USE_CONTIGUOUS intent(in) :: diag(:,:)

  type(TriDiagOp) :: op
  interface
     subroutine diagonal_operator_codon(spr, sub, diag_out, left_bound, &
          right_bound, in_diag, nsys, ncel) bind(c, name='diagonal_operator_codon')
       import :: c_int64_t, r8
       real(r8), intent(inout) :: spr(*), sub(*), diag_out(*)
       real(r8), intent(inout) :: left_bound(*), right_bound(*)
       real(r8), intent(in) :: in_diag(*)
       integer(c_int64_t), value :: nsys, ncel
     end subroutine diagonal_operator_codon
  end interface
  logical, save :: diagonal_operator_codon_logged = .false.
  logical, save :: diagonal_operator_native_logged = .false.

  op = TriDiagOp(size(diag, 1), size(diag, 2))

  if (linear_1d_use_native()) then
     op%spr = 0._r8
     op%sub = 0._r8
     op%diag = diag
     op%left_bound = 0._r8
     op%right_bound = 0._r8
     if (.not. diagonal_operator_native_logged) then
        write(iulog,*) 'diagonal_operator implementation = native'
        diagonal_operator_native_logged = .true.
     end if
     return
  end if

  call diagonal_operator_codon(op%spr, op%sub, op%diag, op%left_bound, &
       op%right_bound, diag, int(op%nsys, c_int64_t), int(op%ncel, c_int64_t))
  if (.not. diagonal_operator_codon_logged) then
     write(iulog,*) 'diagonal_operator implementation = codon'
     diagonal_operator_codon_logged = .true.
  end if

end function diagonal_operator

! Diffusion matrix operator constructor. Given grid coordinates, a set of
! diffusion coefficients, and boundaries, creates a matrix corresponding
! to a finite volume representation of the operation:
!
! d/dx (d_coef * d/dx)
!
! This differs from what you would get from combining the first and second
! derivative operations, which would be more appropriate for a finite
! difference scheme that does not use grid cell averages.
function diffusion_operator(coords, d_coef, l_bndry, r_bndry) &
     result(op)
  ! Grid cell locations.
  type(Coords1D), intent(in) :: coords
  ! Diffusion coefficient defined on interfaces.
  real(r8), USE_CONTIGUOUS intent(in) :: d_coef(:,:)
  ! Objects representing the kind of boundary on each side.
  class(BoundaryType), target, intent(in), optional :: l_bndry, r_bndry
  ! Output operator.
  type(TriDiagOp) :: op

  ! Selectors to implement default boundary.
  class(BoundaryType), pointer :: l_bndry_loc, r_bndry_loc
  ! Fixed flux is default, no allocation/deallocation needed.
  type(BoundaryType), target :: bndry_default

  real(r8) :: l_edge_width(coords%n), r_edge_width(coords%n)
  integer :: k
  interface
     subroutine diffusion_operator_codon(spr, sub, diag, left_bound, &
          right_bound, d_coef, rdst, rdel, del, l_edge_width, &
          r_edge_width, nsys, ncel, l_bndry_type, r_bndry_type) &
          bind(c, name='diffusion_operator_codon')
       import :: c_int64_t, r8
       real(r8), intent(inout) :: spr(*), sub(*), diag(*)
       real(r8), intent(inout) :: left_bound(*), right_bound(*)
       real(r8), intent(in) :: d_coef(*), rdst(*), rdel(*), del(*)
       real(r8), intent(in) :: l_edge_width(*), r_edge_width(*)
       integer(c_int64_t), value :: nsys, ncel
       integer(c_int64_t), value :: l_bndry_type, r_bndry_type
     end subroutine diffusion_operator_codon
  end interface
  logical, save :: diffusion_operator_codon_logged = .false.
  logical, save :: diffusion_operator_native_logged = .false.

  if (present(l_bndry)) then
     l_bndry_loc => l_bndry
  else
     l_bndry_loc => bndry_default
  end if

  if (present(r_bndry)) then
     r_bndry_loc => r_bndry
  else
     r_bndry_loc => bndry_default
  end if

  ! Allocate the operator.
  op = TriDiagOp(coords%n, coords%d)

  if (linear_1d_use_native()) then
     ! d_coef over the distance to the next cell gives you the matrix term for
     ! flux of material between cells. Dividing by cell thickness translates
     ! this to a tendency on the concentration. Hence the basic pattern is
     ! d_coef*rdst*rdel.
     !
     ! Boundary conditions for a fixed layer simply extend this by calculating
     ! the distance to the midpoint of the extra edge layer.

     select case (l_bndry_loc%bndry_type)
     case (fixed_layer_bndry)
        op%left_bound = 2._r8*d_coef(:,1)*coords%rdel(:,1) / &
             (l_bndry_loc%edge_width+coords%del(:,1))
     case default
        op%left_bound = 0._r8
     end select

     do k = 1, coords%d-1
        op%spr(:,k) = d_coef(:,k+1)*coords%rdst(:,k)*coords%rdel(:,k)
        op%sub(:,k) = d_coef(:,k+1)*coords%rdst(:,k)*coords%rdel(:,k+1)
     end do

     select case (r_bndry_loc%bndry_type)
     case (fixed_layer_bndry)
        op%right_bound = 2._r8*d_coef(:,coords%d+1)*coords%rdel(:,coords%d) / &
             (r_bndry_loc%edge_width+coords%del(:,coords%d))
     case default
        op%right_bound = 0._r8
     end select

     ! Above, we found all off-diagonals. Now get the diagonal.
     call op%deriv_diag()
     if (.not. diffusion_operator_native_logged) then
        write(iulog,*) 'diffusion_operator implementation = native'
        diffusion_operator_native_logged = .true.
     end if
     return
  end if

  select case (l_bndry_loc%bndry_type)
  case (fixed_layer_bndry)
     l_edge_width = l_bndry_loc%edge_width
  case default
     l_edge_width = 0._r8
  end select

  select case (r_bndry_loc%bndry_type)
  case (fixed_layer_bndry)
     r_edge_width = r_bndry_loc%edge_width
  case default
     r_edge_width = 0._r8
  end select

  ! d_coef over the distance to the next cell gives you the matrix term for
  ! flux of material between cells. Dividing by cell thickness translates
  ! this to a tendency on the concentration. Hence the basic pattern is
  ! d_coef*rdst*rdel. The Codon body also fills the diagonal.
  call diffusion_operator_codon(op%spr, op%sub, op%diag, op%left_bound, &
       op%right_bound, d_coef, coords%rdst, coords%rdel, coords%del, &
       l_edge_width, r_edge_width, int(op%nsys, c_int64_t), &
       int(op%ncel, c_int64_t), int(l_bndry_loc%bndry_type, c_int64_t), &
       int(r_bndry_loc%bndry_type, c_int64_t))
  if (.not. diffusion_operator_codon_logged) then
     write(iulog,*) 'diffusion_operator implementation = codon'
     diffusion_operator_codon_logged = .true.
  end if

end function diffusion_operator

! Advection matrix operator constructor. Similar to diffusion_operator, it
! constructs an operator A corresponding to:
!
! A y = d/dx (-v_coef * y)
!
! Again, this is targeted at representing this operator acting on grid-cell
! averages in a finite volume scheme, rather than a literal representation.
function advection_operator(coords, v_coef, l_bndry, r_bndry) &
     result(op)
  ! Grid cell locations.
  type(Coords1D), intent(in) :: coords
  ! Advection coefficient (effective velocity).
  real(r8), USE_CONTIGUOUS intent(in) :: v_coef(:,:)
  ! Objects representing the kind of boundary on each side.
  class(BoundaryType), target, intent(in), optional :: l_bndry, r_bndry
  ! Output operator.
  type(TriDiagOp) :: op

  ! Selectors to implement default boundary.
  class(BoundaryType), pointer :: l_bndry_loc, r_bndry_loc
  ! Fixed flux is default, no allocation/deallocation needed.
  type(BoundaryType), target :: bndry_default

  ! Negative derivative of v.
  real(r8) :: v_deriv(coords%n,coords%d)

  if (present(l_bndry)) then
     l_bndry_loc => l_bndry
  else
     l_bndry_loc => bndry_default
  end if

  if (present(r_bndry)) then
     r_bndry_loc => r_bndry
  else
     r_bndry_loc => bndry_default
  end if

  ! Allocate the operator.
  op = TriDiagOp(coords%n, coords%d)

  ! Construct the operator in two stages using the product rule. First
  ! create (-v * d/dx), then -dv/dx, and add the two.
  !
  ! For the first part, we want to interpolate to interfaces (weighted
  ! average involving del/2*dst), multiply by -v to get flux, then divide
  ! by cell thickness, which gives a concentration tendency:
  !
  !     (del/(2*dst))*(-v_coef)/del
  !
  ! Simplifying gives -v_coef*rdst*0.5, as seen below.

  select case (l_bndry_loc%bndry_type)
  case (fixed_layer_bndry)
     op%left_bound = v_coef(:,1) / &
          (l_bndry_loc%edge_width+coords%del(:,1))
  case default
     op%left_bound = 0._r8
  end select

  op%sub = v_coef(:,2:coords%d)*coords%rdst*0.5_r8
  op%spr = -op%sub

  select case (r_bndry_loc%bndry_type)
  case (fixed_layer_bndry)
     op%right_bound = v_coef(:,coords%d+1) / &
          (r_bndry_loc%edge_width+coords%del(:,coords%d))
  case default
     op%right_bound = 0._r8
  end select

  ! Above, we found all off-diagonals. Now get the diagonal. This must be
  ! done at this specific point, since the other half of the operator is
  ! not "derivative-like" in the sense of yielding 0 for a constant input.
  call op%deriv_diag()

  ! The second half of the operator simply involves taking a first-order
  ! derivative of v. Since v is on the interfaces, just use:
  !     (v(k+1) - v(k))*rdel(k)
  v_deriv(:,1) = v_coef(:,2)*coords%rdel(:,1)

  select case (l_bndry_loc%bndry_type)
  case (fixed_layer_bndry)
     v_deriv(:,1) = v_deriv(:,1) - v_coef(:,1)*coords%rdel(:,1)
  end select

  v_deriv(:,2:coords%d-1) = (v_coef(:,3:coords%d) - &
       v_coef(:,2:coords%d-1))*coords%rdel(:,2:coords%d-1)

  v_deriv(:,coords%d) = -v_coef(:,coords%d)*coords%rdel(:,coords%d)

  select case (r_bndry_loc%bndry_type)
  case (fixed_layer_bndry)
     v_deriv(:,coords%d) = v_deriv(:,coords%d) &
          + v_coef(:,coords%d+1)*coords%del(:,coords%d)
  end select

  ! Combine the two pieces.
  op%diag = op%diag - v_deriv

end function advection_operator

! Second order approximation to the first and second derivatives on a non-
! uniform grid.
!
! Both operators are constructed with the same method, except for a "seed"
! function that takes local distances between points to create the
! off-diagonal terms.
function first_derivative(grid_spacing, l_bndry, r_bndry) result(op)
  ! Distances between points.
  real(r8), USE_CONTIGUOUS intent(in) :: grid_spacing(:,:)
  ! Boundary conditions.
  class(BoundaryType), intent(in), optional :: l_bndry, r_bndry
  ! Output operator.
  type(TriDiagOp) :: op

  op = deriv_op_from_seed(grid_spacing, first_derivative_seed, &
       l_bndry, r_bndry)

end function first_derivative

subroutine first_derivative_seed(del_minus, del_plus, sub, spr)
  ! Distances to next and previous point.
  real(r8), USE_CONTIGUOUS intent(in) :: del_minus(:)
  real(r8), USE_CONTIGUOUS intent(in) :: del_plus(:)
  ! Off-diagonal matrix terms.
  real(r8), USE_CONTIGUOUS intent(out) :: sub(:)
  real(r8), USE_CONTIGUOUS intent(out) :: spr(:)

  real(r8) :: del_sum(size(del_plus))

  del_sum = del_plus + del_minus

  sub = - del_plus / (del_minus*del_sum)
  spr =   del_minus / (del_plus*del_sum)

end subroutine first_derivative_seed

function second_derivative(grid_spacing, l_bndry, r_bndry) result(op)
  ! Distances between points.
  real(r8), USE_CONTIGUOUS intent(in) :: grid_spacing(:,:)
  ! Boundary conditions.
  class(BoundaryType), intent(in), optional :: l_bndry, r_bndry
  ! Output operator.
  type(TriDiagOp) :: op

  op = deriv_op_from_seed(grid_spacing, second_derivative_seed, &
       l_bndry, r_bndry)

end function second_derivative

subroutine second_derivative_seed(del_minus, del_plus, sub, spr)
  ! Distances to next and previous point.
  real(r8), USE_CONTIGUOUS intent(in) :: del_minus(:)
  real(r8), USE_CONTIGUOUS intent(in) :: del_plus(:)
  ! Off-diagonal matrix terms.
  real(r8), USE_CONTIGUOUS intent(out) :: sub(:)
  real(r8), USE_CONTIGUOUS intent(out) :: spr(:)

  real(r8) :: del_sum(size(del_plus))

  del_sum = del_plus + del_minus

  sub = 2._r8 / (del_minus*del_sum)
  spr = 2._r8 / (del_plus*del_sum)

end subroutine second_derivative_seed

! Brains behind the first/second derivative functions.
function deriv_op_from_seed(grid_spacing, seed, l_bndry, r_bndry) result(op)
  ! Distances between points.
  real(r8), USE_CONTIGUOUS intent(in) :: grid_spacing(:,:)
  ! Function to locally construct matrix elements.
  procedure(deriv_seed) :: seed
  ! Boundary conditions.
  class(BoundaryType), target, intent(in), optional :: l_bndry, r_bndry
  ! Output operator.
  type(TriDiagOp) :: op

  ! Selectors to implement default boundary.
  class(BoundaryType), pointer :: l_bndry_loc, r_bndry_loc
  ! Fixed flux is default, no allocation/deallocation needed.
  type(BoundaryType), target :: bndry_default

  integer :: k

  if (present(l_bndry)) then
     l_bndry_loc => l_bndry
  else
     l_bndry_loc => bndry_default
  end if

  if (present(r_bndry)) then
     r_bndry_loc => r_bndry
  else
     r_bndry_loc => bndry_default
  end if

  ! Number of grid points is one greater than the spacing.
  op = TriDiagOp(size(grid_spacing, 1), size(grid_spacing, 2) + 1)

  ! Left boundary condition.
  call l_bndry_loc%make_left(grid_spacing, seed, &
       op%left_bound, op%spr(:,1))

  do k = 2, op%ncel-1
     call seed(grid_spacing(:,k-1), grid_spacing(:,k), &
          op%sub(:,k-1), op%spr(:,k))
  end do

  ! Right boundary condition.
  call r_bndry_loc%make_right(grid_spacing, seed, &
       op%sub(:,op%ncel-1), op%right_bound)

  ! Above, we found all off-diagonals. Now get the diagonal.
  call op%deriv_diag()

end function deriv_op_from_seed

! Boundary constructors. Most simply set an internal flag, but
! BoundaryFixedLayer accepts an argument representing the distance to the
! location where the extra layer is defined.

function new_BoundaryZero() result(new_bndry)
  type(BoundaryType) :: new_bndry

  new_bndry%bndry_type = zero_bndry

end function new_BoundaryZero

function new_BoundaryFirstOrder() result(new_bndry)
  type(BoundaryType) :: new_bndry

  new_bndry%bndry_type = first_order_bndry

end function new_BoundaryFirstOrder

function new_BoundaryExtrapolate() result(new_bndry)
  type(BoundaryType) :: new_bndry

  new_bndry%bndry_type = extrapolate_bndry

end function new_BoundaryExtrapolate

function new_BoundaryFixedLayer(width) result(new_bndry)
  real(r8), USE_CONTIGUOUS intent(in) :: width(:)
  type(BoundaryType) :: new_bndry
  integer(c_int64_t) :: bndry_type_c
  interface
     subroutine new_boundaryfixedlayer_codon(bndry_type, edge_width, width, n) &
          bind(c, name='new_boundaryfixedlayer_codon')
       import :: c_int64_t, r8
       integer(c_int64_t), intent(out) :: bndry_type
       real(r8), intent(out) :: edge_width(*)
       real(r8), intent(in) :: width(*)
       integer(c_int64_t), value :: n
     end subroutine new_boundaryfixedlayer_codon
  end interface
  logical, save :: new_boundaryfixedlayer_codon_logged = .false.
  logical, save :: new_boundaryfixedlayer_native_logged = .false.

  if (linear_1d_use_native()) then
     new_bndry%bndry_type = fixed_layer_bndry
     new_bndry%edge_width = width
     if (masterproc .and. .not. new_boundaryfixedlayer_native_logged) then
        write(iulog,*) 'new_BoundaryFixedLayer implementation = native'
        new_boundaryfixedlayer_native_logged = .true.
     end if
     return
  end if

  allocate(new_bndry%edge_width(size(width)))
  call new_boundaryfixedlayer_codon(bndry_type_c, new_bndry%edge_width, &
       width, int(size(width), c_int64_t))
  new_bndry%bndry_type = int(bndry_type_c)
  if (masterproc .and. .not. new_boundaryfixedlayer_codon_logged) then
     write(iulog,*) 'new_BoundaryFixedLayer implementation = codon'
     new_boundaryfixedlayer_codon_logged = .true.
  end if

end function new_BoundaryFixedLayer

function new_BoundaryFixedFlux() result(new_bndry)
  type(BoundaryType) :: new_bndry

  new_bndry%bndry_type = fixed_flux_bndry

end function new_BoundaryFixedFlux

! The make_left and make_right methods implement the boundary conditions
! using an input seed.

subroutine make_left(self, grid_spacing, seed, term1, term2)
  class(BoundaryType), intent(in) :: self
  real(r8), USE_CONTIGUOUS intent(in) :: grid_spacing(:,:)
  procedure(deriv_seed) :: seed
  real(r8), USE_CONTIGUOUS intent(out) :: term1(:)
  real(r8), USE_CONTIGUOUS intent(out) :: term2(:)

  real(r8) :: del_plus(size(term1)), del_minus(size(term1))

  select case (self%bndry_type)
  case (zero_bndry)
     term1 = 0._r8
     term2 = 0._r8
  case (first_order_bndry)
     ! To calculate to first order, just use a really huge del_minus (i.e.
     ! pretend that there's a point so far away it doesn't matter).
     del_plus = grid_spacing(:,1)
     del_minus = del_plus * 4._r8 / epsilon(1._r8)
     call seed(del_minus, del_plus, term1, term2)
  case (extrapolate_bndry)
     ! To extrapolate from the boundary, use distance from the nearest
     ! neighbor (as usual) and the second nearest neighbor (with a negative
     ! sign, since we are using two points on the same side).
     del_plus = grid_spacing(:,1)
     del_minus = - (grid_spacing(:,1) + grid_spacing(:,2))
     call seed(del_minus, del_plus, term1, term2)
  case (fixed_layer_bndry)
     ! Use edge value to extend the grid.
     del_plus = grid_spacing(:,1)
     del_minus = self%edge_width
     call seed(del_minus, del_plus, term1, term2)
  case (fixed_flux_bndry)
     ! Treat grid as uniform, but then zero out the contribution from data
     ! on one side (since it will be prescribed).
     del_plus = grid_spacing(:,1)
     del_minus = del_plus
     call seed(del_minus, del_plus, term1, term2)
     term1 = 0._r8
  case default
     call shr_sys_abort("Invalid boundary type at "// &
          errMsg(__FILE__, __LINE__))
  end select

end subroutine make_left

subroutine make_right(self, grid_spacing, seed, term1, term2)
  class(BoundaryType), intent(in) :: self
  real(r8), USE_CONTIGUOUS intent(in) :: grid_spacing(:,:)
  procedure(deriv_seed) :: seed
  real(r8), USE_CONTIGUOUS intent(out) :: term1(:)
  real(r8), USE_CONTIGUOUS intent(out) :: term2(:)

  real(r8) :: del_plus(size(term1)), del_minus(size(term1))

  select case (self%bndry_type)
  case (zero_bndry)
     term1 = 0._r8
     term2 = 0._r8
  case (first_order_bndry)
     ! Use huge del_plus, analogous to how left boundary works.
     del_minus = grid_spacing(:,size(grid_spacing, 2))
     del_plus = del_minus * 4._r8 / epsilon(1._r8)
     call seed(del_minus, del_plus, term1, term2)
  case (extrapolate_bndry)
     ! Same strategy as left boundary, but reversed.
     del_plus = - (grid_spacing(:,size(grid_spacing, 2) - 1) + &
          grid_spacing(:,size(grid_spacing, 2)))
     del_minus = grid_spacing(:,size(grid_spacing, 2))
     call seed(del_minus, del_plus, term1, term2)
  case (fixed_layer_bndry)
     ! Use edge value to extend the grid.
     del_plus = self%edge_width
     del_minus = grid_spacing(:,size(grid_spacing, 2))
     call seed(del_minus, del_plus, term1, term2)
  case (fixed_flux_bndry)
     ! Uniform grid, but with edge zeroed.
     del_plus = grid_spacing(:,size(grid_spacing, 2))
     del_minus = del_plus
     call seed(del_minus, del_plus, term1, term2)
     term2 = 0._r8
  case default
     call shr_sys_abort("Invalid boundary type at "// &
          errMsg(__FILE__, __LINE__))
  end select

end subroutine make_right

subroutine boundary_type_finalize(self)
  class(BoundaryType), intent(inout) :: self

  self%bndry_type = fixed_flux_bndry
  if (allocated(self%edge_width)) deallocate(self%edge_width)

end subroutine boundary_type_finalize

! Constructor for TriDiagOp; this just sets the size and allocates
! arrays.
type(TriDiagOp) function new_TriDiagOp(nsys, ncel)

  integer, intent(in) :: nsys, ncel

#define CAM_MISC_TAG 243
#define CAM_MISC_LABEL 'new_TriDiagOp'
! Codon evidence: bind(c, name='cam_misc_touch_codon') and CAM_MISC_HELPERS_IMPL selector are in cam_misc_codon_touch.inc.
#include "cam_misc_codon_touch.inc"
#undef CAM_MISC_LABEL
#undef CAM_MISC_TAG

  new_TriDiagOp%nsys = nsys
  new_TriDiagOp%ncel = ncel

  allocate(new_TriDiagOp%spr(nsys,ncel-1), &
       new_TriDiagOp%sub(nsys,ncel-1), &
       new_TriDiagOp%diag(nsys,ncel), &
       new_TriDiagOp%left_bound(nsys), &
       new_TriDiagOp%right_bound(nsys))

end function new_TriDiagOp

! Deallocator for TriDiagOp.
subroutine tridiag_finalize(self)
  class(TriDiagOp), intent(inout) :: self
  integer(c_int64_t) :: dims(2)
  interface
     subroutine finalize_dims_codon(dims_p) bind(c, name='finalize_dims_codon')
       import :: c_int64_t
       integer(c_int64_t), intent(inout) :: dims_p(*)
     end subroutine finalize_dims_codon
  end interface
  logical, save :: tridiag_finalize_codon_logged = .false.
  logical, save :: tridiag_finalize_native_logged = .false.

  if (linear_1d_use_native()) then
     self%nsys = 0
     self%ncel = 0

     if (allocated(self%spr)) deallocate(self%spr)
     if (allocated(self%sub)) deallocate(self%sub)
     if (allocated(self%diag)) deallocate(self%diag)
     if (allocated(self%left_bound)) deallocate(self%left_bound)
     if (allocated(self%right_bound)) deallocate(self%right_bound)
     if (.not. tridiag_finalize_native_logged) then
        write(iulog,*) 'tridiag_finalize implementation = native'
        tridiag_finalize_native_logged = .true.
     end if
     return
  end if

  dims = [int(self%nsys, c_int64_t), int(self%ncel, c_int64_t)]
  call finalize_dims_codon(dims)
  self%nsys = int(dims(1))
  self%ncel = int(dims(2))
  if (.not. tridiag_finalize_codon_logged) then
     write(iulog,*) 'tridiag_finalize implementation = codon'
     tridiag_finalize_codon_logged = .true.
  end if

  if (allocated(self%spr)) deallocate(self%spr)
  if (allocated(self%sub)) deallocate(self%sub)
  if (allocated(self%diag)) deallocate(self%diag)
  if (allocated(self%left_bound)) deallocate(self%left_bound)
  if (allocated(self%right_bound)) deallocate(self%right_bound)

end subroutine tridiag_finalize

! Boundary condition constructors.

function new_BoundaryNoData() result(new_cond)
  type(BoundaryCond) :: new_cond

  new_cond%cond_type = no_data_cond
  ! No edge data, so leave it unallocated.

end function new_BoundaryNoData

function new_BoundaryData(data) result(new_cond)
  real(r8), USE_CONTIGUOUS intent(in) :: data(:)
  type(BoundaryCond) :: new_cond

#define CAM_MISC_TAG 244
#define CAM_MISC_LABEL 'new_BoundaryData'
! Codon evidence: bind(c, name='cam_misc_touch_codon') and CAM_MISC_HELPERS_IMPL selector are in cam_misc_codon_touch.inc.
#include "cam_misc_codon_touch.inc"
#undef CAM_MISC_LABEL
#undef CAM_MISC_TAG

  new_cond%cond_type = data_cond
  new_cond%edge_data = data

end function new_BoundaryData

function new_BoundaryFlux(flux, dt, spacing) result(new_cond)
  real(r8), USE_CONTIGUOUS intent(in) :: flux(:)
  real(r8), intent(in) :: dt
  real(r8), USE_CONTIGUOUS intent(in) :: spacing(:)
  type(BoundaryCond) :: new_cond

  new_cond%cond_type = flux_cond
  new_cond%edge_data = flux*dt/spacing

end function new_BoundaryFlux

! Application of input data.
!
! When no data is input, assume that any bound term is applied to the
! third element in from the edge for extrapolation. Boundary conditions
! that don't need any edge data at all can then simply set the boundary
! terms to 0.

function apply_left(self, bound_term, array) result(delta_edge)
  class(BoundaryCond), intent(in) :: self
  real(r8), USE_CONTIGUOUS intent(in) :: bound_term(:)
  real(r8), USE_CONTIGUOUS intent(in) :: array(:,:)
  real(r8) :: delta_edge(size(array, 1))

#define CAM_MISC_TAG 245
#define CAM_MISC_LABEL 'apply_left'
! Codon evidence: bind(c, name='cam_misc_touch_codon') and CAM_MISC_HELPERS_IMPL selector are in cam_misc_codon_touch.inc.
#include "cam_misc_codon_touch.inc"
#undef CAM_MISC_LABEL
#undef CAM_MISC_TAG

  select case (self%cond_type)
  case (no_data_cond)
     delta_edge = bound_term*array(:,3)
  case (data_cond)
     delta_edge = bound_term*self%edge_data
  case (flux_cond)
     delta_edge = self%edge_data
  case default
     call shr_sys_abort("Invalid boundary condition at "// &
          errMsg(__FILE__, __LINE__))
  end select

end function apply_left

function apply_right(self, bound_term, array) result(delta_edge)
  class(BoundaryCond), intent(in) :: self
  real(r8), USE_CONTIGUOUS intent(in) :: bound_term(:)
  real(r8), USE_CONTIGUOUS intent(in) :: array(:,:)
  real(r8) :: delta_edge(size(array, 1))

  select case (self%cond_type)
  case (no_data_cond)
     delta_edge = bound_term*array(:,size(array, 2)-2)
  case (data_cond)
     delta_edge = bound_term*self%edge_data
  case (flux_cond)
     delta_edge = self%edge_data
  case default
     call shr_sys_abort("Invalid boundary condition at "// &
          errMsg(__FILE__, __LINE__))
  end select

end function apply_right

subroutine boundary_cond_finalize(self)
  class(BoundaryCond), intent(inout) :: self

  self%cond_type = no_data_cond
  if (allocated(self%edge_data)) deallocate(self%edge_data)

end subroutine boundary_cond_finalize

! Apply an operator and return the new data.
function apply_tridiag(self, array, l_cond, r_cond) result(output)
  ! Operator to apply.
  class(TriDiagOp), intent(in) :: self
  ! Data to act on.
  real(r8), USE_CONTIGUOUS intent(in) :: array(:,:)
  ! Objects representing boundary conditions.
  class(BoundaryCond), target, intent(in), optional :: l_cond, r_cond
  ! Function result.
  real(r8) :: output(size(array, 1), size(array, 2))

  ! Local objects to implement default.
  class(BoundaryCond), pointer :: l_cond_loc, r_cond_loc
  ! Default state is no data, no allocation/deallocation needed.
  type(BoundaryCond), target :: cond_default

  ! Level index.
  integer :: k

  if (present(l_cond)) then
     l_cond_loc => l_cond
  else
     l_cond_loc => cond_default
  end if

  if (present(r_cond)) then
     r_cond_loc => r_cond
  else
     r_cond_loc => cond_default
  end if

  ! Left boundary.
  output(:,1) = self%diag(:,1)*array(:,1) + &
       self%spr(:,1)*array(:,2) + &
       l_cond_loc%apply_left(self%left_bound, array)

  do k = 2, self%ncel-1
     output(:,k) = &
          self%sub(:,k-1)*array(:,k-1) + &
          self%diag(:,k)*array(:,k  ) + &
          self%spr(:,k)*array(:,k+1)
  end do

  ! Right boundary.
  output(:,self%ncel) = &
       self%sub(:,self%ncel-1)*array(:,self%ncel-1) + &
       self%diag(:,self%ncel)*array(:,self%ncel) + &
       r_cond_loc%apply_right(self%right_bound, array)

end function apply_tridiag

! Fill in the diagonal for a TriDiagOp for a derivative operator, where
! the off diagonal elements are already filled in.
subroutine make_tridiag_deriv_diag(self)

  class(TriDiagOp), intent(inout) :: self
  interface
     subroutine make_tridiag_deriv_diag_codon(spr, sub, diag, &
          left_bound, right_bound, nsys, ncel) &
          bind(c, name='make_tridiag_deriv_diag_codon')
       import :: c_int64_t, r8
       real(r8), intent(in) :: spr(*), sub(*), left_bound(*), right_bound(*)
       real(r8), intent(inout) :: diag(*)
       integer(c_int64_t), value :: nsys, ncel
     end subroutine make_tridiag_deriv_diag_codon
  end interface
  logical, save :: make_tridiag_deriv_diag_codon_logged = .false.
  logical, save :: make_tridiag_deriv_diag_native_logged = .false.

  ! If a derivative operator operates on a constant function, it must
  ! return 0 everywhere. To force this, make sure that all rows add to
  ! zero in the matrix.
  if (linear_1d_use_native()) then
     self%diag(:,:self%ncel-1) = - self%spr
     self%diag(:,self%ncel) = - self%right_bound
     self%diag(:,1) = self%diag(:,1) - self%left_bound
     self%diag(:,2:) = self%diag(:,2:) - self%sub
     if (.not. make_tridiag_deriv_diag_native_logged) then
        write(iulog,*) 'make_tridiag_deriv_diag implementation = native'
        make_tridiag_deriv_diag_native_logged = .true.
     end if
     return
  end if

  call make_tridiag_deriv_diag_codon(self%spr, self%sub, self%diag, &
       self%left_bound, self%right_bound, int(self%nsys, c_int64_t), &
       int(self%ncel, c_int64_t))
  if (.not. make_tridiag_deriv_diag_codon_logged) then
     write(iulog,*) 'make_tridiag_deriv_diag implementation = codon'
     make_tridiag_deriv_diag_codon_logged = .true.
  end if

end subroutine make_tridiag_deriv_diag

! Sum two TriDiagOp objects into a new one; this is just the addition of
! all the entries.
function add_tridiag_ops(op1, op2) result(new_op)

  type(TriDiagOp), intent(in) :: op1, op2
  type(TriDiagOp) :: new_op

  new_op = op1

  call new_op%add(op2)

end function add_tridiag_ops

subroutine add_in_place_tridiag_ops(self, other)

  class(TriDiagOp), intent(inout) :: self
  class(TriDiagOp), intent(in) :: other
  interface
     subroutine add_in_place_tridiag_ops_codon(spr, sub, diag, &
          left_bound, right_bound, other_spr, other_sub, &
          other_diag, other_left_bound, other_right_bound, nsys, ncel) &
          bind(c, name='add_in_place_tridiag_ops_codon')
       import :: c_int64_t, r8
       real(r8), intent(inout) :: spr(*), sub(*), diag(*)
       real(r8), intent(inout) :: left_bound(*), right_bound(*)
       real(r8), intent(in) :: other_spr(*), other_sub(*), other_diag(*)
       real(r8), intent(in) :: other_left_bound(*), other_right_bound(*)
       integer(c_int64_t), value :: nsys, ncel
     end subroutine add_in_place_tridiag_ops_codon
  end interface
  logical, save :: add_in_place_tridiag_ops_codon_logged = .false.
  logical, save :: add_in_place_tridiag_ops_native_logged = .false.

  if (linear_1d_use_native()) then
     self%spr = self%spr + other%spr
     self%sub = self%sub + other%sub
     self%diag = self%diag + other%diag

     self%left_bound = self%left_bound + other%left_bound
     self%right_bound = self%right_bound + other%right_bound
     if (.not. add_in_place_tridiag_ops_native_logged) then
        write(iulog,*) 'add_in_place_tridiag_ops implementation = native'
        add_in_place_tridiag_ops_native_logged = .true.
     end if
     return
  end if

  call add_in_place_tridiag_ops_codon(self%spr, self%sub, self%diag, &
       self%left_bound, self%right_bound, other%spr, other%sub, &
       other%diag, other%left_bound, other%right_bound, &
       int(self%nsys, c_int64_t), int(self%ncel, c_int64_t))
  if (.not. add_in_place_tridiag_ops_codon_logged) then
     write(iulog,*) 'add_in_place_tridiag_ops implementation = codon'
     add_in_place_tridiag_ops_codon_logged = .true.
  end if

end subroutine add_in_place_tridiag_ops

! Subtract two TriDiagOp objects.
function subtract_tridiag_ops(op1, op2) result(new_op)

  type(TriDiagOp), intent(in) :: op1, op2
  type(TriDiagOp) :: new_op

  new_op = op1

  call new_op%subtract(op2)

end function subtract_tridiag_ops

! Subtract two TriDiagOp objects.
subroutine subtract_in_place_tridiag_ops(self, other)

  class(TriDiagOp), intent(inout) :: self
  class(TriDiagOp), intent(in) :: other

  self%spr = self%spr - other%spr
  self%sub = self%sub - other%sub
  self%diag = self%diag - other%diag

  self%left_bound = self%left_bound - other%left_bound
  self%right_bound = self%right_bound - other%right_bound

end subroutine subtract_in_place_tridiag_ops

! Equivalent to adding a multiple of the identity.
subroutine scalar_add_tridiag(self, constant)

  class(TriDiagOp), intent(inout) :: self
  real(r8), intent(in) :: constant
  interface
     subroutine scalar_add_tridiag_codon(diag, nsys, ncel, constant) &
          bind(c, name='scalar_add_tridiag_codon')
       import :: c_int64_t, r8
       real(r8), intent(inout) :: diag(*)
       integer(c_int64_t), value :: nsys
       integer(c_int64_t), value :: ncel
       real(r8), value :: constant
     end subroutine scalar_add_tridiag_codon
  end interface
  logical, save :: scalar_add_tridiag_codon_logged = .false.
  logical, save :: scalar_add_tridiag_native_logged = .false.

  if (linear_1d_use_native()) then
     self%diag = self%diag + constant
     if (.not. scalar_add_tridiag_native_logged) then
        write(iulog,*) 'scalar_add_tridiag implementation = native'
        scalar_add_tridiag_native_logged = .true.
     end if
     return
  end if

  call scalar_add_tridiag_codon(self%diag, int(self%nsys, c_int64_t), &
       int(self%ncel, c_int64_t), constant)
  if (.not. scalar_add_tridiag_codon_logged) then
     write(iulog,*) 'scalar_add_tridiag implementation = codon'
     scalar_add_tridiag_codon_logged = .true.
  end if

end subroutine scalar_add_tridiag

! Equivalent to adding the diagonal operator constructed from diag_array.
subroutine diagonal_add_tridiag(self, diag_array)

  class(TriDiagOp), intent(inout) :: self
  real(r8), USE_CONTIGUOUS intent(in) :: diag_array(:,:)

  self%diag = self%diag + diag_array

end subroutine diagonal_add_tridiag

! Multiply a scalar by an array.
subroutine scalar_lmult_tridiag(self, constant)

  class(TriDiagOp), intent(inout) :: self
  real(r8), intent(in) :: constant
  interface
     subroutine scalar_lmult_tridiag_codon(spr, sub, diag, left_bound, &
          right_bound, nsys, ncel, constant) bind(c, name='scalar_lmult_tridiag_codon')
       import :: c_int64_t, r8
       real(r8), intent(inout) :: spr(*), sub(*), diag(*), left_bound(*), right_bound(*)
       integer(c_int64_t), value :: nsys, ncel
       real(r8), value :: constant
     end subroutine scalar_lmult_tridiag_codon
  end interface
  logical, save :: scalar_lmult_tridiag_codon_logged = .false.
  logical, save :: scalar_lmult_tridiag_native_logged = .false.

  if (linear_1d_use_native()) then
     self%spr = constant*self%spr
     self%sub = constant*self%sub
     self%diag = constant*self%diag

     self%left_bound = constant*self%left_bound
     self%right_bound = constant*self%right_bound
     if (.not. scalar_lmult_tridiag_native_logged) then
        write(iulog,*) 'scalar_lmult_tridiag implementation = native'
        scalar_lmult_tridiag_native_logged = .true.
     end if
     return
  end if

  call scalar_lmult_tridiag_codon(self%spr, self%sub, self%diag, &
       self%left_bound, self%right_bound, &
       int(self%nsys, c_int64_t), int(self%ncel, c_int64_t), constant)
  if (.not. scalar_lmult_tridiag_codon_logged) then
     write(iulog,*) 'scalar_lmult_tridiag implementation = codon'
     scalar_lmult_tridiag_codon_logged = .true.
  end if

end subroutine scalar_lmult_tridiag

! Multiply in an array as if it contained the entries of a diagonal matrix
! being multiplied from the left.
subroutine diagonal_lmult_tridiag(self, diag_array)

  class(TriDiagOp), intent(inout) :: self
  real(r8), USE_CONTIGUOUS intent(in) :: diag_array(:,:)

  self%spr = self%spr * diag_array(:,:self%ncel-1)
  self%sub = self%sub * diag_array(:,2:)
  self%diag = self%diag * diag_array(:,:)

  self%left_bound = self%left_bound * diag_array(:,1)
  self%right_bound = self%right_bound * diag_array(:,self%ncel)

end subroutine diagonal_lmult_tridiag

! Decomposition constructor
!
! The equation to be solved later (with left_div) is:
!     - A(k)*q(k+1) + B(k)*q(k) - C(k)*q(k-1) = D(k)
!
! The solution (effectively via LU decomposition) has the form:
!     E(k) = C(k) / (B(k) - A(k)*E(k+1))
!     F(k) = (D(k) + A(k)*F(k+1)) / (B(k) - A(k)*E(k+1))
!     q(k) = E(k) * q(k-1) + F(k)
!
! Unlike Richtmyer and Morton, E and F are defined by iterating backward
! down to level 1, and then q iterates forward.
!
! E can be calculated and stored now, without knowing D.
! To calculate F later, we store A and the denominator.
function new_TriDiagDecomp(op, graft_decomp) result(decomp)
  type(TriDiagOp), intent(in) :: op
  type(TriDiagDecomp), intent(in), optional :: graft_decomp

  type(TriDiagDecomp) :: decomp

  interface
     subroutine new_tridiagdecomp_codon(ca, dnom, ze, op_spr, op_sub, &
          op_diag, op_left_bound, op_right_bound, nsys, ncel) &
          bind(c, name='new_tridiagdecomp_codon')
       import :: c_int64_t, r8
       real(r8), intent(inout) :: ca(*), dnom(*), ze(*)
       real(r8), intent(in) :: op_spr(*), op_sub(*), op_diag(*)
       real(r8), intent(in) :: op_left_bound(*), op_right_bound(*)
       integer(c_int64_t), value :: nsys, ncel
     end subroutine new_tridiagdecomp_codon
     subroutine new_tridiagdecomp_graft_codon(ca, dnom, ze, op_spr, &
          op_sub, op_diag, op_left_bound, op_right_bound, graft_ca, &
          graft_dnom, graft_ze, nsys, op_ncel, decomp_ncel) &
          bind(c, name='new_tridiagdecomp_graft_codon')
       import :: c_int64_t, r8
       real(r8), intent(inout) :: ca(*), dnom(*), ze(*)
       real(r8), intent(in) :: op_spr(*), op_sub(*), op_diag(*)
       real(r8), intent(in) :: op_left_bound(*), op_right_bound(*)
       real(r8), intent(in) :: graft_ca(*), graft_dnom(*), graft_ze(*)
       integer(c_int64_t), value :: nsys, op_ncel, decomp_ncel
     end subroutine new_tridiagdecomp_graft_codon
  end interface
  logical, save :: new_tridiagdecomp_codon_logged = .false.
  logical, save :: new_tridiagdecomp_native_logged = .false.
  integer :: k

  if (present(graft_decomp)) then
     decomp%nsys = graft_decomp%nsys
     decomp%ncel = graft_decomp%ncel
  else
     decomp%nsys = op%nsys
     decomp%ncel = op%ncel
  end if

  ! Simple allocation with no error checking.
  allocate(decomp%ca(decomp%nsys,decomp%ncel))
  allocate(decomp%dnom(decomp%nsys,decomp%ncel))
  allocate(decomp%ze(decomp%nsys,decomp%ncel))

  if (linear_1d_use_native()) then
     ! decomp%ca is simply the negative of the tridiagonal's superdiagonal.
     decomp%ca(:,:op%ncel-1) = -op%spr
     decomp%ca(:,op%ncel) = -op%right_bound

     if (present(graft_decomp)) then
        ! Copy in graft_decomp beyond op%ncel.
        decomp%ca(:,op%ncel+1:) = graft_decomp%ca(:,op%ncel+1:)
        decomp%dnom(:,op%ncel+1:) = graft_decomp%dnom(:,op%ncel+1:)
        decomp%ze(:,op%ncel+1:) = graft_decomp%ze(:,op%ncel+1:)
        ! Fill in dnom edge value.
        decomp%dnom(:,op%ncel) = 1._r8 / (op%diag(:,op%ncel) - &
             decomp%ca(:,op%ncel)*decomp%ze(:,op%ncel+1))
     else
        ! If no grafting, the edge value of dnom comes from the diagonal.
        decomp%dnom(:,op%ncel) = 1._r8 / op%diag(:,op%ncel)
     end if

     do k = op%ncel - 1, 1, -1
        decomp%ze(:,k+1)   = - op%sub(:,k) * decomp%dnom(:,k+1)
        decomp%dnom(:,k) = 1._r8 / &
             (op%diag(:,k) - decomp%ca(:,k)*decomp%ze(:,k+1))
     end do

     ! Don't multiply edge level by denom, because we want to leave it up to
     ! the BoundaryCond object to decide what this means in left_div.
     decomp%ze(:,1) = -op%left_bound
     if (.not. new_tridiagdecomp_native_logged) then
        write(iulog,*) 'new_tridiagdecomp implementation = native'
        new_tridiagdecomp_native_logged = .true.
     end if
     return
  end if

  if (present(graft_decomp)) then
     call new_tridiagdecomp_graft_codon(decomp%ca, decomp%dnom, &
          decomp%ze, op%spr, op%sub, op%diag, op%left_bound, &
          op%right_bound, graft_decomp%ca, graft_decomp%dnom, &
          graft_decomp%ze, int(decomp%nsys, c_int64_t), &
          int(op%ncel, c_int64_t), int(decomp%ncel, c_int64_t))
  else
     call new_tridiagdecomp_codon(decomp%ca, decomp%dnom, decomp%ze, &
          op%spr, op%sub, op%diag, op%left_bound, op%right_bound, &
          int(decomp%nsys, c_int64_t), int(decomp%ncel, c_int64_t))
  end if
  if (.not. new_tridiagdecomp_codon_logged) then
     write(iulog,*) 'new_tridiagdecomp implementation = codon'
     new_tridiagdecomp_codon_logged = .true.
  end if

end function new_TriDiagDecomp

! Left-division (multiplication by inverse) using a decomposed operator.
!
! See the comment above for the constructor for a quick explanation of the
! intermediate variables. The "q" argument is "D(k)" on input and "q(k)" on
! output.
subroutine decomp_left_div(decomp, q, l_cond, r_cond)

  ! Decomposed matrix.
  class(TriDiagDecomp), intent(in) :: decomp
  ! Data to left-divide by the matrix.
  real(r8), USE_CONTIGUOUS intent(inout) :: q(:,:)
  ! Objects representing boundary conditions.
  class(BoundaryCond), intent(in), optional :: l_cond, r_cond

  ! "F" from the equation above.
  real(r8) :: zf(decomp%nsys,decomp%ncel)
  real(r8) :: l_edge_data(decomp%nsys), r_edge_data(decomp%nsys)

  integer(c_int64_t) :: has_l_cond, has_r_cond, l_cond_type, r_cond_type
  integer :: k
  interface
     subroutine decomp_left_div_codon(ca, dnom, ze, q, zf, l_edge_data, &
          r_edge_data, nsys, ncel, has_l_cond, l_cond_type, &
          has_r_cond, r_cond_type) bind(c, name='decomp_left_div_codon')
       import :: c_int64_t, r8
       real(r8), intent(in) :: ca(*), dnom(*), ze(*)
       real(r8), intent(inout) :: q(*), zf(*)
       real(r8), intent(in) :: l_edge_data(*), r_edge_data(*)
       integer(c_int64_t), value :: nsys, ncel
       integer(c_int64_t), value :: has_l_cond, l_cond_type
       integer(c_int64_t), value :: has_r_cond, r_cond_type
     end subroutine decomp_left_div_codon
  end interface
  logical, save :: decomp_left_div_codon_logged = .false.
  logical, save :: decomp_left_div_native_logged = .false.

  if (linear_1d_use_native()) then
     ! Include boundary conditions.
     if (present(l_cond)) then
        q(:,1) = q(:,1) + l_cond%apply_left(decomp%ze(:,1), q)
     end if

     if (present(r_cond)) then
        q(:,decomp%ncel) = q(:,decomp%ncel) + &
             r_cond%apply_right(decomp%ca(:,decomp%ncel), q)
     end if

     zf(:,decomp%ncel) = q(:,decomp%ncel) * decomp%dnom(:,decomp%ncel)

     do k = decomp%ncel - 1, 1, -1
        zf(:,k) = (q(:,k) + decomp%ca(:,k)*zf(:,k+1)) * decomp%dnom(:,k)
     end do

     ! Perform back substitution

     q(:,1) = zf(:,1)

     do k = 2, decomp%ncel
        q(:,k) = zf(:,k) + decomp%ze(:,k)*q(:,k-1)
     end do
     if (.not. decomp_left_div_native_logged) then
        write(iulog,*) 'decomp_left_div implementation = native'
        decomp_left_div_native_logged = .true.
     end if
     return
  end if

  has_l_cond = 0_c_int64_t
  has_r_cond = 0_c_int64_t
  l_cond_type = 0_c_int64_t
  r_cond_type = 0_c_int64_t
  l_edge_data = 0._r8
  r_edge_data = 0._r8

  if (present(l_cond)) then
     has_l_cond = 1_c_int64_t
     l_cond_type = int(l_cond%cond_type, c_int64_t)
     if (allocated(l_cond%edge_data)) l_edge_data = l_cond%edge_data
  end if

  if (present(r_cond)) then
     has_r_cond = 1_c_int64_t
     r_cond_type = int(r_cond%cond_type, c_int64_t)
     if (allocated(r_cond%edge_data)) r_edge_data = r_cond%edge_data
  end if

  call decomp_left_div_codon(decomp%ca, decomp%dnom, decomp%ze, q, zf, &
       l_edge_data, r_edge_data, int(decomp%nsys, c_int64_t), &
       int(decomp%ncel, c_int64_t), has_l_cond, l_cond_type, has_r_cond, &
       r_cond_type)
  if (.not. decomp_left_div_codon_logged) then
     write(iulog,*) 'decomp_left_div implementation = codon'
     decomp_left_div_codon_logged = .true.
  end if

end subroutine decomp_left_div

! Decomposition deallocation.
subroutine decomp_finalize(decomp)
  class(TriDiagDecomp), intent(inout) :: decomp
  integer(c_int64_t) :: dims(2)
  interface
     subroutine finalize_dims_codon(dims_p) bind(c, name='finalize_dims_codon')
       import :: c_int64_t
       integer(c_int64_t), intent(inout) :: dims_p(*)
     end subroutine finalize_dims_codon
  end interface
  logical, save :: decomp_finalize_codon_logged = .false.
  logical, save :: decomp_finalize_native_logged = .false.

  if (linear_1d_use_native()) then
     decomp%nsys = 0
     decomp%ncel = 0

     if (allocated(decomp%ca)) deallocate(decomp%ca)
     if (allocated(decomp%dnom)) deallocate(decomp%dnom)
     if (allocated(decomp%ze)) deallocate(decomp%ze)
     if (.not. decomp_finalize_native_logged) then
        write(iulog,*) 'decomp_finalize implementation = native'
        decomp_finalize_native_logged = .true.
     end if
     return
  end if

  dims = [int(decomp%nsys, c_int64_t), int(decomp%ncel, c_int64_t)]
  call finalize_dims_codon(dims)
  decomp%nsys = int(dims(1))
  decomp%ncel = int(dims(2))
  if (.not. decomp_finalize_codon_logged) then
     write(iulog,*) 'decomp_finalize implementation = codon'
     decomp_finalize_codon_logged = .true.
  end if

  if (allocated(decomp%ca)) deallocate(decomp%ca)
  if (allocated(decomp%dnom)) deallocate(decomp%dnom)
  if (allocated(decomp%ze)) deallocate(decomp%ze)

end subroutine decomp_finalize

end module linear_1d_operators
