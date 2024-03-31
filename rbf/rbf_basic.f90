submodule (rbf_mod) rbf_basic
contains

module subroutine init (this, field, func, norm)
  ! dummy vars
  class(rbf_interp), intent(inout) :: this
  class(scalar_field), intent(in) :: field
  class(rbf_function), intent(in), optional :: func
  logical, intent(in), optional :: norm
  ! local
  type(multiquadric) :: default
  !
  this%field = field

  if (present(func)) then
    this%func = func
  else
    this%func = default
  end if
  if (present(norm)) this%norm = norm

  call this%calc_weights()
end subroutine init


module function interp (this, pt) result(val)
  ! dummy vars
  class(rbf_interp), intent(in) :: this
  real(c_double), dimension(3), intent(in) :: pt
  ! result
  real(c_double) :: val
  ! local
  integer(c_int) :: i
  real(c_double) :: sum, fval
  real(c_double), dimension(3) :: p2
  !
  sum = 0._c_double
  val = 0._c_double
  do i=1,this%field%count()
    p2 = this%field%get_vertex(i)
    fval = this%func%rbf( this%radius(pt,p2) )
    val = val + this%weights(i)*fval
    sum = sum + fval
  end do
  if (this%norm) val = val/sum
end function interp


module subroutine calc_weights (this) 
  use lapack95
  ! dummy vars
  class(rbf_interp), intent(inout) :: this
  ! local
  integer(c_int) :: i, j
  ! Note that these are 'normal' fortran integers not c_int sized.
  integer :: info, npts 
  integer, dimension(:), allocatable :: ipiv
  real(c_double) :: sum
  real(c_double), dimension(3) :: p1, p2
  real(c_double), dimension(:,:), allocatable :: rbf, rhs
  !
  npts = this%field%count()

  allocate( rhs(npts,1) )
  allocate( ipiv(npts) )
  allocate( rbf(npts,npts) )

  do i=1, npts      
    p1 = this%field%get_vertex(i)

    sum = 0.0_c_double
    do j=1, npts
      p2 = this%field%get_vertex(j)
      rbf(i,j) = this%func%rbf( this%radius(p1,p2) )
      sum = sum + rbf(i,j)
    end do

    rhs(i,1) = this%field%get_value(i)
    if (this%norm) rhs(i,1) = rhs(i,1)*sum
  end do
  
  ! Solve linear system RBF*weights = RHS
  call getrf(rbf, ipiv, info)
  call getrs(rbf, ipiv, rhs, info=info)
  this%weights = rhs(:,1)
  write(*,*)
  write(*,*) "weights"
  write(*,*) this%weights
end subroutine calc_weights

end submodule rbf_basic