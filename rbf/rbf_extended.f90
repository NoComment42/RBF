submodule (rbf_mod) rbf_extended
contains

  module subroutine init_rbfe (this, field, degree, func, norm)
    ! dummy vars
    class(rbf_interp_ext), intent(inout) :: this
    class(scalar_field), intent(in) :: field
    integer(c_int), intent(in) :: degree
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

    this%degree = degree 
    allocate( this%poly_weights(this%degree + 1) )
    this%poly_weights = 0.0_c_double

    if (present(norm)) this%norm = norm

    call this%calc_weights()
  end subroutine init_rbfe

  module subroutine calc_weights_rbfe (this)
    use lapack95
    ! dummy vars
    class(rbf_interp_ext), intent(inout) :: this
    ! local 
    integer(c_int) :: i, j
    ! Note that these are 'normal' fortran integers not c_int sized.
    integer :: info, npts, dim
    integer, dimension(:), allocatable :: ipiv
    real(c_double) :: sum
    real(c_double), dimension(3) :: p1, p2
    real(c_double), dimension(:,:), allocatable :: rbf, poly, rhs
    !
    npts = this%field%count()

    dim = npts + 3*this%degree + 1

    allocate( rhs(dim, 1) )
    allocate( ipiv(dim) )
    allocate( rbf(dim,dim) )   

    rhs = 0.0_c_double
    rbf = 0.0_c_double

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

      do j=1, this%degree + 1
        rbf(i,j + npts) = p1**(j-1)
        rbf(j + npts,i) = rbf(i,j + npts)
      end do
    end do

    write(*,'(<dim>(ES15.7,X))') ( rbf(i,:), i=1,dim)

    ! Solve linear system RBF*weights = RHS
    call getrf(rbf, ipiv, info)
    call getrs(rbf, ipiv, rhs, info=info)
    this%weights = rhs(:,1)

  end subroutine calc_weights_rbfe

  module function eval_poly_rbfe (this, x) result(res)
    ! dummy vars
    class(rbf_interp_ext), intent(in) :: this
    real(c_double), dimension(3), intent(in) :: x
    ! result
    real(c_double) :: res
    ! local
    integer :: i, j
    !

    res = this%poly_weights(1)
    do i=2, this%degree + 1
      res = res + this%poly_weights(j)*x(1)**(i-1)
      res = res + this%poly_weights(j+1)*x(2)**(i-1)
      res = res + this%poly_weights(j+2)*x(3)**(i-1)
      j = j + 3
    end do
  end function eval_poly_rbfe



end submodule rbf_extended