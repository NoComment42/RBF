submodule (rbf_mod) rbf_extended
contains

  module subroutine calc_weights_rbfe (this)
    use lapack95
    ! dummy vars
    class(rbf_interp_ext), intent(inout) :: this
    ! local 
    integer(c_int) :: i, j
    ! Note that these are 'normal' fortran integers not c_int sized.
    integer :: info, npts 
    integer, dimension(:), allocatable :: ipiv
    real(c_double) :: sum
    real(c_double), dimension(3) :: p1, p2
    real(c_double), dimension(:,:), allocatable :: rbf, poly, rhs
    !
    npts = this%field%count()

    allocate( rhs(npts,1) )
    allocate( ipiv(npts) )
    allocate( rbf(npts,npts) )   
  end subroutine

end submodule rbf_extended