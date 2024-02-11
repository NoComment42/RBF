module rbf_mod
  use iso_c_binding
  use iso_fortran_env
  use field_mod
  implicit none
  private

  public :: rbf_function, rbf_interp, rbf_multiquadric

  type, abstract :: rbf_function
  contains
    procedure(rbf_def), pass, deferred :: rbf
  end type rbf_function

  interface
    function rbf_def (this, dist) result(val)
      use iso_c_binding
      import rbf_function
      implicit none
      ! dummy vars
      class(rbf_function), intent(in) :: this
      real(c_double), intent(in) :: dist
      ! result
      real(c_double) :: val
    end function
  end interface

  type, extends(rbf_function) :: rbf_multiquadric
    private
    real(c_double) :: r0 = 1.0_c_double
  contains
    procedure, pass :: set_r0 => mq_set_r0
    procedure, pass :: rbf => mq_rbf
  end type

  type rbf_interp
    private
    class(rbf_function), allocatable :: func
    type(scalar_field) :: field
    real(c_double), dimension(:), allocatable :: weights
  contains
    procedure, pass :: init
    procedure, pass :: interp
    procedure, private, pass :: calc_weights
    procedure, private, nopass :: radius
  end type rbf_interp

contains


  subroutine init (this, field, func)
    ! dummy vars
    class(rbf_interp), intent(inout) :: this
    class(scalar_field), intent(in) :: field
    class(rbf_function), intent(in), optional :: func
    ! local
    type(rbf_multiquadric) :: default
    !
    this%field = field

    if (present(func)) then
      this%func = func
    else
      this%func = default
    end if

    call this%calc_weights()
  end subroutine init


  function interp (this, pt) result(val)
    ! dummy vars
    class(rbf_interp), intent(in) :: this
    real(c_double), dimension(3), intent(in) :: pt
    ! result
    real(c_double) :: val
    !
    val = 0._c_double
  end function interp


  subroutine calc_weights (this) 
    ! dummy vars
    class(rbf_interp), intent(inout) :: this
    ! local
    integer(c_int) :: i, j, npts
    real(c_double), dimension(3) :: p1, p2
    real(c_double), dimension(:), allocatable :: rhs
    real(c_double), dimension(:,:), allocatable :: rbf
    !
    npts = this%field%count()

    allocate( rhs(npts) )
    allocate( rbf(npts,npts) )

    do i=1, npts
      p1 = this%field%get_vertex(i)

      do j=1, npts
        p2 = this%field%get_vertex(j)
        rbf(i,j) = this%func%rbf( this%radius(p1,p2) )
        !write(*,*) i, j, this%radius(p1,p2), rbf(i,j)
      end do

      rhs(i) = this%field%get_value(i)
    end do
    
    ! Solve linear system RBF*weights = RHS
    
  end subroutine calc_weights


  function radius (p1, p2) result(rad)
    ! dummy vars
    real(c_double), dimension(3), intent(in) :: p1, p2
    ! result
    real(c_double) :: rad
    !
    rad = norm2(p2-p1)
  end function radius

!*************** Multiquadric Radial Basis Function ***************!

  subroutine mq_set_r0 (this,r0)
    ! dummy vars
    class(rbf_multiquadric), intent(inout) :: this
    real(c_double), intent(in) :: r0
    !
    this%r0 = r0
  end subroutine mq_set_r0

  function mq_rbf (this, dist) result(val)
    ! dummy vars
    class(rbf_multiquadric), intent(in) :: this
    real(c_double), intent(in) :: dist
    ! result
    real(c_double) :: val
    !
    val = sqrt(dist**2 + this%r0**2)
  end function
end module rbf_mod