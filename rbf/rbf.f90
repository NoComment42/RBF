module rbf_mod
  use iso_c_binding
  use iso_fortran_env
  use field_mod
  implicit none
  private

  public :: rbf_function, rbf_interp, multiquadric

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

  type, extends(rbf_function) :: multiquadric
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
    logical :: norm = .false.
  contains
    procedure, pass :: init
    procedure, pass :: interp
    procedure, private, pass :: calc_weights
    procedure, private, nopass :: radius
  end type rbf_interp
  type, extends(rbf_interp) :: rbf_interp_ext
    private
    real(c_double), dimension(:), allocatable :: poly_weights
  contains
    procedure, pass :: calc_weights => calc_weights_rbfe
  end type

  interface
    module subroutine init (this, field, func, norm)
      import, all
      ! dummy vars
      class(rbf_interp), intent(inout) :: this
      class(scalar_field), intent(in) :: field
      class(rbf_function), intent(in), optional :: func
      logical, intent(in), optional :: norm
    end subroutine init

    module function interp (this, pt) result(val)
      import, all
      ! dummy vars
      class(rbf_interp), intent(in) :: this
      real(c_double), dimension(3), intent(in) :: pt
      ! result
      real(c_double) :: val
    end function interp

    module subroutine calc_weights (this) 
      use lapack95
      import, all
      ! dummy vars
      class(rbf_interp), intent(inout) :: this
    end subroutine calc_weights

    module subroutine calc_weights_rbfe (this)
      use lapack95
      import, all
      ! dummy vars
      class(rbf_interp_ext), intent(inout) :: this
    end subroutine calc_weights_rbfe    

  end interface

contains

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
    class(multiquadric), intent(inout) :: this
    real(c_double), intent(in) :: r0
    !
    this%r0 = r0
  end subroutine mq_set_r0

  function mq_rbf (this, dist) result(val)
    ! dummy vars
    class(multiquadric), intent(in) :: this
    real(c_double), intent(in) :: dist
    ! result
    real(c_double) :: val
    !
    val = sqrt(dist**2 + this%r0**2)
  end function



!*************** Test procedure for MKL **************************!
  subroutine solve ()
    use lapack95
    real(c_double), dimension(5,5) :: &
    a=reshape( (/ 6.80,-2.11, 5.66, 5.97, 8.23, &
                 -6.05,-3.30, 5.36,-4.44, 1.08, &
                 -0.45, 2.58,-2.70, 0.27, 9.04, &
                  8.32, 2.71, 4.35,-7.17, 2.14, &
                 -9.67,-5.14,-7.26, 6.08,-6.87 /), &
                 (/5,5/) )
    real(c_double), dimension(5,3) :: &
    b=reshape( (/ 4.02, 6.19,-8.22,-7.57,-3.03, &
                 -1.56, 4.00,-8.67, 1.75, 2.86, &
                  9.81,-4.09,-4.57,-8.61, 8.99 /), &
                  (/5,3/) )
    integer, dimension(5) :: ipiv
    integer :: info, j

    !call dgesv(5, 3, a, 5, ipiv, b, 5, info)
    call getrf(a, ipiv, info)

    write(*,*) info
    if (info .gt. 0) then
      write(*,'( "u(" ,I3, ",", I3, ") is zero" )') info, info
    end if
    write(*,'(5(i3, " "))') ipiv

    call getrs(a, ipiv, b, info=info)

    write(*,*) info
    if (info .gt. 0) then
      write(*,'( "u(" ,I3, ",", I3, ") is zero" )') info, info
    end if
    write(*,'(5(i3, " "))') ipiv

    write(*,'(3(ES15.7, " "))') ( b(j,:), j=1,5)
  end subroutine solve

  end module rbf_mod