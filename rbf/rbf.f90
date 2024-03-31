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

contains


  subroutine init (this, field, func, norm)
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


  function interp (this, pt) result(val)
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


  subroutine calc_weights (this) 
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