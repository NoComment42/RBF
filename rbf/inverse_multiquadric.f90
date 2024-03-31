module inverse_multiquadric_mod
  use rbf_mod
  use iso_c_binding
  implicit none
  private

  public :: inverse_multiquadric

  type, extends(rbf_function) :: inverse_multiquadric
    private
    real(c_double) :: r0 = 1.0_c_double
  contains
    procedure, pass :: set_r0 => imq_set_r0
    procedure, pass :: rbf => imq_rbf
  end type

contains

  subroutine imq_set_r0 (this,r0)
    ! dummy vars
    class(inverse_multiquadric), intent(inout) :: this
    real(c_double), intent(in) :: r0
    !
    this%r0 = r0
  end subroutine imq_set_r0

  function imq_rbf (this, dist) result(val)
    ! dummy vars
    class(inverse_multiquadric), intent(in) :: this
    real(c_double), intent(in) :: dist
    ! result
    real(c_double) :: val
    !
    val = sqrt(dist**2 + this%r0**2)**(-1)
  end function imq_rbf

end module