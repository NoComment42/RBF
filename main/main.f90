program interpolation
  use iso_c_binding
  use iso_fortran_env, only: stdout => output_unit
  use field_mod
  use rbf_mod
  use inverse_multiquadric_mod
  use helper_mod
  implicit none

  interface
    subroutine checkRBF(rbf, npts, pc, vals)
      use iso_c_binding
      use rbf_mod
      implicit none
      type(rbf_interp), intent(in) :: rbf
      integer, intent(in) :: npts
      real(c_double), dimension(:,:), allocatable, intent(out) :: pc
      real(c_double), dimension(:), allocatable, intent(out) :: vals
    end subroutine
  end interface

  integer :: npts
  real(c_double), dimension(:,:), allocatable :: box, pc
  real(c_double), dimension(:), allocatable :: scalar, vals
  type(scalar_field) :: field
  type(rbf_interp) :: rbf, rbf_normed
  type(inverse_multiquadric) :: inv_mq


  write(stdout,*) "Starting interpolation test bed!"


  write(stdout,*) "Make regularly space box of points....."
  call make_regular_box(2, [0._c_double, 1._c_double], box)

  write(stdout,*) "Add scalar field to box....."
  write(stdout,*) size(box,1), size(box,2)
  allocate( scalar(size(box,2)) )
  scalar = 1.0_c_double

  call write_points_wscalar("box.csv", box, scalar)

  write(stdout,*) "Initialize field from box....."
  call field%init(box,scalar)
  call write_points_wscalar("field.csv", field%get_vertices(), field%get_values())

  write(stdout,*) "Setup radial basis function interpolator....."
  call rbf%init(field, func=inv_mq , norm=.false.)
  call rbf_normed%init(field, func=inv_mq, norm=.true.)

  write(stdout,*) "Test interpolators"
  
  npts = 11
  call checkRBF(rbf, npts, pc, vals)
  call write_points_wscalar("interped.csv", pc, vals)

  call checkRBF(rbf_normed, npts, pc, vals)
  call write_points_wscalar("interped-normed.csv", pc, vals)

end program

subroutine checkRBF(rbf, npts, pc, vals)
  use iso_c_binding
  use rbf_mod
  implicit none
  type(rbf_interp), intent(in) :: rbf
  integer, intent(in) :: npts
  real(c_double), dimension(:,:), allocatable, intent(inout) :: pc
  real(c_double), dimension(:), allocatable, intent(out) :: vals
  !
  real(c_double), dimension(3) :: pt, step
  integer :: i
  !

  if(allocated(pc)) deallocate(pc)
  if(allocated(vals)) deallocate(vals)
  allocate( pc(3,npts) )
  allocate( vals(npts) )
  pt = 0.0_c_double
  step = 0.1_c_double
  pt = pt - step
  do i=1,npts
    pt = pt + step
    vals(i) = rbf%interp(pt)
    pc(:,i) = pt
  end do

end subroutine

