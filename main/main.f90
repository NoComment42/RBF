program interpolation
  use iso_c_binding
  use iso_fortran_env, only: stdout => output_unit
  use field_mod
  use rbf_mod
  use helper_mod
  implicit none

  real(c_double), dimension(:,:), allocatable :: box
  real(c_double), dimension(:), allocatable :: scalar
  type(scalar_field) :: field
  type(rbf_interp) :: rbf

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

  write(stdout,*) "Setupd radial basis function interpolator....."
  call rbf%init(field)

end program


