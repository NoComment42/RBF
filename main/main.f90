program interpolation
  use iso_c_binding
  use iso_fortran_env, only: stdout => output_unit
  use field_mod
  use helper_mod
  implicit none

  real(c_double), dimension(:,:), allocatable :: box

  write(stdout,*) "Starting interpolation test bed!"


  write(stdout,*) "Make regularly space box of points....."
  call make_regular_box(10, [-1._c_double, 1._c_double], box)
  call write_points("box.csv", box)

end program


