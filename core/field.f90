module field_mod
  use iso_c_binding
  implicit none

  type field_t
    real(c_double), dimension(:), allocatable :: scalar
  end type

contains

end module field_mod