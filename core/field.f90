module field_mod
  use iso_c_binding
  implicit none
  private

  public :: scalar_field

  type scalar_field
    private
    real(c_double), dimension(:,:), allocatable :: vertices
    real(c_double), dimension(:), allocatable :: scalar
    integer(c_int) :: num_points
  contains
    procedure, pass :: init => sf_init
    procedure, pass :: count => sf_count
    procedure, pass :: get_vertex => sf_get_vertex
    procedure, pass :: get_vertices => sf_get_vertices
    procedure, pass :: get_value => sf_get_value
    procedure, pass :: get_values => sf_get_values
  end type

  interface scalar_field
    module procedure new_scalar_field
  end interface

contains

  !> Constructor function that returns an instance of a scalar_field object.
  !!
  !! @param[in] vertices (c_double,c_double) 2D array (3,npoints) that contains the coordinates of the vertices.
  !! @param[in] scalar (c_double) 1D array (npoints) of the scalar field associated with the vertices.
  !! @result scalar_field Scalar_field object.
  function new_scalar_field (vertices, scalar) result(field)
    ! dummy vars
    real(c_double), dimension(:,:), intent(in) :: vertices
    real(c_double), dimension(size(vertices,1)), intent(in) :: scalar
    ! result
    type(scalar_field) :: field
    !
    call field%init(vertices, scalar)
  end function

  !> Initiailizes a scalar_field object given its vertex locations and field value.
  !!
  !! @param[in] vertices (c_double,c_double) 2D array (3,npoints) that contains the coordinates of the vertices.
  !! @param[in] scalar (c_double) 1D array (npoints) of the scalar field associated with the vertices.
  subroutine sf_init (this, vertices, scalar)
    ! dummy vars
    class(scalar_field), intent(inout) :: this
    real(c_double), dimension(:,:), intent(in) :: vertices
    real(c_double), dimension(:), intent(in) :: scalar
    !
    allocate( this%vertices, source=vertices)
    this%scalar = scalar
    this%num_points = size(vertices,2)
  end subroutine sf_init

  !> Returns the number of points in the field.
  !!
  !! @result (int) The number of points in this field.
  function sf_count (this) result(count)
    ! dummy vars
    class(scalar_field), intent(in) :: this
    ! result
    integer(c_int) :: count
    !
    count = this%num_points
  end function sf_count


  function sf_get_vertex (this, i) result(vertex)
    ! dummy vars
    class(scalar_field), intent(in) :: this
    integer(c_int), intent(in) :: i
    ! result
    real(c_double), dimension(3) :: vertex
    !
    vertex = this%vertices(:,i)
  end function sf_get_vertex


  function sf_get_vertices (this) result(vertices)
    ! dummy vars
    class(scalar_field), intent(in) :: this
    ! result
    real(c_double), dimension(:,:), allocatable :: vertices
    !
    vertices = this%vertices
  end function sf_get_vertices


  function sf_get_value (this, i) result(val)
    ! dummy vars
    class(scalar_field), intent(in) :: this
    integer(c_int), intent(in) :: i
    ! result
    real(c_double) :: val
    !
    val = this%scalar(i)
  end function sf_get_value


  function sf_get_values (this) result(vals)
    ! dummy vars
    class(scalar_field), intent(in) :: this
    ! result
    real(c_double), dimension(size(this%scalar)) :: vals
    !
    vals = this%scalar
  end function sf_get_values

end module field_mod