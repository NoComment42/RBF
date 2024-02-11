module helper_mod
  use iso_c_binding
  use iso_fortran_env
  implicit none
  private

  public :: make_regular_box, write_points, write_points_wscalar

contains

  subroutine make_regular_box (npts, dims, box)
    use iso_c_binding
    implicit none
    ! dummy vars
    integer(c_int), intent(in) :: npts
    real(c_double), dimension(2), intent(in) :: dims
    real(c_double), dimension(:,:), allocatable, intent(out) :: box
    ! local
    integer(c_int) :: i, j, k, p
    real(c_double) :: dx
    !
    dx = (dims(2) - dims(1))/real(npts-1,c_double)

    allocate( box(3, npts**3) )

    p=0
    do i=0,npts-1
      do j=0,npts-1
        do k=0,npts-1
          p = p + 1
          box(1,p) = dims(1) + real(i,c_double)*dx
          box(2,p) = dims(1) + real(j,c_double)*dx
          box(3,p) = dims(1) + real(k,c_double)*dx
        end do
      end do
    end do

  end subroutine make_regular_box


  subroutine write_points (filename, points)
    ! dummy vars
    character(len=*), intent(in) :: filename
    real(c_double), dimension(:,:), intent(in) :: points
    ! local
    integer :: fid, i, npts
    !
    open(newunit=fid, file=trim(filename), status='unknown', action='write')

    npts = size(points,2)

    write(fid,'("x,y,z")')
    write(fid,'(2(ES15.7,",",X),ES15.7)') (points(:,i), i=1,npts)
    close(fid)
  end subroutine write_points

  subroutine write_points_wscalar (filename, points, scalar)
    ! dummy vars
    character(len=*), intent(in) :: filename
    real(c_double), dimension(:,:), intent(in) :: points
    real(c_double), dimension(:), intent(in) :: scalar
    ! local
    integer :: fid, i, npts
    !
    open(newunit=fid, file=trim(filename), status='unknown', action='write')

    npts = size(points,2)

    write(fid,'("x,y,z,var1")')
    do i=1,npts
      write(fid,'(3(ES15.7,",",X),ES15.7)') points(:,i), scalar(i)
    end do
    close(fid)
  end subroutine write_points_wscalar

end module helper_mod
