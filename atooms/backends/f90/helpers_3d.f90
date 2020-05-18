module helpers
  
  implicit none

contains

  pure subroutine pbc(r,box,hbox)
    double precision, intent(inout) :: r(:)
    double precision, intent(in)    :: box(:), hbox(:)
    if (abs(r(1)) > hbox(1)) r(1) = r(1) - sign(box(1),r(1))
    if (abs(r(2)) > hbox(2)) r(2) = r(2) - sign(box(2),r(2))
    if (abs(r(3)) > hbox(3)) r(3) = r(3) - sign(box(3),r(3))
  end subroutine pbc

  pure subroutine dot(r1,r2,out)
    double precision, intent(in)  :: r1(:), r2(:)
    double precision, intent(out) :: out
    out = r1(1)*r2(1) + r1(2)*r2(2) + r1(3)*r2(3)
  end subroutine dot

  pure subroutine distance(i,j,pos,rij)
    integer, intent(in) :: i, j
    double precision, intent(in)    :: pos(:,:)
    double precision, intent(inout) :: rij(:)
    rij(1) = pos(1,i) - pos(1,j)
    rij(2) = pos(2,i) - pos(2,j)
    rij(3) = pos(3,i) - pos(3,j)
  end subroutine distance

  pure function outer_product(x,y) result (out)
    double precision, intent(in) :: x(:), y(:)
    double precision :: out(size(x),size(x))
    integer :: i, j
    do i = 1,size(x)
       do j = 1,size(x)
          out(i,j) = x(i) * y(j)
       end do
    end do
  end function outer_product

end module helpers
