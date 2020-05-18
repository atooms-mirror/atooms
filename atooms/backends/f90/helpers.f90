module helpers
  
  implicit none

contains

  pure subroutine pbc(r,box,hbox)
    double precision, intent(inout) :: r(:)
    double precision, intent(in)    :: box(:), hbox(:)
    where (abs(r) > hbox)
       r = r - sign(box,r)
    end where
  end subroutine pbc

  pure subroutine distance(i,j,pos,rij)
    integer, intent(in) :: i, j
    double precision, intent(in)    :: pos(:,:)
    double precision, intent(inout) :: rij(:)
    rij = pos(:,i) - pos(:,j)
  end subroutine distance

  pure subroutine dot(r1,r2,out)
    double precision, intent(in)  :: r1(:), r2(:)
    double precision, intent(out) :: out
    out = dot_product(r1,r2)
  end subroutine dot
  
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
