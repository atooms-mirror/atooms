 
module neighbor_list

  use helpers
  
  implicit none

contains

  subroutine compute(box,pos,ids,rcut,neighbors,number_neighbors,error)
    double precision, intent(in)    :: box(:)
    double precision, intent(in)    :: pos(:,:), rcut(:,:)
    integer,          intent(in)    :: ids(:)
    integer,          intent(inout) :: neighbors(:,:), number_neighbors(size(pos,2))
    logical,          intent(out)   :: error
    double precision                :: rij(size(pos,1)), rijsq, hbox(size(pos,1))
    integer                         :: i, j, isp, jsp
    error = .false.
    number_neighbors = 0
    hbox = box / 2
    do i = 1,size(pos,2)
       isp = ids(i)
       do j = i+1,size(pos,2)
          jsp = ids(j)
          call distance(i,j,pos,rij)
          call pbc(rij,box,hbox)
          call dot(rij,rij,rijsq)
          if (rijsq <= rcut(isp,jsp)**2) then
             number_neighbors(i) = number_neighbors(i) + 1
             if (number_neighbors(i) <= size(neighbors,1)) then
                neighbors(number_neighbors(i),i) = j
             else
                error = .true.
             end if
          end if
       end do
    end do
  end subroutine compute
  
end module neighbor_list

