 module neighbor_list

  ! use helpers
  
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
  
  subroutine compute(box,pos,ids,rcut,neighbors,number_neighbors,error)
    !! Compute neighbor lists using III Newton law
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

  logical function need_update_largest(displacement, skin)
    !! Update is needed if the largest displacement exceeds 1/2
    !! of the Verlet skin (from Allen & Tildesley)
    real(8), intent(in) :: displacement(:,:), skin
    real(8) :: dr(size(displacement,2))
    integer ::i
    dr = [(dot_product(displacement(:,i),displacement(:,i)), i=1,size(dr))]
    need_update_largest = maxval(dr)>(0.25d0 * skin**2)
  end function need_update_largest

  subroutine add_displacement(pos,pos_old,box,displacement)
    !! Add displacements of particles for a subset of particles specified by index.
    !! Assume that PBC has been applied, hence we need to refold them.
    real(8), intent(in)    :: pos(:,:)
    real(8), intent(inout) :: pos_old(:,:), displacement(:,:)
    real(8), intent(in)    :: box(:)
    real(8) :: hbox(size(box))
    integer :: i
    hbox = box / 2
    displacement(:,:) = displacement(:,:) + pos(:,:) - pos_old(:,:)
    do i = 1,size(pos,2)
       call pbc(displacement(:,i),box,hbox)
    end do
    pos_old(:,:) = pos(:,:)
  end subroutine add_displacement
  
end module neighbor_list

