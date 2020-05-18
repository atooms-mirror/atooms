 
module interaction

  use helpers
  use cutoff  !, only: is_zero, smooth, tailor
  use potential  !, only: compute
  
  implicit none

contains

  subroutine forces(box,pos,ids,neighbors,number_neighbors,for,epot,virial)
    double precision, intent(in)    :: box(:)
    double precision, intent(in)    :: pos(:,:)
    integer,          intent(in)    :: ids(:)
    integer,          intent(in)    :: neighbors(:,:), number_neighbors(:)
    double precision, intent(inout) :: for(size(pos,1),size(pos,2))
    double precision, intent(out)   :: epot, virial
    double precision                :: rij(size(pos,1)), rijsq, uij, wij, hbox(size(pos,1)), hij
    integer                         :: i, j, isp, jsp, jn
    logical                         :: zero_ij
    ! TODO: it should be possible to tailor the cutoff from python, but the compute interface is not parsed correctly
    call tailor(compute)
    for = 0.0d0
    epot = 0.0d0
    virial = 0.0d0
    hbox = box / 2
    do i = 1,size(pos,2)
       isp = ids(i)
       do jn = 1,number_neighbors(i)
          j = neighbors(jn)
          !if (newton) then
          !   if (j<i) cycle
          !end if
          jsp = ids(j)
          call distance(i,j,pos,rij)
          call pbc(rij,box,hbox)
          call dot(rij,rij,rijsq)
          call is_zero(isp,jsp,rijsq,zero_ij)          
          if (.not.zero_ij) then
             ! TODO: remove hij via interface             
             call compute(isp,jsp,rijsq,uij,wij,hij) ! wij = -(du/dr)/r
             call smooth(isp,jsp,rijsq,uij,wij,hij) ! wij = -(du/dr)/r
             !print*, isp, jsp, rijsq**0.5, uij, wij
             epot = epot + uij
             virial = virial + wij * rijsq
             for(:,i) = for(:,i) + wij * rij
             for(:,j) = for(:,j) - wij * rij
          end if
       end do
    end do
  end subroutine forces

  subroutine hessian(box,pos,ids,neighbors,number_neighbors,hes)
    double precision, intent(in)    :: box(:)
    double precision, intent(in)    :: pos(:,:)
    integer,          intent(in)    :: ids(:)
    integer,          intent(in)    :: ids(:)
    integer,          intent(in)    :: neighbors(:,:), number_neighbors(:)
    double precision, intent(inout) :: hes(size(pos,1),size(pos,2),size(pos,1),size(pos,2))
    double precision                :: rij(size(pos,1)), rijsq, uij, wij, wwij, hbox(size(pos,1))
    double precision                :: unity(size(pos,1),size(pos,1)), mij(size(pos,1),size(pos,1)), mmij(size(pos,1),size(pos,1))
    integer                         :: i, j, isp, jsp
    logical                         :: zero_ij
    call tailor(compute)
    hes = 0.0d0
    unity = 0.0d0
    hbox = box / 2
    do i = 1,size(unity,1)
       unity(i,i) = 1.0d0
    end do
    loop_i: do i = 1,size(pos,2)
       isp = ids(i)
       loop_j: do jn = 1,number_neighbors(i)
          j = neighbors(jn)
          jsp = ids(j)
          rij = pos(:,i) - pos(:,j)
          call pbc(rij,box,hbox)
          rijsq = dot_product(rij,rij)
          call is_zero(isp,jsp,rijsq,zero_ij)
          if (.not.zero_ij) then
             call compute(isp,jsp,rijsq,uij,wij,wwij)
             call smooth(isp,jsp,rijsq,uij,wij,wwij)
             mij = unity(:,:) * wij
             mmij = outer_product(rij,rij) * wwij
             ! Diagonal in i,j - diagonal in x,y
             hes(:,i,:,i) = hes(:,i,:,i) - mij
             hes(:,j,:,j) = hes(:,j,:,j) - mij
             ! Off-diagonal in i,j - diagonal in x,y
             hes(:,i,:,j) = hes(:,i,:,j) + mij
             hes(:,j,:,i) = hes(:,j,:,i) + mij
             ! Diagonal in i,j - off-diagonal in x,y
             hes(:,i,:,i) = hes(:,i,:,i) + mmij
             hes(:,j,:,j) = hes(:,j,:,j) + mmij
             ! Off-diagonal in i,j - off-diagonal in x,y
             hes(:,i,:,j) = hes(:,i,:,j) - mmij
             hes(:,j,:,i) = hes(:,j,:,i) - mmij
          end if
       end do loop_j
    end do loop_i
  end subroutine hessian
  
end module interaction

