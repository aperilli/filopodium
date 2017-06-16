subroutine new_cell_list(nmon,rc,pos,nx,ny,nz,ll,head)

use common_variables

  implicit none
  integer, dimension(nfil), intent(in) :: nmon
  integer, intent(in) :: nx,ny,nz
  real(kind=8), intent(in) :: rc 
  real(kind=8), dimension(3,max_size_fil,nfil), intent(in) :: pos
  integer, intent(out) :: ll(nfil*max_size_fil)
  integer, intent(out) :: head(1:nx,1:ny,1:nz)
  integer :: i,n,m,icell(3),ilist
  real(kind=8), dimension(3) :: rn
  real(kind=8) :: xmin,ymin,ymax,zmin,zmax

  head(:,:,:)=0

  xmin=minval(pos(1,1,:))-0.5d0*delta_inter
  ymin=minval(pos(2,1,:))-0.5d0*delta_inter
  ymax=maxval(pos(2,1,:))+0.5d0*delta_inter
  zmin=minval(pos(3,1,:))-0.5d0*delta_inter
  zmax=maxval(pos(3,1,:))+0.5d0*delta_inter
  rn(1)=(box_size(1)-xmin)/(int((box_size(1)-xmin)/rc)+1)
  rn(2)=(ymax-ymin)/(int((ymax-ymin)/rc)+1)
  rn(3)=(zmax-zmin)/(int((zmax-zmin)/rc)+1)
  do n=1,nfil
    do i=1,nmon(n)
      icell(1)=int((pos(1,i,n)-xmin)/rn(1))+1
      icell(2)=int((pos(2,i,n)-ymin)/rn(2))+1
      icell(3)=int((pos(3,i,n)-zmin)/rn(3))+1
      ilist=i
      do m=1,n-1
        ilist=ilist+nmon(m)
      end do
      ll(ilist)=head(icell(1),icell(2),icell(3))
      head(icell(1),icell(2),icell(3))=ilist
    end do
  end do

return

end subroutine new_cell_list
