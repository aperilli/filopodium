subroutine new_combined_list(rlist,nmon,pos,poslist,nlist,list)

use common_variables

  implicit none
  integer, dimension(nfil), intent(in) :: nmon
  real(kind=8), intent(in) :: rlist
  real(kind=8), dimension(3,max_size_fil,nfil), intent(in) :: pos
  real(kind=8), dimension(3,max_size_fil,nfil), intent(out) :: poslist
  integer, intent(out) :: list(nfil*max_size_fil,20),nlist(nfil*max_size_fil)
  integer, allocatable :: head(:,:,:),linked_list(:)
  integer :: ll(nfil*max_size_fil),ncellx,ncelly,ncellz
  integer :: n,i,j,jlist,ilist,icell(3),jcellx,jcelly,jcellz,m
  real(kind=8) :: dr(3),dist,xmin,ymin,ymax,zmin,zmax,rn(3)

  xmin=minval(pos(1,1,:))-0.5d0*delta_inter
  ymin=minval(pos(2,1,:))-0.5d0*delta_inter
  ymax=maxval(pos(2,1,:))+0.5d0*delta_inter
  zmin=minval(pos(3,1,:))-0.5d0*delta_inter
  zmax=maxval(pos(3,1,:))+0.5d0*delta_inter
  ncellx=int(abs(box_size(1)-xmin)/rcut_inter)+1
  ncelly=int(abs(ymax-ymin)/rcut_inter)+1
  ncellz=int(abs(zmax-zmin)/rcut_inter)+1
  allocate(head(ncellx,ncelly,ncellz))
  allocate(linked_list(sum(nmon(:))))
  call new_cell_list(nmon,rcut_inter,pos,ncellx,ncelly,ncellz,linked_list,head)
  do i=1,ncellx
    do j=1,ncelly
      write(20,*)i,j,head(:,j,i)
      write(20,*)
      write(20,*)
    end do
  end do
  write(30,*)'list:',linked_list(1:sum(nmon(1:nfil)))

!  call new_cell_list(nmon,rlist,pos,nx,ny,nz,ll,head)
  nlist(:)=0
  list(:,:)=0
  poslist(:,:,:)=pos(:,:,:)
  do n=1,nfil
    do i =1,nmon(n)
      icell(1)=int((pos(1,i,n)-xmin)/rn(1))+1
      icell(2)=int((pos(2,i,n)-ymin)/rn(2))+1
      icell(3)=int((pos(3,i,n)-zmin)/rn(3))+1
      ilist=i
      do m=1,n-1
        ilist=ilist+nmon(m)
      end do
      do jcellx=max(icell(1)-1,1),max(icell(1)+1,ncellx)
        do jcelly=max(icell(2)-1,1),max(icell(2)+1,ncelly)
          do jcellz=max(icell(3)-1,1),max(icell(3)+1,ncellz)
            jlist=head(jcellx,jcelly,jcellz)
            do while (jlist.ne.0)
              if(jlist.ne.ilist) then
                m=0
                j=jlist
                do while (j.gt.0)
                  m=m+1
                  j=j-nmon(m)
                end do
                j=j+nmon(m)
                dr(:)=pos(:,i,n)-pos(:,j,m)
                dist=sqrt(sum(dr(:)**2))
                if (dist.lt.rlist) then
                  nlist(ilist)=nlist(ilist)+1
                  list(ilist,nlist(ilist))=j
                end if
                jlist=linked_list(jlist)
              end if
            end do
          end do
        end do
      end do
    end do
  end do

return

end subroutine new_combined_list
