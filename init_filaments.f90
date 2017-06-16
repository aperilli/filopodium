subroutine init_filaments(qid,nmon,pos,vel,Lwall,vel_wall,t_in)

use common_variables
use st_mpi
use mpi_routines

  implicit none
  character(len=30), intent(in) :: qid 
  integer, dimension(nfil), intent(out) :: nmon 
  real(kind=8), dimension(nproc), intent(out):: t_in
  real(kind=8), dimension(3,max_size_fil,nfil), intent(out) :: pos,vel
  real(kind=8), intent(out) :: Lwall,vel_wall

  real(kind=8), dimension(3) :: rnd
  integer :: nn,iy,iz,m,icf,idum,l
  real(kind=8), allocatable :: y_ret(:),z_ret(:)
  real(kind=8) :: lcell,delta,dl,ltot,nfilsq,eta(2),z,dpoint,tdum
  real(kind=8) :: rng,rng_cm,Lwall_dum
  real(kind=8), dimension(3,max_size_fil,nfil) :: pos_dum,vel_dum
  external :: rng,rng_cm

  integer :: npoint,nshell,i,j,k,ipoint,ishell,n,nmon_dum(1:nfil)
  integer, allocatable :: vacancy(:)
  logical :: ifexist
  character(len=50) :: filen

  pos(:,:,:)=0.d0
  vel(:,:,:)=0.d0
  Lwall=0.d0
  vel_wall=0.d0
  t_in=0.d0

  filen=TRIM(qid)//'.rs'

  if(itype.ne.0) then
    inquire(file=TRIM(filen),exist=ifexist)
    if(ifexist) then
      open(2,file=TRIM(filen),form='formatted',status='old')
      do l=1,9999
        read(2,*,end=211) icf,tdum
      end do
211   if(root) write(*,*)'restart file found with ',l-1,'configurations'
      if(l-1.ne.nproc) then
        if(root) write(*,*) 'nconfs /= nproc',l-1,nproc
        if(root) write(*,*) 'if nconf > nproc continue but CHECK!!!'
        if (l-1<nproc) call parallel_abort()
      endif
      rewind(2)
    else
      if(root) write(*,*) 'itypes 1, 2 and 3 require a restart file'
      call parallel_abort()
    endif

    do idum=1,nproc
      read(2,*)icf,tdum,Lwall_dum,nmon_dum(1:nfil), &
            (((pos_dum(i,j,k),i=1,3), j=1,nmon_dum(k)),k=1,nfil), &
            (((vel_dum(i,j,k),i=1,3), j=1,nmon_dum(k)),k=1,nfil)
      if (icf.eq.myrank+1) then
        if(itype.eq.3)t_in(myrank+1)=tdum
!        print*,'initial time',t_in(myrank+1),myrank+1
        Lwall=Lwall_dum
        nmon(:)=nmon_dum(:)
        do j=1,nfil
          m=nmon(j)
          do i=1,m
            pos(:,i,j)=pos_dum(:,i,j)
            vel(:,i,j)=vel_dum(:,i,j)
          end do
        end do
      end if
    end do

!    print*,myrank+1,tdum,Lwall,nmon(1:nfil)

    close(2)

  elseif(itype.eq.0) then
    nmon(:)=3
    Lwall=nfil*kBT*log(rho1)/dl0/ktrap

    nshell=0
    npoint=0
    do while(npoint.lt.nfil)
      nshell=nshell+1
      npoint=3*nshell*(nshell-1)+1
    end do
    allocate(y_ret(1:npoint),z_ret(1:npoint),vacancy(1:npoint))
    m=0
    do i=-(nshell-1),(nshell-1)
      z=sqrt(3.d0)*i*delta_inter/2.d0
      do j=1,2*nshell-1-iabs(i)
        m=m+1
        z_ret(m)=z
        y_ret(m)=(-(2.d0*nshell-1.d0-abs(i))+(2.d0*j-1.d0))*delta_inter/2.d0
      end do
    end do

!    eta(1)=rng()
!    eta(2)=rng()
    pos(2,1:2,1)=y_ret((npoint+1)/2)!+(eta(1)*2.d0-1.d0)*dl0/2.d0
    pos(3,1:2,1)=z_ret((npoint+1)/2)!+(eta(2)*2.d0-1.d0)*dl0/2.d0
    pos(1,1,1)=(1.d0-0.5d0*(nfil+1.d0))*dl0/nfil
    pos(1,2,1)=pos(1,1,1)+dl0

    j=2
    do ishell=2,nshell
      do ipoint=1,npoint
        dpoint=y_ret(ipoint)**2+z_ret(ipoint)**2
        dpoint=sqrt(dpoint)
        if(dpoint.gt.(ishell-2.d0)*delta_inter &
         .and.dpoint.le.(ishell-1.d0)*delta_inter) then
!          eta(1)=rng()
!          eta(2)=rng()
          pos(2,1:2,j)=y_ret(ipoint)!+(eta(1)*2.d0-1.d0)*dl0/2.d0
          pos(3,1:2,j)=z_ret(ipoint)!+(eta(2)*2.d0-1.d0)*dl0/2.d0
          pos(1,1,j)=(j-0.5d0*(nfil+1.d0))*dl0/nfil
          pos(1,2,j)=pos(1,1,j)+dl0 
          vel(:,1,j)=0.d0
          vel(:,2,j)=0.d0
          j=j+1
        end if
        if(j.eq.nfil+1) goto 200
      end do
    end do

200  do j=1,nfil
      m=nmon(j)
      do i=3,m
        pos(1,i,j)=pos(1,i-1,j)+dl0
        pos(2:3,i,j)=pos(2:3,1,j)
        vel(:,i,j)=(2.d0*rng()-1.d0)*vmxd
      end do
    end do

    vel_wall=0.d0
  end if

return
end subroutine init_filaments


