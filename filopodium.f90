program filopodium

use common_variables
use compute_energy

  implicit none
  integer :: size_fil
  integer :: i,j,itime,k,m,nhist
  
  real(kind=8), allocatable :: pos(:,:,:),posj(:,:)
  real(kind=8), allocatable :: vel(:,:,:)
  real(kind=8), allocatable :: force(:,:,:)
  real(kind=8) :: energy,time,epot,ekin,en_bond,en_bend,en_intra,en_wall,en_i
  real(kind=8), allocatable :: ksi(:,:,:),eeta(:,:,:)
  real(kind=8), allocatable :: j_hist(:),j_histav(:),j_histnorm(:),fwall(:,:)
  real(kind=8), dimension(3) :: pos_new,vel_new
  real(kind=8) :: a1,cnorm,ymin,ymax,zmin,zmax,xmin

  integer :: values(1:8), kran, size_vmd, delta_n,n
  integer, dimension(:), allocatable :: seed,nmon,n_old
!  integer, allocatable :: linked_list(:),head(:,:,:)
  integer :: nattempt_p,nattempt_d,pol,depol,npol,ndepol,ncellx,ncelly,ncellz

  real(kind=8) :: boxx,boxy,boxz
  character(len=2) :: nfs

  call read_input(size_fil)

  write(*,*) "Number of filaments      :",nfil
  write(*,*) "itype                    :",itype
  write(*,*) "Size of the box          :",box_size(1:3)
  write(*,*) "LJ: epsilon_intra        :",eps_intra
  write(*,*) "LJ: rcutoff_intra        :",rcut_intra
  write(*,*) "LJ: epsilon_inter        :",eps_inter
  write(*,*) "LJ: rcutoff_inter        :",rcut_inter
   write(*,*) "Persistence length       :",lp
  write(*,*) "Bonding pot. constant    :",k_bond
  write(*,*) "Bending pot. constant    :",k_bend
  write(*,*) "Monomer mass             :",mass
  write(*,*) "Timestep                 :",timestep
  write(*,*) "Number of steps          :",nstep
  write(*,*) "U0, W0                   :",U0,W0
  write(*,*) "LJ 9-3: epsilon_wall     :",eps_wall
  write(*,*) "LJ 9-3: rcutoff_wall     :",rcut_wall 

  write(nfs,'(i2.2)') nfil

  allocate(pos(3,max_size_fil,nfil))
  allocate(posj(3,max_size_fil))
  allocate(vel(3,max_size_fil,nfil))
  allocate(force(3,max_size_fil,nfil))
  allocate(ksi(3,max_size_fil,nfil),eeta(3,max_size_fil,nfil))
  allocate(fwall(max_size_fil,nfil))
  allocate(nmon(nfil))
  allocate(n_old(nfil))

  nhist=int(box_size(1)/dl0)+5
  allocate(j_hist(nhist),j_histav(2*nhist),j_histnorm(2*nhist))
  j_hist=0.d0
  j_histav=0.d0
  j_histnorm=0.d0

  call date_and_time(values=values)

  call random_seed(size=kran)
  allocate(seed(1:kran))
  seed(:) = values(8)
  call random_seed(put=seed)

  npol=0
  ndepol=0
  nattempt_p=0
  nattempt_d=0
  pol=0
  depol=0

  pos(:,:,:)=0.d0
  vel(:,:,:)=0.d0
  call init_filaments(nmon,pos,vel)

  open(20,file='test.xyz')

  size_vmd=int((box_size(1)+10.d0)/dl0)
  !write(20,*)nfil*size_vmd
  !write(20,*)
  do j=1,nfil
    m=nmon(j)
    posj(:,1:m)=pos(:,1:m,j)
    do i=1,m
      call bonding_energy(posj,m,i,en_bond)
      call bending_energy(posj,m,i,en_bend)
      call intra_energy(posj,m,i,en_intra)
      call filwall_energy(posj,m,i,en_wall)
  !    write(20,'(a4,3f15.8)') 'Ar',pos(:,i,j)
    end do
    epot=0.5d0*en_bond+en_bend/3.0d0+0.5d0*en_intra+en_wall-(m-1)*eps_0
  !  if(m.le.size_vmd) then
  !    do i=m+1,size_vmd
  !      write(20,'(a4,3f15.8)') 'Ar',0.d0,0.d0,0.d0
  !    end do
  !  end if
  end do
  write(20,*)
  write(20,*)

  call compute_ekin(vel,nmon,ekin)
!  call compute_forces(pos,nmon,force,fwall)
  energy=ekin+epot
  print*,'initial energy: ',energy

  do itime=1,nstep
    time=itime*timestep
    n_old(:)=nmon(:)
    call attempt_reaction(nmon,pos,vel,pol,depol)
    nattempt_p=nattempt_p+pol
    nattempt_d=nattempt_d+depol
    call integrate_filament(time,nmon,pos,vel)
    delta_n=sum(nmon(:)-n_old(:))
    if(delta_n.ne.0.d0) then
      write(200,'(g20.10,'//nfs//'i7)')time,nmon(:)
      if(delta_n.lt.0) ndepol=ndepol+1
      if(delta_n.gt.0) npol=npol+1
    endif
    do j=1,nfil
      m=nmon(j)
      j_hist(m)=j_hist(m)+1.d0
    end do
  !  if(mod(itime,500).eq.0)then
  !    do n=1,nfil
  !      m=nmon(n)
  !      do i=1,m
  !        write(20,'(a4,3f15.8)') 'Ar',pos(:,i,n)
  !      end do
  !      if(m.le.size_vmd)then
  !        do i=m+1,size_vmd
  !          write(20,'(a4,3f15.8)') 'Ar',0.d0,0.d0,0.d0
  !        end do
  !      end if
  !    end do
  !    write(20,*)
  !    write(20,*)
  !  end if
  end do

  open(4,file='restart.dat')
  do n=1,nfil
    write(4,*)nmon(n)
  end do
  close(4)


  j_hist(:)=j_hist(:)/sum(j_hist(:))
  call cumul(j_hist,j_histav,j_histnorm,nhist)
  open(4,file='jdistribution.dat')
  cnorm=max(j_histnorm(1)-1.d0,1.d0)
  k=0
  do j=1,nhist
    k=k+2
    a1=sqrt(j_histav(k)/cnorm)
    write(4,'(4g20.10)') j,j_hist(j),j_histav(k-1), a1
  end do
  close(4)

  print*,npol,' polymerizations over ',nattempt_p,' attempts'
  print*,ndepol,' depolymerizations over ',nattempt_d,' attempts'
  print*,nattempt_p+nattempt_d,' attempted reactions over ',nstep,' steps'

  close(20)




end program filopodium

 
