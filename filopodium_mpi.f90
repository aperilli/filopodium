subroutine filopodium_mpi(qid)

use common_variables
use compute_energy
use st_mpi
use mpi_routines

  implicit none
  integer :: size_fil,ntraj
  integer :: i,j,itime,k,m,nhist,jp,ln,cont,ico
  
  real(kind=8), allocatable :: pos(:,:,:),posj(:,:)
  real(kind=8), allocatable :: vel(:,:,:)
  real(kind=8), allocatable :: force(:,:,:)
  real(kind=8) :: energy,time,epot,ekin,en_bond,en_bend,en_intra,en_wall,en_i
  real(kind=8), allocatable :: ksi(:,:,:),eeta(:,:,:),t_in(:),rperp(:)
  real(kind=8), allocatable :: pk(:),pk_all(:)
  real(kind=8), dimension(3) :: pos_new,vel_new
  real(kind=8) :: a1,cnorm,ymin,ymax,zmin,zmax,xmin,s_pk
  
  integer:: start_time_array(8),curr_time_array(8)
  real*8 :: start_time,curr_time,elapsed_time
 
  real(kind=8) :: fwall,vel_wall,Lwall,posjx

  real(kind=8), dimension(3) :: pos_send,vel_send

  real(kind=8), dimension(1) :: time_array,Lwall_array

  integer :: values(1:8), kran, size_vmd, delta_n,n
  integer, dimension(:), allocatable :: seed,nmon,n_old
!  integer, allocatable :: linked_list(:),head(:,:,:)
  integer :: nattempt_p,nattempt_d,pol,depol,npol,ndepol,ncellx,ncelly,ncellz

  integer :: idebug
  logical, parameter :: reseed=.false.,debug=.false.,verbose=.false.
  logical :: ife

  character(len=30) :: qid,qid_local
  character(len=70) :: fchkrng

  integer :: rng_pack_length=40000
  integer :: rng_size,get_rng_pack_length,rnk
  integer, dimension(1) :: rng_size_array
  real(kind=8) :: rng,rng_cm
  external :: rng,rng_cm
  external :: pack_rng_stream,pack_rng_cm_stream,get_rng_pack_length
  external :: unpack_rng_cm_stream,unpack_rng_stream
  character, allocatable :: store_rng(:)

  real(kind=8) :: boxx,boxy,boxz
  character(len=2) :: nfs
  character(len=5) :: nposvel

  call date_and_time(values=start_time_array)
  start_time=start_time_array(5)*3600 &
            +start_time_array(6)*60 &
            +start_time_array(7) &
            +start_time_array(8)*0.001

  allocate(store_rng(rng_pack_length))

  iseed=777+111*myrank
  idebug=0
  if (root) print *,'comm_id =',comm_id
  call easy_init_rng(iseed,myrank,nproc,idebug)
  write(*,*)'iseed: ',myrank,myrank_world,iseed,rng()

  if(reseed) then
    read(*,*) iseed
    call easy_init_rng(iseed,myrank,nproc,idebug)
    if (debug) write(3770+myrank,*)'first Random#:',rng()          ! this is local to the specific core
    if (debug) write(3770+myrank,*)'first CM_Random#:',rng_cm()    ! this is global to the pool of cores
  endif
  
  call read_input(qid,ntraj)
  allocate(t_in(ntraj))

  if(root) then
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
  end if

  write(nfs,'(i2.2)') nfil

  allocate(pos(3,max_size_fil,nfil))
  allocate(posj(3,max_size_fil))
  allocate(vel(3,max_size_fil,nfil))
  allocate(force(3,max_size_fil,nfil))
  allocate(ksi(3,max_size_fil,nfil),eeta(3,max_size_fil,nfil))
  allocate(nmon(nfil),n_old(nfil))
  allocate(rperp(nfil))

  nhist=250!int(box_size(1)/dl0)+20
  allocate(pk(-nhist:nhist),pk_all(-nhist:nhist))
  pk=0.d0

  fchkrng=TRIM(qid)//'.rng'
  inquire(file=TRIM(fchkrng),exist=ife)
  if(ife.eqv..true.) then
    open(37,file=TRIM(fchkrng),form='unformatted',status='old')
    read(37) rng_size
    rewind(37)
    read(37) rng_size,store_rng(1:rng_size)
    call unpack_rng_cm_stream(store_rng)

111 read(37,end=222) rnk,rng_size,store_rng(1:rng_size)
    if(rnk.ne.myrank) goto 111
222 call unpack_rng_stream(store_rng)
    if(debug.eqv..true.) write(3770+myrank,*) 'local rng initialized',rng_size,rng()
    close(37)
  else
    if (root) write(6,*) 'rng file not found'
    if (itype==1) call parallel_abort()
  endif


  npol=0
  ndepol=0
  nattempt_p=0
  nattempt_d=0
  pol=0
  depol=0

  call init_filaments(qid,nmon,pos,vel,Lwall,vel_wall,t_in)

  epot=0.d0
  energy=0.d0
  do j=1,nfil
    do i=2,nmon(j),2
      call bonding_energy(pos(1,1,j),nmon(j),i,en_bond)
      epot=epot+en_bond
    end do
    do i=3,nmon(j),3
      call bending_energy(pos(1,1,j),nmon(j),i,en_bend)
      epot=epot+en_bend
    end do
    do i=1,nmon(j)
      call intra_energy(pos(1,1,j),nmon(j),i,en_intra)
      call filwall_energy(Lwall,pos(1,i,j),en_wall)
      epot=epot+0.5d0*en_intra+en_wall
    end do
  end do
  call compute_ekin(vel,nmon,ekin)
  energy=ekin+epot+0.5d0*(mass_wall*vel_wall**2+ktrap*Lwall**2)
  if(root)print*,'initial energy: ',energy

  ln=index(qid,' ')
  qid_local=qid
  write(qid_local(ln:ln+2),'(i3.3)')myrank+1
  if(itype.ne.3)open(70+myrank,file=TRIM(qid_local)//'.out',form='formatted',status='unknown')
  if(itype.eq.3)open(70+myrank,file=TRIM(qid_local)//'.out',form='formatted',access='append',status='unknown')
  
  ico=1
  if(itype.eq.3)ico=floor(t_in(myrank+1)/freq_out)+1

  do itime=1,nstep
    time=itime*timestep+t_in(myrank+1)
    n_old(:)=nmon(:)
    call attempt_reaction(Lwall,nmon,pos,vel,pol,depol)
    nattempt_p=nattempt_p+pol
    nattempt_d=nattempt_d+depol
    delta_n=sum(nmon(:)-n_old(:))
    if(delta_n.ne.0.d0) then
      if(delta_n.eq.-1) then
        ndepol=ndepol+1
      elseif(delta_n.eq.1) then 
        npol=npol+1
      else
        print*,myrank+1,'something went wrong: delta_n', delta_n
        call parallel_abort() 
      end if
    end if
    do j=1,nfil
      m=nmon(j)-int(Lwall/dl0)-1
      pk(m)=pk(m)+1.d0
    end do
    call integrate_filament(time,Lwall,nmon,pos,vel)
    if(itype.eq.2.or.itype.eq.3) call integrate_wall(time,nmon,pos,Lwall,vel_wall)
    if(time.ge.(ico*freq_out))then
      do j=1,nfil
        rperp(j)=sqrt((pos(2,nmon(j),j)-pos(2,1,j))**2+&
                      (pos(3,nmon(j),j)-pos(3,1,j))**2)
      enddo
      ico=ico+1
      write(70+myrank,'(2g20.10,'//nfs//'i7,4'//nfs//'g20.10)') &
                time,Lwall,nmon(:),(pos(1:3,nmon(j),j),j=1,nfil),(rperp(j),j=1,nfil)
      call date_and_time(values=curr_time_array)
      curr_time=curr_time_array(5)*3600 &
               +curr_time_array(6)*60 &
               +curr_time_array(7) &
               +curr_time_array(8)*0.001
      elapsed_time=curr_time-start_time
      if(start_time_array(3).ne.curr_time_array(3))elapsed_time=24*3600-start_time+curr_time
      if(elapsed_time.ge.81600)goto 3000
    end if
  end do
  
3000 continue

  close(70+myrank)

  print*, myrank+1,'time loop finished'

  call sum_reduce_array_comm(pk,pk_all,nhist,my_comm_local)
  if (root) then
    s_pk=sum(pk_all)
    open(3,file=TRIM(qid)//'.pk',form='formatted',status='unknown')
    do j=1,nhist
      if(pk_all(j).gt.0)write(3,'(i4,f10.6)')j,pk_all(j)/s_pk
    enddo
    close(3)
  endif

  if (root) then  ! receive and write the restart file
    open(4,file=TRIM(qid)//'.chk',form='formatted',status='unknown')
    write(nposvel,'(i5.5)')sum(nmon(:)*6)
    write(4,'(i4,2f14.4,'//nfs//'i4,'//nposvel//'f14.4)') &
      myrank+1,time,Lwall,nmon(1:nfil), &
      (((pos(i,j,k),i=1,3), j=1,nmon(k)),k=1,nfil), &
      (((vel(i,j,k),i=1,3), j=1,nmon(k)),k=1,nfil)
    do jp=1,nproc-1
      call recv_dble_id_scalar(Lwall,200+jp,jp)
      call recv_dble_id_scalar(time,1200+jp,jp)
      call recv_int_id(nmon,nfil,2200+jp,jp)
      print*,jp,nmon(:)
      cont=0
      do j=1,nfil
        do i=1,nmon(j)
          cont=cont+1
          call recv_dble_id(pos_send,3,3200+jp+cont,jp)
          cont=cont+1
          call recv_dble_id(vel_send,3,5200+jp+cont,jp)
          pos(:,i,j)=pos_send(:)
          vel(:,i,j)=vel_send(:)
        end do
        if(debug) write(3770+myrank,*)'jp=',myrank,'filament',j,'received'
      end do
      if(debug) write(3770+myrank,*) 'jp=',jp,'time=',time,' received'
      write(nposvel,'(i5.5)')sum(nmon(:)*6)
      print*,jp,nposvel,'receive'
      write(4,'(i6,2f12.4,'//nfs//'i4,'//nposvel//'f14.4)') &
        jp+1,time,Lwall,nmon(1:nfil), &
        (((pos(i,j,k),i=1,3), j=1,nmon(k)),k=1,nfil), &
        (((vel(i,j,k),i=1,3), j=1,nmon(k)),k=1,nfil)
    enddo
    close(4)
  else
    call send_dble_scalar(Lwall,0,200+myrank)
    call send_dble_scalar(time,0,1200+myrank)
    call send_int(nmon,nfil,0,2200+myrank)
    if(debug) write(3770+myrank,*) 'jp=',myrank,'time=',time,' send'
    cont=0
    do j=1,nfil
      do i=1,nmon(j)
        pos_send(:)=pos(:,i,j)
        vel_send(:)=vel(:,i,j)
        cont=cont+1
        call send_dble(pos(:,i,j),3,0,3200+myrank+cont)
        cont=cont+1
        call send_dble(vel(:,i,j),3,0,5200+myrank+cont)
      end do
      if(debug) write(3770+myrank,*)'jp=',myrank,'filament',j,'sent'
    end do
  endif


  if (root) then         
! only the master core perform this for the common string
    open(43,file=fchkrng,form='unformatted') ! checkpoint file
    rng_size = get_rng_pack_length()
    rng_size_array(1)=rng_size
    if(rng_size.gt.rng_pack_length) then
      write(6,*)'RNG_PACK_LENGTH not large enough.'
      write(6,*)'RNG_PACK_LENGTH,MINIMUM SIZE:',rng_pack_length,rng_size
      call flush(6)
      call parallel_abort()
    endif
    call pack_rng_cm_stream(rng_size,store_rng)
    write(43) rng_size,store_rng(1:rng_size)
  endif

! each core should executes this for the local string
! the way I have implemented in my code is that each core send the seed to the "root" (myrank==0) core

  call pack_rng_stream(rng_size,store_rng)
  if (root) then
! write its own local
    write(43) myrank,rng_size,store_rng(1:rng_size)
! receive from other cores
    do jp=1,nproc-1
      call recv_int_id(rng_size_array,1,5000+jp,jp)
      call recv_char_id(store_rng,rng_size,500+jp,jp)
      write(43) jp,rng_size,store_rng(1:rng_size)
    enddo
    close(43)
  else
! send to root 
    call send_int(rng_size_array,1,0,5000+myrank)
    call send_char(store_rng,rng_size,0,500+myrank)
  endif

  print*,'Trajectory: ', myrank+1
  print*,npol,' polymerizations over ',nattempt_p,' attempts'
  print*,ndepol,' depolymerizations over ',nattempt_d,' attempts'
  print*,nattempt_p+nattempt_d,' attempted reactions over ',nstep,' steps'

return

end subroutine filopodium_mpi

 
