program filopod_wrapper

use st_mpi
use mpi_routines

implicit none

  integer :: mngrp=800,mnprm=20
  character(len=30) :: fs,qid_global
  character(len=30), allocatable :: qid(:),p(:)
  integer :: lpx,i,j,n,cnt,procsum,tmp1,cnt1,split_key
  integer, allocatable :: pnum(:)
  integer :: rng_pack_length,ipickoff,intread

  real(kind=8) :: rng,rng_cm
  external :: rng,rng_cm
  external :: pack_rng_stream,pack_rng_cm_stream,get_rng_pack_length
  external :: ipickoff,intread
  logical :: tmp

  call parallel_init()
  allocate(qid(nproc_world),p(mnprm),pnum(nproc_world))
    
!  if (nproc_world.gt.maxproc_world) then
!    if (root_world) write (*,*) 'nproc_world.gt.maxproc_world' &
!      ,maxproc_world
!    call parallel_abort()
!  endif     questo non dovrebbe servire

  qid_global='run'
  lpx=index(qid_global,' ')-1
  fs=qid_global(1:lpx)//'.in'
  if (root_world) write (*,*) ' trying to open file ',fs
  open(1,file=fs,status='old',form='formatted')

  procsum=0
  cnt=0
  tmp=.false.
901   j=ipickoff(1,p,n,mnprm,tmp)
  print *,'myrank_world',myrank_world,j,n,p(1:n)
  if(j.ne.0) go to 902
  if(n.le.0) go to 901

  cnt=cnt+1
  if(n.lt.2) then
    if(root_world) write(*,*)'insufficient parameters in input file'
    call parallel_abort()
  endif

!  if(cnt.gt.mngrp) then !too many runs 
!  if (root_world) write (*,*) 'Too many runs:',mngrp
!       call parallel_abort()
!      endif

  qid(cnt)=p(1)
  pnum(cnt)=intread(p(2))

  if(pnum(cnt).gt.maxproc) then !too many procs in one run 
    if (root_world) write (*,*) 'Too many procs in a single run:' &
      ,pnum(cnt),maxproc
       call parallel_abort()
  endif

  procsum=procsum+pnum(cnt)

  go to 901

! finished reading file,
! make sure nproc_world is equal to the sum of processors in each job
902 if(procsum.ne.nproc_world) then
    if (root_world) write (*,*) 'Number of procs does not agree with &
        number in input file' ,procsum,nproc_world
    call parallel_abort()
  endif

! assign keys 
  split_key=0
  cnt1=1
  tmp1=0
903 if(cnt1.gt.cnt) then
    if (root_world) write (*,*) 'Problem assigning keys.'
    call parallel_abort()
  endif
  tmp1=tmp1+pnum(cnt1)
  if(myrank_world.lt.tmp1) then
    split_key=cnt1
    go to 904
  else
    cnt1=cnt1+1
    go to 903
  endif

904 call split_comm(split_key)   !creates local communicators

  kill_early=.false.   !fixed for now

  call filopodium_mpi(qid(split_key))   ! call filopodium_mpi with corresponding qid

  if(kill_early) then
    call parallel_abort()
  else
    call parallel_end()
  endif


end program filopod_wrapper
