module mpi_routines 
! various parallel routines  MPI

use st_mpi
include 'mpif.h'

contains

subroutine barrier_comm(comm)
implicit none
 
  integer :: ierr,comm

  call MPI_Barrier(comm,ierr)
      
end subroutine barrier_comm 

!---------------------------------------------------------------
!mmorales: subroutines for immediate communication

subroutine Isend_int(a,n,iproc,itag,handle)
implicit none

  integer :: n,iproc,itag,a(n),ierr,handle

  call MPI_Isend(a,n,MPI_INTEGER,iproc,itag &
                     ,my_comm_local,handle,ierr)
      
end subroutine Isend_int

!---------------------------------------------------------------

subroutine Wait_isend(request)
implicit none
       
  integer :: request,status(MPI_STATUS_SIZE),err

  call MPI_Wait(request,status,err)

end subroutine Wait_isend

!---------------------------------------------------------------

subroutine Test_mpi(handle,flag)
implicit none

  logical :: flag
  integer :: handle,err,status(MPI_STATUS_SIZE)

  call MPI_Test(handle,flag,status,err)

end subroutine Test_mpi

!---------------------------------------------------------------

subroutine Irecv_int(a,n,iproc,itag,handle)
implicit none

  integer :: n,itag,a(n),ierr,handle,iproc

  call MPI_Irecv(a,n,MPI_INTEGER,iproc,itag &
                ,my_comm_local,handle,ierr)
      
end subroutine Irecv_int


!mmorales

!---------------------------------------------------------------

subroutine parallel_init()
implicit none
         
  integer :: ierr,icolor,ikey

  call MPI_Init(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank_world,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc_world,ierr)
  root_world= .false.
  if (myrank_world.eq.0) root_world = .true.
!        group1 = .false.
!        group2 = .false.
!        if (myrank.lt.(nproc/2)) then
!              icolor = 1
!              ikey = myrank
!              group1 = .true.
!        endif
!        if (myrank.ge.(nproc/2)) then
!              icolor = 2
!              ikey = myrank - nproc/2
!              group2 = .true.
!        endif
!        call MPI_COMM_SPLIT(MPI_COMM_WORLD,icolor,ikey,my_comm,ierr)

end subroutine parallel_init

!---------------------------------------------------------------

subroutine split_comm(skey)    ! create communicators for separate runs
implicit none

  integer :: ierr,skey

  comm_world=MPI_COMM_WORLD
  call MPI_COMM_SPLIT(comm_world,skey,myrank_world,my_comm_local,ierr)

  call MPI_COMM_RANK(my_comm_local,myrank,ierr)
  call MPI_COMM_SIZE(my_comm_local,nproc,ierr)
  root= .false.
  if (myrank.eq.0) root = .true.

!         write(*,*) 'myrank_world,myrank,my_comm_local'
!         write(*,*) myrank_world,myrank,my_comm_local
end subroutine split_comm

!---------------------------------------------------------------

subroutine split_mycomm(skey,OldComm,OldRank,NewComm,NewRank,NewNprocs,Newroot)
! used to split local communicators even further 
! OldComm is a preexisting communication from which NewComm is created
implicit none

  integer :: ierr,skey,OldComm,NewComm,OldRank,NewRank
  integer :: NewNprocs
  logical :: Newroot

  call MPI_COMM_SPLIT(OldComm,skey,OldRank,NewComm,ierr)

  call MPI_COMM_RANK(NewComm,NewRank,ierr)
  call MPI_COMM_SIZE(NewComm,NewNprocs,ierr)
  Newroot= .false.
  if (NewRank.eq.0) Newroot = .true.

end subroutine split_mycomm

!---------------------------------------------------------------

subroutine parallel_end()
implicit none
       
  integer :: ierr
  call MPI_Finalize(ierr)
         
stop
end subroutine parallel_end

!---------------------------------------------------------------
!mmorales with the wrapper, the code should be aborted when problems are found

subroutine parallel_abort()
implicit none

  integer :: ierr,code

  if(root) call flush(6)
  if(root) call flush(0)
  call flush(myrank_world+70)

  call MPI_Abort(MPI_COMM_WORLD,1,ierr)
         
stop
end subroutine parallel_abort

!---------------------------------------------------------------

subroutine bcast_array_comm(a,n,myid,comm)
implicit none
         
  integer :: n,myid
  real(kind=8) :: a(n)
  integer :: ierr,comm
         
  call MPI_BCast(a,n,MPI_DOUBLE_PRECISION,myid,comm,ierr)
      
end subroutine bcast_array_comm

!---------------------------------------------------------------

subroutine bcast_array(a,n,myid)
implicit none

  integer :: n,myid
  real(kind=8) :: a(n)
  integer :: ierr

  call MPI_BCast(a,n,MPI_DOUBLE_PRECISION,myid,my_comm_local,ierr)
      
end subroutine bcast_array

!---------------------------------------------------------------

subroutine bcast_int_comm(int,n,myid,comm)
implicit none

  integer :: n,int,myid
  integer :: ierr,comm
       
  call MPI_BCast(int,n,MPI_INTEGER,myid,comm,ierr)
      
end subroutine bcast_int_comm

!---------------------------------------------------------------

subroutine bcast_int(integ,n,myid)
implicit none

  integer :: n,myid
  integer :: integ(n),ierr
       
  call MPI_BCast(integ,n,MPI_INTEGER,myid,my_comm_local,ierr)
      
end subroutine bcast_int

!---------------------------------------------------------------

subroutine bcast_char_comm(ch,n,myid,comm)
implicit none

  integer :: n,ierr,myid,comm
  character(len=20) :: ch
      
  call mpi_bcast(ch,20,mpi_character,myid,comm,ierr)
      
end subroutine bcast_char_comm

!---------------------------------------------------------------

subroutine bcast_char(ch,n,myid)
implicit none

  integer n,ierr,myid
  character(len=20) :: ch
      
  call mpi_bcast(ch,20,mpi_character,myid,my_comm_local,ierr)
      
end subroutine bcast_char

!---------------------------------------------------------------

subroutine sum_reduce(r)
implicit none

  integer :: ierr
  real(kind=8) :: r,temp
       
  temp =r
  call MPI_AllReduce(temp,r,1,MPI_DOUBLE_PRECISION,MPI_SUM &
                       ,my_comm_local,ierr)
      
end subroutine sum_reduce

!---------------------------------------------------------------

subroutine sum_reduce_comm(r,comm)
implicit none

  integer :: ierr,comm
  real(kind=8) :: r,temp
       
  temp =r
  call MPI_AllReduce(temp,r,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      
end subroutine sum_reduce_comm

!---------------------------------------------------------------

subroutine sum_reduce_int(i)
implicit none

  integer :: ierr
  integer :: i,temp
       
  temp =i
  call MPI_AllReduce(temp,i,1,MPI_INTEGER,MPI_SUM,my_comm_local,ierr)

end subroutine sum_reduce_int

!---------------------------------------------------------------

subroutine sum_reduce_int_comm(i,comm)
implicit none

  integer :: ierr
  integer :: i,temp,comm

  temp =i
  call MPI_AllReduce(temp,i,1,MPI_INTEGER,MPI_SUM,comm,ierr)

end subroutine sum_reduce_int_comm

!---------------------------------------------------------------

subroutine sum_reduce_array_comm(v,tmp,n,comm)
implicit none

  integer :: ierr,n,i,comm
  real(kind=8), intent(in) :: v(n)
  real(kind=8), intent(out):: tmp(n)
       
  do i=1,n
    tmp(i) =v(i)
  enddo
  call MPI_AllReduce(tmp,v,n,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      
end subroutine sum_reduce_array_comm

!---------------------------------------------------------------

subroutine sum_reduce_array(v,tmp,n)
implicit none

  integer :: ierr,n,i
  real(kind=8) :: v(n),tmp(n)
       
  do i=1,n
    tmp(i) =v(i)
  enddo
       
  call MPI_AllReduce(tmp,v,n,MPI_DOUBLE_PRECISION,MPI_SUM &
                        ,my_comm_local,ierr)
      
end subroutine sum_reduce_array

!---------------------------------------------------------------

subroutine sum_reduce_darray_comm(recv,send,n,comm)
implicit none

  integer :: ierr,n,i,comm
  real(kind=8) :: send(n),recv(n)
       
  call MPI_AllReduce(send,recv,n,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      
end subroutine sum_reduce_darray_comm

!---------------------------------------------------------------

subroutine sum_reduce_darray_comm_repl(array,n,comm)
implicit none

  integer :: ierr,n,i,comm
  real(kind=8) :: tmp(n),array(n)
       
  do i=1,n
    tmp(i) =array(i)
  enddo

  call MPI_AllReduce(tmp,array,n,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      
end subroutine sum_reduce_darray_comm_repl

!---------------------------------------------------------------

subroutine sum_reduce_int_array(v,tmp,n)
implicit none

  integer :: ierr,n,i
  integer :: v(n),tmp(n)
         
  do i=1,n        
    tmp(i) =v(i)
  enddo
         
  call MPI_AllReduce(tmp,v,n,MPI_INTEGER,MPI_SUM,my_comm_local,ierr)
      
end subroutine sum_reduce_int_array

!---------------------------------------------------------------
!  morales: additional communication routines

subroutine sendrecv_replace_dble(buf,cnt,dest,stag,src,rtag)
implicit none  
       
  integer :: cnt, dest, stag, src, rtag
  integer :: status(MPI_STATUS_SIZE), ierr
  real(kind=8) :: buf(cnt)

  call MPI_SENDRECV_REPLACE(buf,cnt,MPI_DOUBLE_PRECISION &
                ,dest,stag,src,rtag,my_comm_local,status,ierr)

end subroutine sendrecv_replace_dble

!---------------------------------------------------------------

subroutine sendrecv_replace_cmpl(buf,cnt,dest,stag,src,rtag)
implicit none

  integer :: cnt, dest, stag, src, rtag
  integer :: status(MPI_STATUS_SIZE), ierr
  complex :: buf(cnt)

  call MPI_SENDRECV_REPLACE(buf,cnt,MPI_DOUBLE_COMPLEX &
                 ,dest,stag,src,rtag,my_comm_local,status,ierr)

end subroutine sendrecv_replace_cmpl

!---------------------------------------------------------------

subroutine sendrecv_dble(sbuf,scnt,dest,stag,rbuf,rcnt,src,rtag,nrecv)
implicit none

  integer :: scnt,rcnt, dest, stag, src, rtag
  integer :: status(MPI_STATUS_SIZE), ierr,nrecv
  real(kind=8) :: sbuf(scnt),rbuf(rcnt)

  CALL MPI_SENDRECV(sbuf,scnt,MPI_DOUBLE_PRECISION,dest,stag &
                 ,rbuf,rcnt,MPI_DOUBLE_PRECISION,src,rtag, &
                  my_comm_local,status,ierr)
  CALL MPI_GET_COUNT(status,MPI_DOUBLE_PRECISION,nrecv,ierr)

end subroutine sendrecv_dble

!---------------------------------------------------------------

subroutine sendrecv_cmpl(sbuf,scnt,dest,stag,rbuf,rcnt,src,rtag)
implicit none

  integer :: scnt,rcnt, dest, stag, src, rtag
  integer :: status(MPI_STATUS_SIZE), ierr
  complex :: sbuf(scnt),rbuf(rcnt)

  call MPI_SENDRECV(sbuf,scnt,MPI_DOUBLE_COMPLEX &
                 ,dest,stag,rbuf,rcnt,MPI_DOUBLE_COMPLEX,src &
                 ,rtag,my_comm_local,status,ierr)

end subroutine sendrecv_cmpl

!---------------------------------------------------------------

subroutine send_char(a,n,iproc,itag)
implicit none

  integer :: n,iproc,itag
  character :: a(n)
  integer :: ierr

  call MPI_Send(a,n,MPI_CHARACTER,iproc,itag ,my_comm_local,ierr)

end subroutine send_char

!---------------------------------------------------------------

subroutine send_dble(a,n,iproc,itag)
implicit none

  integer :: n,iproc,itag
  real(kind=8) :: a(n)
  integer :: ierr

  call MPI_Send(a,n,MPI_DOUBLE_PRECISION,iproc,itag,my_comm_local,ierr)

end subroutine send_dble

!---------------------------------------------------------------

subroutine send_dble_scalar(a,iproc,itag)
implicit none

  integer :: iproc,itag
  real(kind=8) :: a
  integer :: ierr

  call MPI_Send(a,1,MPI_DOUBLE_PRECISION,iproc,itag,my_comm_local,ierr)

end subroutine send_dble_scalar


!---------------------------------------------------------------

subroutine send_cmp(a,n,iproc,itag)
implicit none

  integer :: n,iproc,itag
  complex :: a(n)
  integer :: ierr

  call MPI_Send(a,n,MPI_DOUBLE_COMPLEX,iproc,itag,my_comm_local,ierr)

end subroutine send_cmp

!---------------------------------------------------------------

subroutine send_int(a,n,iproc,itag)
implicit none

  integer :: n,iproc,itag,a(n),ierr

  call MPI_Send(a,n,MPI_INTEGER,iproc,itag,my_comm_local,ierr)

end subroutine send_int

!---------------------------------------------------------------

subroutine recv_dble(a,n,itag)
implicit none

  integer :: n,itag
  real(kind=8) :: a(n)
  integer :: ierr,status(MPI_STATUS_SIZE)

  call MPI_Recv(a,n,MPI_DOUBLE_PRECISION,0,itag,my_comm_local,status,ierr)

end subroutine recv_dble

!---------------------------------------------------------------

subroutine recv_int(a,n,itag)
implicit none

  integer :: n,itag,a(n),ierr,status(MPI_STATUS_SIZE)

  call MPI_Recv(a,n,MPI_INTEGER,0,itag,my_comm_local,status,ierr)

end subroutine recv_int

!---------------------------------------------------------------

subroutine recv_dble_any(a,n,itag,iproc)
implicit none

  integer :: n,itag,iproc
  real(kind=8) :: a(n)
  integer :: ierr,status(MPI_STATUS_SIZE)

  call MPI_Recv(a,n,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE &
                     ,itag,my_comm_local,status,ierr)
  iproc=status(MPI_SOURCE)
!  itag=status(MPI_TAG)

end subroutine recv_dble_any

!---------------------------------------------------------------

subroutine recv_int_any(a,n,itag,iproc)
implicit none

  integer :: n,itag,a(n),ierr,status(MPI_STATUS_SIZE),iproc

  call MPI_Recv(a,n,MPI_INTEGER,MPI_ANY_SOURCE,itag,my_comm_local,status,ierr)
  iproc=status(MPI_SOURCE)
!  itag=status(MPI_TAG)

end subroutine recv_int_any

!---------------------------------------------------------------

subroutine recv_dble_id(a,n,itag,iproc)
implicit none

  integer :: n,itag,iproc,ip
  real(kind=8) :: a(n)
  integer :: ierr,status(MPI_STATUS_SIZE)

  call MPI_Recv(a,n,MPI_DOUBLE_PRECISION,iproc &
                    ,itag,my_comm_local,status,ierr)
  ip=status(MPI_SOURCE)
  if (ip.ne.iproc)write (6,*) 'recv_dble_id error: iproc,ip',iproc,ip
!  itag=status(MPI_TAG)

end subroutine recv_dble_id

!---------------------------------------------------------------

subroutine recv_dble_id_scalar(a,itag,iproc)
implicit none
    
  integer :: itag,iproc,ip
  real(kind=8) :: a
  integer :: ierr,status(MPI_STATUS_SIZE)
  
  call MPI_Recv(a,1,MPI_DOUBLE_PRECISION,iproc &
                    ,itag,my_comm_local,status,ierr)
  ip=status(MPI_SOURCE)
  if (ip.ne.iproc)write (6,*) 'recv_dble_id error: iproc,ip',iproc,ip
!  itag=status(MPI_TAG)
      
end subroutine recv_dble_id_scalar

!---------------------------------------------------------------

subroutine recv_char_id(a,n,itag,iproc)
implicit none

  integer :: n,itag,ierr,status(MPI_STATUS_SIZE),iproc
  character :: a(n)

  call MPI_Recv(a,n,MPI_CHARACTER,iproc &
                     ,itag,my_comm_local,status,ierr)
!  iproc=status(MPI_SOURCE)
!  itag=status(MPI_TAG)

end subroutine recv_char_id

!---------------------------------------------------------------

subroutine recv_cmp_id(a,n,itag,iproc)
implicit none

  integer :: n,itag,iproc,ip
  complex :: a(n)
  integer ierr,status(MPI_STATUS_SIZE)

  call MPI_Recv(a,n,MPI_DOUBLE_COMPLEX,iproc &
                     ,itag,my_comm_local,status,ierr)
  ip=status(MPI_SOURCE)
  if (ip.ne.iproc) write (6,*) 'recv_dble_id error: iproc,ip',iproc,ip
!  itag=status(MPI_TAG)

end subroutine recv_cmp_id

!---------------------------------------------------------------

subroutine recv_int_id(a,n,itag,iproc)
implicit none

  integer :: n,itag,a(n),ierr,status(MPI_STATUS_SIZE),iproc

  call MPI_Recv(a,n,MPI_INTEGER,iproc &
                 ,itag,my_comm_local,status,ierr)
!  iproc=status(MPI_SOURCE)
!  itag=status(MPI_TAG)

end subroutine recv_int_id

!---------------------------------------------------------------

subroutine second(ttt)
implicit none

  real(kind=8) :: ttt

  ttt=MPI_WTIME()

end subroutine second

end module
