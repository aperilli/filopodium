module st_mpi
implicit none

integer :: nproc_world, myrank_world
integer :: my_comm_local,comm_world
integer :: comm_id,comm_id2
integer, parameter :: maxproc_world=512
     
integer :: nproc,myrank,iseed
integer, parameter ::maxproc=512
integer :: istatmpi
logical :: root,root_world,kill_early
      
common/mpi/nproc,myrank,iseed,nproc_world,myrank_world &
     ,istatmpi(maxproc),my_comm_local,comm_world &
     ,comm_id,comm_id2
common/lmpi/root,root_world,kill_early


end module st_mpi

