subroutine attempt_reaction(Lwall,nmon,pos,vel,nattempt_p,nattempt_d)

use common_variables
use st_mpi

implicit none

  real(kind=8), intent(in) :: Lwall
  real(kind=8), dimension(3,max_size_fil,nfil), intent(inout) :: pos,vel
  integer, dimension(nfil), intent(inout) :: nmon
  integer, intent(out) :: nattempt_p,nattempt_d
  real(kind=8) :: rng
  external :: rng

  real(kind=8), dimension(3) :: pos_nm1,pos_n,pos_new,vel_new
  real(kind=8), allocatable :: pos_reactive(:,:)
  real(kind=8) :: xi,psi
  integer :: reactive_fil,sizefil

  pos_new=0.d0
  vel_new=0.d0

  nattempt_p=0
  nattempt_d=0

  xi=rng()
  if(xi.ge.nfil*timestep*W0.and.xi.lt.nfil*timestep*(W0+U0)) then !polymerization
    psi=rng()
    reactive_fil=1+floor(psi*nfil)
    sizefil=nmon(reactive_fil)
    if(sizefil.lt.max_size_fil) then
      call polymerization(Lwall,pos(1,1,reactive_fil),sizefil,pos_new,vel_new)
      nattempt_p=nattempt_p+1
    endif
    if(sizefil.eq.nmon(reactive_fil)+1) then
      nmon(reactive_fil)=sizefil
      pos(:,sizefil,reactive_fil)=pos_new(:)
      vel(:,sizefil,reactive_fil)=vel_new(:)
    endif
  elseif(xi.lt.nfil*timestep*W0) then !depolymerization
    psi=rng()
    reactive_fil=1+floor(psi*nfil)
    sizefil=nmon(reactive_fil)
    if(sizefil.gt.2) then
      call depolymerization(Lwall,pos(1,1,reactive_fil),sizefil)
      nattempt_d=nattempt_d+1
    endif
    if(sizefil.eq.nmon(reactive_fil)-1) then
      nmon(reactive_fil)=sizefil
      pos(:,sizefil+1,reactive_fil)=0.d0
      vel(:,sizefil+1,reactive_fil)=0.d0
    end if
!  elseif(xi(1).lt.timestep*U0.and.xi(2).lt.timestep*W0) then
!    print*,'double attempted de/pol, traj.', myrank+1
!    n_double=1
!    if(xi(1).lt.xi(2)) then !polymerization
!      psi=rng()
!      reactive_fil=1+floor(psi*nfil)
!      sizefil=nmon(reactive_fil)
!      if(sizefil.lt.max_size_fil) then
!        call polymerization(Lwall,pos(1,1,reactive_fil),sizefil,pos_new,vel_new)
!        nattempt_p=nattempt_p+1
!      endif
!      if(sizefil.eq.nmon(reactive_fil)+1) then
!        nmon(reactive_fil)=sizefil
!        pos(:,sizefil,reactive_fil)=pos_new(:)
!        vel(:,sizefil,reactive_fil)=vel_new(:)
!      endif
!    elseif(xi(2).lt.xi(1)) then !depolymerization
!      psi=rng()
!      reactive_fil=1+floor(psi*nfil)
!      sizefil=nmon(reactive_fil)
!      if(sizefil.gt.2) then
!        call depolymerization(Lwall,pos(1,1,reactive_fil),sizefil)
!        nattempt_d=nattempt_d+1
!      endif
!      if(sizefil.eq.nmon(reactive_fil)-1) then
!        nmon(reactive_fil)=sizefil
!        pos(:,sizefil+1,reactive_fil)=0.d0
!        vel(:,sizefil+1,reactive_fil)=0.d0
!      end if
!    endif
  endif
  

return

end subroutine attempt_reaction
