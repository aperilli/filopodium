module compute_energy

use common_variables
use mpi_routines
use st_mpi

contains

subroutine bonding_energy(pos,nmon,mon,en_bond)
! pos: position of each monomer in the filament
! nmon: number of monomers in the filament
! mon: monomer for which we calculate the energy
 
  implicit none
  real(kind=8), dimension (3,nmon), intent(in) :: pos
  integer, intent(in) :: nmon, mon 
  real(kind=8), intent(out) :: en_bond  

  integer :: i,j
  real(kind=8), dimension(3) :: dr
  real(kind=8) :: r2,rr
 
  en_bond=0.d0
  
  if(mon-1.ge.1) then
    dr(:)=pos(:,mon)-pos(:,mon-1)
    r2=sum(dr(:)**2)
    rr=sqrt(r2)
    en_bond=en_bond+0.5d0*k_bond*(rr-dl0)**2
  end if
  if(mon+1.le.nmon) then
    dr(:)=pos(:,mon+1)-pos(:,mon)
    r2=sum(dr(:)**2)
    rr=sqrt(r2)
    en_bond=en_bond+0.5d0*k_bond*(rr-dl0)**2
  end if

return 
end subroutine bonding_energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bending_energy(pos,nmon,mon,en_bend)

  implicit none
  real(kind=8), dimension (3,nmon), intent(in) :: pos
  integer, intent(in) :: nmon, mon 
  real(kind=8), intent(out) :: en_bend  

  real(kind=8), dimension(3) :: dr,dq
  real(kind=8) :: r2,rr,q2,qq,costheta
 
  en_bend=0.d0

  if(mon-2.ge.1) then
    dr(:)=pos(:,mon-1)-pos(:,mon-2)
    r2=sum(dr(:)**2)
    rr=sqrt(r2)
    dq(:)=pos(:,mon)-pos(:,mon-1)
    q2=sum(dq(:)**2)
    qq=sqrt(q2)

    costheta=sum(dr(:)*dq(:))
    costheta=costheta/(rr*qq)

    en_bend=en_bend+k_bend*(1-costheta)
  end if
  
  if(mon-1.ge.1.and.mon+1.le.nmon) then
    dr(:)=pos(:,mon+1)-pos(:,mon)
    r2=sum(dr(:)**2)
    rr=sqrt(r2)
    dq(:)=pos(:,mon)-pos(:,mon-1)
    q2=sum(dq(:)**2)
    qq=sqrt(q2)

    costheta=sum(dr(:)*dq(:))
    costheta=costheta/(rr*qq)

    en_bend=en_bend+k_bend*(1-costheta)
  end if

  if(mon+2.le.nmon) then
    dr(:)=pos(:,mon+2)-pos(:,mon+1)
    r2=sum(dr(:)**2)
    rr=sqrt(r2)
    dq(:)=pos(:,mon+1)-pos(:,mon)
    q2=sum(dq(:)**2)
    qq=sqrt(q2)

    costheta=sum(dr(:)*dq(:))
    costheta=costheta/(rr*qq)

    en_bend=en_bend+k_bend*(1-costheta)
  end if

return 
end subroutine bending_energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine intra_energy(pos,nmon,mon,en_intra)

  implicit none
  real(kind=8), dimension (3,nmon), intent(in) :: pos
  integer, intent(in) :: nmon, mon 
  real(kind=8), intent(out) :: en_intra
  integer :: i,j
  real(kind=8), dimension(3) :: dr
  real(kind=8) :: r2,rr,lr,lr2,lr6
 
  en_intra=0.d0

  if(mon-3.ge.1) then
    do i=1,mon-3
      dr(:)=pos(:,mon)-pos(:,i)
      r2=sum(dr(:)**2)
      rr=sqrt(r2)
      if(rr.eq.0.d0) then
        print*,'Self LJ interaction'
        call parallel_abort()
      end if        
      if(rr.lt.rcut_intra) then
        lr=sigma_intra/rr
        lr2=lr**2
        lr6=lr2**3
        en_intra=en_intra+4.d0*eps_intra*(lr6**2-lr6)-vcut_intra
      end if
    end do
  elseif(mon+3.le.nmon) then
    do i=mon+3,nmon
      dr(:)=pos(:,i)-pos(:,mon)
      r2=sum(dr(:)**2)
      rr=sqrt(r2)
      if(rr.eq.0.d0) then
        print*,'Self LJ interaction'
        call parallel_abort()
      end if        
      if(rr.lt.rcut_intra) then
        lr=sigma_intra/rr
        lr2=lr**2
        lr6=lr2**3
        en_intra=en_intra+4.d0*eps_intra*(lr6**2-lr6)-vcut_intra
      end if
    end do
  end if

return 

end subroutine intra_energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
subroutine filwall_energy(Lwall,pos,en_wall)

  implicit none
  real(kind=8), intent(in) :: Lwall
  real(kind=8), intent(in) :: pos
  real(kind=8), intent(out) :: en_wall  

  real(kind=8) :: lw,lw3,distwall
 
  en_wall=0.d0

  distwall=Lwall-pos
  if(distwall.lt.0.d0) then
    print*,myrank+1,'A filament overlaps with the wall!'
    call parallel_abort()
  end if
  if(distwall.lt.rcut_wall) then 
    lw=sigma_wall/distwall
    lw3=lw**3
    en_wall=en_wall+3.d0*sqrt(3.d0)/2.d0*eps_wall*(lw3**3-lw3)-vcut_wall
  end if
  
return
end subroutine filwall_energy


end module compute_energy
