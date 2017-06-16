subroutine compute_forces(Lwall,pos,nmon,force,fwall)

use common_variables
use st_mpi

  implicit none
  integer, dimension(nfil), intent(in) :: nmon
  real(kind=8), intent(in) :: Lwall 
  real(kind=8), dimension(3,max_size_fil,nfil), intent(in) :: pos
  real(kind=8), dimension(3,max_size_fil,nfil), intent(out) :: force
  real(kind=8), intent(out) :: fwall

  real(kind=8), dimension(3) :: dr,dq
  real(kind=8) :: r2,rr,q2,qq,lr,lr2,lr6,lw,lw3,distwall
  real(kind=8) :: add1,add2,add3,fr,costheta
  integer :: i,j,k,j3,m

  force(:,:,:)=0.d0
  fwall=0.d0

  do k=1,nfil
    m=nmon(k)
      !-- Bonding Potential Energy with each bond (1)
    do j=1,m-1
      dr(:)=pos(:,j+1,k)-pos(:,j,k)
 
      r2=sum(dr(:)**2)
      rr=sqrt(r2)
      fr=k_bond*(1.d0-dl0/rr)
      force(:,j  ,k)=force(:,j  ,k)+dr(:)*fr
      force(:,j+1,k)=force(:,j+1,k)-dr(:)*fr
    end do

    !-- Angular Potential Energy with each bend between two bonds (2)
    !-- 3-body
    do j=2,m-1
    !r(j-1,j)
      dr(:)=pos(:,j,k)-pos(:,j-1,k)
  
      r2=sum(dr(:)**2)!drx**2+dry**2+drz**2
      rr=sqrt(r2)

    !r(j,j+1)
      dq(:)=pos(:,j+1,k)-pos(:,j,k)
 
      q2=sum(dq(:)**2)
      qq=sqrt(q2)

      costheta=sum(dr(:)*dq(:))
      costheta=costheta/(rr*qq)

      do i=1,3
      add1=-k_bend*(dq(i)/(rr*qq)-costheta*(dr(i)/r2))
      force(i,j-1,k)=force(i,j-1,k)+add1
      add2=k_bend*((dq(i)-dr(i))/(rr*qq)-costheta*((dr(i)/r2)+(-dq(i)/q2)))
      force(i,j  ,k)=force(i,j ,k)+add2
      add3=k_bend*(dr(i)/(rr*qq)-costheta*dq(i)/q2)
      force(i,j+1,k)=force(i,j+1,k)+add3
      end do
    end do

    !--intramolecular LJ potential 
    do j=1,m-3
      j3=j+3
      do i=j3,m
        dr(:)=pos(:,j,k)-pos(:,i,k)

        r2=sum(dr(:)**2)
        rr=sqrt(r2)

        if(rr.lt.rcut_intra) then
          lr=sigma_intra/rr
          lr2=lr**2
          lr6=lr2**3
          fr=-24.d0*eps_intra*(lr6-2.d0*lr6**2)/r2
          force(:,j,k)=force(:,j,k)+dr(:)*fr
          force(:,i,k)=force(:,i,k)-dr(:)*fr
        end if
      end do
    end do

  !--interaction with the wall (LJ 9-3) -- force: x component only 

    do j=1,m 
      distwall=Lwall-pos(1,j,k)
      if(distwall.lt.rcut_wall) then 
        lw=sigma_wall/distwall
        lw3=lw**3
        fr=-9.d0*sqrt(3.d0)/2.d0*eps_wall*(3.d0*lw3**3-lw3)/distwall
        force(1,j,k)=force(1,j,k)+fr
        fwall=fwall-fr
      end if
    end do
  end do

  return
        
end subroutine compute_forces

