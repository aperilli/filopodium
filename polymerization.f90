subroutine polymerization(Lwall,pos_reactive,sizefil,pos_new,vel_new)

use common_variables
use products
use compute_energy

implicit none

  real(kind=8), dimension(3), intent(out) :: pos_new,vel_new
  integer, intent(inout) :: sizefil
  real(kind=8), intent(in) :: pos_reactive(1:3,1:sizefil),Lwall 

  real(kind=8), dimension(3,sizefil+1) :: pos_trial
  real(kind=8), dimension(3) :: pos_nm1,pos_n
  real(kind=8) :: r,costheta,phi,rnd(2),theta,rlocal(3)
  real(kind=8), dimension(12) :: rgauss
  real(kind=8) :: th_nm1,phi_nm1,tnm1
  real(kind=8) :: energy_new,acceptance,en_bond,en_bend,en_intra,en_wall
  real(kind=8) :: dr(3),rr,r2,lr,lr2,lr6,distwall,lw,lw3
  real(kind=8) :: eta
  real(kind=8), dimension(3,3) :: rot_z,rot_y,rot_yz
  real(kind=8) :: rng
  external :: rng

  energy_new=0.d0
  pos_new=0.d0
  vel_new=0.d0

!  pos_nm1(:)=pos_reactive(:,sizefil-1)
!  pos_n(:)=pos_reactive(:,sizefil)

! r, theta, phi: bond length, bending angle and angle phi with respect to the 
! direction of the last bond

  call randomgaussian(rgauss)
  r=dl0+rgauss(1)*sqrt(kBT/2.d0/k_bond)

  costheta=1.d0+kBT*log(exp(-2.0d0*k_bend/kBT)+rng()*(1.0-exp(-k_bend/kBT)))/k_bend
  phi=2.d0*pi*rng()

  theta=acos(costheta)

! coordinates of the new monomer with respect to a local RS with the x axis along the last bond

  rlocal(1)=r*costheta
  rlocal(2)=r*sin(theta)*cos(phi)
  rlocal(3)=r*sin(theta)*sin(phi)

! tnm1: vector joining the last two monomers
! th_nm1 and phi_nm1: polar and azimuthal angle of tnm1 in the RS centered on the nth monomer

  tnm1=sum((pos_reactive(:,sizefil)-pos_reactive(:,sizefil-1))**2)
  tnm1=sqrt(tnm1)
  th_nm1=acos((pos_reactive(3,sizefil)-pos_reactive(3,sizefil-1))/tnm1)
  phi_nm1=atan((pos_reactive(2,sizefil)-pos_reactive(2,sizefil-1))/ &
               (pos_reactive(1,sizefil)-pos_reactive(1,sizefil-1)))

  rot_y(1,1)=sin(th_nm1)
  rot_y(2,1)=0.d0
  rot_y(3,1)=-cos(th_nm1)
  rot_y(1,2)=0.d0
  rot_y(2,2)=1.d0
  rot_y(3,2)=0.d0
  rot_y(1,3)=cos(th_nm1)
  rot_y(2,3)=0.d0
  rot_y(3,3)=sin(th_nm1) 

  rot_z(1,1)=cos(phi_nm1)
  rot_z(2,1)=sin(phi_nm1)
  rot_z(3,1)=0.d0
  rot_z(1,2)=-sin(phi_nm1)
  rot_z(2,2)=cos(phi_nm1)
  rot_z(3,2)=0.d0
  rot_z(1,3)=0.d0
  rot_z(2,3)=0.d0
  rot_z(3,3)=1.d0 

  call prod_matrix(rot_y,rot_z,rot_yz)
 ! position of the new monomer in the global RS

  call matrix_vector(rot_yz,rlocal,pos_new)
  pos_new(:)=pos_reactive(:,sizefil)+pos_new(:)
! acceptance test

  if(pos_new(1).gt.Lwall) then
    !print*,myrank+1, 'attempted polymerization -- overlap with the wall'
    return
  elseif(pos_new(1).lt.Lwall) then
    pos_trial(:,1:sizefil)=pos_reactive(:,:)
    pos_trial(:,sizefil+1)=pos_new(:)
    call intra_energy(pos_trial,sizefil+1,sizefil+1,en_intra)
    call filwall_energy(Lwall,pos_new(1),en_wall)
    acceptance=exp(-(en_intra+en_wall)/kBT)
    eta=rng()
    if(eta.lt.acceptance) then
      sizefil=sizefil+1
      if(sizefil.gt.max_size_fil) then
        print*,'Error: maximum filament size exceeded!'    
        call parallel_abort()
      end if
      call randomgaussian(rgauss)
      vel_new(1:3)=rgauss(1:3)*sqrt(kBT/mass)
    end if
  end if

return

end subroutine polymerization
