subroutine integrate_velocity(ksi,eeta,force,nmon,vel)

use common_variables, only : nfil,fric_mon,mass,timestep,kBT,max_size_fil

  implicit none

  integer, dimension(nfil), intent(in) :: nmon
  real(kind=8), dimension(3,max_size_fil,nfil), intent(in) :: ksi, eeta
  real(kind=8), dimension(3,max_size_fil,nfil), intent(in) :: force
  real(kind=8), dimension(3,max_size_fil,nfil), intent(inout) :: vel

  real(kind=8) :: f(3),sq3,h2g,h2,hg,h32sg,h12s,h32s,sigma,gammaa
  integer :: i,j

  gammaa=fric_mon!/mass
  sigma=sqrt(2.d0*kBT*gammaa/mass)
  sq3=sqrt(3.d0)
  h2=0.5d0*timestep
  h2g=0.5d0*timestep**2*gammaa
  hg=0.5d0*timestep*gammaa
  h32sg=0.25d0*sqrt(timestep**3)*sigma*gammaa
  h12s=0.5d0*sqrt(timestep)*sigma

  do j=1,nfil
    do i=3,nmon(j)
      f(:)=force(:,i,j)/mass
      vel(:,i,j)=vel(:,i,j)+h2*f(:)-hg*vel(:,i,j) &
                 +h12s*ksi(:,i,j)-0.25d0*h2g*(f(:)-gammaa*vel(:,i,j)) &
                 -h32sg*(0.5d0*ksi(:,i,j)+eeta(:,i,j)/sq3)
    end do
  end do

return

end subroutine integrate_velocity


