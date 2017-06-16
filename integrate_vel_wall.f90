subroutine integrate_vel_wall(ksi,eeta,fwall,vel_wall)

use common_variables

  implicit none

  real(kind=8), intent(in) :: ksi, eeta, fwall
  real(kind=8), intent(inout) :: vel_wall

  real(kind=8) :: sq3,h2g,h2,hg,h32sg,h12s,h32s,sigma

  sigma=sqrt(2.d0*kBT*fric_wall/mass_wall)
  sq3=sqrt(3.d0)
  h2=0.5d0*timestep
  h2g=0.5d0*timestep**2*fric_wall
  hg=0.5d0*timestep*fric_wall
  h32sg=0.25d0*sqrt(timestep**3)*sigma*fric_wall
  h12s=0.5d0*sqrt(timestep)*sigma

  vel_wall=vel_wall+h2*fwall/mass_wall-hg*vel_wall &
          +h12s*ksi-0.25d0*h2g*(fwall/mass_wall-fric_wall*vel_wall) &
          -h32sg*(0.5d0*ksi+eeta/sq3)

return

end subroutine integrate_vel_wall

