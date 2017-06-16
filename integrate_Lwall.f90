subroutine integrate_Lwall(eeta,vel_wall,Lwall)

use common_variables

  implicit none

  real(kind=8), intent(in) :: eeta
  real(kind=8), intent(inout) :: vel_wall,Lwall

  real(kind=8) :: sq3,h32s,sigma,dL

  sigma=sqrt(2.d0*kBT*fric_wall/mass_wall)
  sq3=sqrt(3.d0)
  h32s=0.5d0*sqrt(timestep**3)*sigma

  dL=timestep*vel_wall+h32s*eeta/sq3
  Lwall=Lwall+dL

return

end subroutine integrate_Lwall
