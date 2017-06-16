subroutine integrate_position(Lwall,eeta,vel,nmon,pos)

use common_variables, only : timestep, fric_mon, mass, nfil,kBT,max_size_fil
use st_mpi
use mpi_routines

  implicit none

  integer, dimension(nfil), intent(in) :: nmon
  real(kind=8), dimension(3,max_size_fil,nfil), intent(in) :: eeta
  real(kind=8), intent(in) :: Lwall
  real(kind=8), dimension(3,max_size_fil,nfil), intent(inout) :: vel,pos

  real(kind=8) :: sq3,h32s,sigma,gammaa,dr(1:3)
  integer :: i,j

  gammaa=fric_mon!/mass
  sigma=sqrt(2.d0*kBT*gammaa/mass)
  sq3=sqrt(3.d0)
  h32s=0.5d0*sqrt(timestep**3)*sigma

  do j=1,nfil
    do i=3,nmon(j)
      dr(:)=timestep*vel(:,i,j)+h32s*eeta(:,i,j)/sq3
      !if(dr(1).gt.0.5d0)print*,myrank+1,'dx ',pos(1,i,j),dr(1),vel(1,i,j)
      !pos(:,i,j)=pos(:,i,j)+timestep*vel(:,i,j)+h32s*eeta(:,i,j)/sq3
      !if(pos(1,i,j)+dr(1).gt.Lwall) then
      !  print*,myrank+1,'collision',pos(1,i,j)+dr(1),'-->',Lwall-(pos(1,i,j)+dr(1)-Lwall)
      !  pos(1,i,j)=Lwall-(pos(1,i,j)+dr(1)-Lwall)
      !  vel(1,i,j)=-vel(1,i,j)
        !pos(1,i,j)=pos(1,i,j)+timestep*vel(1,i,j)+h32s*eeta(1,i,j)/sq3
        !print*,myrank+1,'overlap with the wall',pos(1,i,j),pos(1,i,j)-timestep*vel(1,i,j)-h32s*eeta(1,i,j)/sq3,vel(1,i,j)
        !call parallel_abort()
      !else  
        pos(1,i,j)=pos(1,i,j)+dr(1)
      !end if
      pos(2:3,i,j)=pos(2:3,i,j)+dr(2:3)
    end do
  end do

return

end subroutine integrate_position


