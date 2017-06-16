subroutine compute_ekin(vel,nmon,ek)

use common_variables, only : nfil,mass,max_size_fil

  implicit none
  integer, dimension(nfil), intent(in) :: nmon
  real(kind=8), dimension(3,max_size_fil,nfil), intent(in) :: vel
  real(kind=8) :: ek
  integer :: i,j

  do j=1,nfil
    do i=1,nmon(j)
      ek=0.5d0*mass*sum(vel(:,i,j)**2)
    end do
  end do

return

end subroutine compute_ekin

