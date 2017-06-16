subroutine randomgaussian(rgauss)

use common_variables, only : pi

  implicit none
  real(kind=8), dimension(12), intent(out) :: rgauss
  real(kind=8) :: rng
  external :: rng
  integer :: i

!  call random_number(rgauss)
  do i=1,12
    rgauss(i)=rng()
  end do
  do i =1,6
    rgauss(i)=sqrt(-2.0*log(rgauss(2*i-1))) * cos(2*pi*rgauss(2*i))
  end do

  return

end subroutine randomgaussian



