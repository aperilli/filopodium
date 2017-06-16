subroutine integrate_filament(time,Lwall,nmon,pos,vel)

use common_variables
use st_mpi
use mpi_routines

  implicit none

  integer, dimension(nfil), intent(in) :: nmon
  real(kind=8), intent(in) :: time,Lwall
  real(kind=8), dimension(3,max_size_fil,nfil), intent(inout) :: pos,vel

  real(kind=8), dimension(3,max_size_fil,nfil) :: force,ksi,eeta
  real(kind=8), dimension(12) :: rgauss
  real(kind=8) :: ekin, energy,fwall
  integer :: i,j,k,m

  call compute_forces(Lwall,pos,nmon,force,fwall) 

  do j=1,nfil
    m=nmon(j)
    do i=3,m
      call randomgaussian(rgauss)
      do k=1,3
        ksi(k,i,j)=rgauss(2*k-1)
        eeta(k,i,j)=rgauss(2*k)
      end do
    end do
  end do
  call integrate_velocity(ksi,eeta,force,nmon,vel)
  call integrate_position(Lwall,eeta,vel,nmon,pos)
  do j=1,nfil
    m=nmon(j)
    do i=3,m
      if(pos(1,i,j).lt.0.d0)then
        print*,myrank+1,fwall,force(1,i,j),vel(1,i,j),pos(1,i,j),eeta(1,i,j),ksi(1,i,j)
        call parallel_abort()
      end if
    end do
  end do
  call compute_forces(Lwall,pos,nmon,force,fwall)
  call integrate_velocity(ksi,eeta,force,nmon,vel)

return

end subroutine integrate_filament

