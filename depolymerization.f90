subroutine depolymerization(Lwall,pos_reactive,sizefil)

use common_variables
use compute_energy
use st_mpi

implicit none

  integer, intent(inout) :: sizefil
  real(kind=8), intent(in) :: Lwall,pos_reactive(1:3,1:sizefil)
  real(kind=8) :: rng
  external :: rng

  real(kind=8) :: energy,acceptance,en_intra,en_wall,eta

  call intra_energy(pos_reactive,sizefil,sizefil,en_intra)
  call filwall_energy(Lwall,pos_reactive(1,sizefil),en_wall)
  acceptance=exp((en_intra+en_wall-eps_0)/kBT)
  eta=rng()
  if(eta.lt.acceptance) then
    sizefil=sizefil-1
  end if
 
return
end subroutine depolymerization
