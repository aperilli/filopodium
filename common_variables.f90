module common_variables 
implicit none

real(kind=8) :: lp,vmxd!,kBT,pi,vmxd,mass
real(kind=8) :: dl0=1.d0
real(kind=8) :: kBT=1.d0
real(kind=8) :: pi=dacos(-1.d0)
!real(kind=8) :: mass=24.1e-20!fric_mon=1.d0
real(kind=8) :: D_wall=1.d0!fric_wall=1.d0
integer :: max_size_fil=200

real(kind=8) :: eps_intra,k_bond,k_bend,timestep,W0,freq_out
real(kind=8) :: eps_wall,rcut_wall,vcut_wall,sigma_wall
real(kind=8) :: rcut_intra,vcut_intra,sigma_intra,fric_wall
real(kind=8) :: eps_inter,rcut_inter,vcut_inter,sigma_inter
real(kind=8) :: rho1,U0,fric_mon,mass,eps_0,delta_inter,ktrap,mass_wall
real(kind=8), dimension(3) :: box_size 
integer :: nfil,itype
integer(kind=8) :: nstep


end module common_variables
