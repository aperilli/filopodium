subroutine read_input(qid,ntraj)

use common_variables
use st_mpi
use mpi_routines

  implicit none
  character(len=30), intent(in) :: qid
  integer, intent(out) :: ntraj
 
  real(kind=8) :: lrc,lrc6,lcw,lcw3
  character(len=256) :: line,keyword,keyword1,inputfile
  logical :: ifex,foundsharp
  integer :: iostat
  integer :: i,j,k

  inquire(file='inputfile.in',exist=ifex)
  if(ifex) then
    open(10,file='inputfile.in',status='old',form='formatted')
  else
    write(*,*) 'input file not found.'
    call parallel_abort() 
  endif
 
  do k=1,100
    read(10,'(a)',iostat=iostat) line
! when the file finishes, exit the loop
    if(iostat/=0) exit
! delete everything past an eventual "#" comment
    foundsharp=.false.
    do i=1,len(line)
      if(line(i:i)=="#") foundsharp=.true.
      if(foundsharp) line(i:i)=" "
    end do
! if the remaining line is empty, skip it
    if(len_trim(line)==0) cycle
! read the first word from line
    read(line,*) keyword
! the second word is then read to the proper variable
    select case(keyword)
    case("nfil")
      read(line,*) keyword1,nfil
    case("timestep")
      read(line,*) keyword1,timestep
    case("box_size")
      read(line,*) keyword1,box_size(1),box_size(2),box_size(3)
    case("eps_intra")
      read(line,*) keyword1,eps_intra
    case("ktrap")
      read(line,*) keyword1,ktrap
    case("mass_wall")
      read(line,*) keyword1,mass_wall
    case("nstep")
      read(line,*) keyword1,nstep
    case("k_bond")
      read(line,*) keyword1,k_bond
    case("lp")
      read(line,*) keyword1,lp
    case("fric_mon")
      read(line,*) keyword1,fric_mon
    case("itype")
      read(line,*) keyword1,itype
    case("rcut_intra")
      read(line,*) keyword1,rcut_intra
    case("sigma_intra")
      read(line,*) keyword1,sigma_intra
    case("eps_wall")
      read(line,*) keyword1,eps_wall
    case("sigma_wall")
      read(line,*) keyword1,sigma_wall
    case("rho1")
      read(line,*) keyword1,rho1
    case("mass")
      read(line,*) keyword1,mass
    case("eps_0")
      read(line,*) keyword1,eps_0
    case("W0")
      read(line,*) keyword1,W0
    case("rcut_inter")
      read(line,*) keyword1,rcut_inter
    case("sigma_inter")
      read(line,*) keyword1,sigma_inter
    case("eps_inter")
      read(line,*) keyword1,eps_inter
    case("delta_inter")
      read(line,*) keyword1,delta_inter
    case("ntraj")
      read(line,*) keyword1,ntraj
    case("freq_out")
      read(line,*) keyword1,freq_out
    case default
! an unknown word will stop the execution
      if(root)write(0,*) "Unknown keyword :",trim(keyword)
      call parallel_abort()
    end select
  end do
  close(10)

  if(ntraj.ne.nproc) then
    if(root) write(0,*) 'run',qid,' : ntraj.ne.nproc',ntraj,nproc
    call parallel_abort()
  end if

  k_bend=kBT*lp/dl0
  fric_mon=fric_mon/mass
  fric_wall=kBT/D_wall
  fric_wall=fric_wall/mass_wall
  vmxd=sqrt(3.d0*kBT/mass)
  U0=W0*rho1
  
  rcut_wall=3.d0**(1.d0/6.d0)*sigma_wall

  lrc=sigma_intra/rcut_intra
  lrc6=lrc**6
  vcut_intra=4.d0*eps_intra*(lrc6**2-lrc6)

  lcw=sigma_wall/rcut_wall
  lcw3=lcw**3
  vcut_wall=3.d0*sqrt(3.d0)/2.d0*eps_wall*(lcw3**3-lcw3)

  lrc=sigma_inter/rcut_inter
  lrc6=lrc**6
  vcut_inter=4.d0*eps_inter*(lrc6**2-lrc6)

return
end subroutine read_input

