subroutine cumul(a,b,c,n)
implicit none

  integer :: n,k,i
  real(kind=8), intent(in) :: a(1:n)
  real(kind=8), intent(inout) :: b(1:2*n),c(1:2*n)

  k=0
  do i=1,n
    k=k+2
    b(k)=(b(k)+b(k-1)**2)*c(k)+a(i)**2
    b(k-1)=b(k-1)*c(k-1)+a(i)
!      if(a(i)/=0.d0) then 
    c(k-1)=c(k-1)+1.d0
    c(k)=c(k)+1.d0
!      endif
    b(k-1)=b(k-1)/max(1.d0,c(k-1))
    b(k)=b(k)/max(1.d0,c(k))-(b(k-1)**2)
  enddo

return

end subroutine cumul
