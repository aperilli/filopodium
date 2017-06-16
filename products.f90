!prodotto matrici 3x3 e vettori di dim.3
module products
implicit none 

contains

subroutine prod_matrix(A,B,C)
implicit none

real(kind=8), intent(in) :: A(3,3),B(3,3)
real(kind=8), intent(out) :: C(3,3)

integer :: i,j,k

C=0.d0
do i=1,3
  do j=1,3
    do k=1,3
      C(i,j)=C(i,j)+A(i,k)*B(k,j)
    end do
  end do
end do
return
end subroutine prod_matrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine matrix_vector(A,v,B)
implicit none

real(kind=8), intent(in) :: A(3,3),v(3)
real(kind=8), intent(out) :: B(3)

integer :: i,j

B=0.d0
do i=1,3
  do j=1,3
    B(i)=B(i)+A(i,j)*v(j)
  end do
end do

return
end subroutine matrix_vector

end module products
