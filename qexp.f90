function qexp(x)
implicit none
real(kind=16) :: qexp
real(kind=8), intent(in) :: x

qexp=exp(x)

end function
