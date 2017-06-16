function intread(p)
implicit none

  integer :: intread,ln
  character(len=*) :: p!(*)
  character(len=25) :: blanks,both
  blanks='                '
  ln=index(p,' ')-1
  both=blanks(1:25-ln)//p(1:ln)
  read(both,'(i25)') intread
  print*,'both,ln',both,ln

return
end

