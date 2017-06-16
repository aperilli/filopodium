function ipickoff(iu,p,n,mn,iroot)
! on entry: iu is unit number for read, mn is the size of array p
! on exit pickoff.eq.1 means eof was encountered
! p contains as characters the parameters, there are n of them
!  blanks or commas separate characters
! line length is maximally 160 and characters after a : , ! or ; are ignored.

use mpi_routines
use st_mpi

implicit none

  integer :: ipickoff,iu,n, mn,ifact,i,k,ld,ist,ln,ls
  integer :: line=160
  character(len=30) :: p(mn)
  character(len=161) :: datas
  logical, intent(in) :: iroot

100   continue
  ifact=0
  n=0 
  i=0
  read (iu,6,END=2) datas
  ipickoff=0
6  format(a160)
! starting from end find how long the line is
  do k=160,1,-1
    if(datas(k:k).ne.' ') go to 401
  enddo
401 ld=k
!     Additional safety check to stop XLF problem
!     If test above fails for all k (due to empty line) then 
!     ld will be initialized to 0, and the next write will crash.
  if(ld<1)return

  if(iu.le.1.and.iroot) write(6,7)'blabla',iu,datas(1:ld)
7  format(i3,' CMD:',a)
! now strip off comment part
  ls=0
  do k=1,ld
    if(datas(k:k).eq.':'.or.datas(k:k).eq.';'.or.datas(k:k).eq.'!') go to 402
    ls=k
  enddo
402  ls=ls+1
  datas(ls:ls)=' '

1  i=i+1
  if(i.gt.ls)then
    if(n.gt.0)then
!          if(root)write (6,*) ls,ld,n,'=n ',(p(i),i=1,n)
      return
    endif
! read a new record if a blank line is encountered
    go to 100
  endif

  if(datas(i:i).eq.' '.or.datas(i:i).eq.',')then
    if(ifact.ne.0) then
! finish the word
      ln=i-ist
      if(ln.ge.len(p(n))) then
        if(iroot)write (6,*) 'space for parameters',len(p(n)) &
                    ,' not large enough in pickoff',ln
        call parallel_abort()
      endif

      p(n)(1:ln)=datas(ist:i-1)
! terminate with  blanks
      do 30 k=ln+1,len(p(n))
30      p(n)(k:k)=' '
        ifact=0
    endif

  else

    if(ifact.eq.0) then
! found a new parameter
      ifact=1
      ist=i
      n=n+1
      if(n.gt.mn) then
        if(iroot)write (6,*)' too many parameters in pickoff',n,mn
        call parallel_abort()
      endif
    endif

  endif
  go to 1

2  ipickoff=1

return
end function ipickoff
