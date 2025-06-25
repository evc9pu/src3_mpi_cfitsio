subroutine waveset(cwav,cfldust,ctherm,nlam,MXLAM)
  !
  ! brilliant method for handling wavelength loop arrays
  !
  ! file created: 1999/11/20
  !
  ! history:
  !
  !
  implicit none

  integer :: MXLAM,nlam,i
  character(len=*) :: cwav(MXLAM),cfldust(MXLAM),ctherm(MXLAM)

  open(unit=55,file='wave.in',status='old')

  i=1
  do
     read(55,*,iostat=ioerr) cwav(i),cfldust(i),ctherm(i)
     if(ioerr.ne.0) exit
     i = i + 1
  end do
  nlam = i-1

  close(55)

  return
end subroutine waveset
