
module test

  implicit none
  
  integer :: id

  integer :: p(2,2)
  integer :: q(2,2)

  !$omp threadprivate(q)

end module test


program main

  use test
  use omp_lib

  implicit none

  integer :: num
  integer :: i,h

  q=1
  q=q*10

  p=0
  
  !$omp parallel private(id,h) reduction(+:p)
  !$omp do

  do i=1,4

     call add_test(i)

  end do

  !$omp end do
  !$omp end parallel

  stop

end program main

subroutine add_test(i)

  use test

  integer,save :: ii
  integer:: i

  print*,i
  if (i.eq.2) stop
!  print*,ran()

end subroutine add_test

  real(8) function ran() result(xi)
    implicit none
    real(8),save :: am
    integer, parameter :: ia=16807,im=2147483647,iq=127773,ir=2836
    integer, save :: ix=-1,iy=-1,k
    integer :: idum
    idum=-2211
    if (idum <= 0 .or. iy < 0) then
       am=nearest(1.0d0,-1.0d0)/im
       iy=ior(ieor(888889999,abs(idum)),1)
       ix=ieor(777755555,abs(idum))
       idum=abs(idum)+1
    end if
    ix=ieor(ix,ishft(ix,13))
    ix=ieor(ix,ishft(ix,-17))
    ix=ieor(ix,ishft(ix,5))
    k=iy/iq
    iy=ia*(iy-k*iq)-ir*k
    if (iy < 0) iy=iy+im
    xi=am*ior(iand(im,ieor(ix,iy)),1)
    return
  end function ran
