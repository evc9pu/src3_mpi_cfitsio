subroutine AVcontour(fname,Av)

  use grid_mod
  use dust_mod
  use tts_mod
  use opacin_mod
  use constants

  implicit none

  real(8) :: thet,rad,taud,dr,Av

  integer :: ir,it,ip,id,iit,nthet
  character(len=30) :: fname

  open(unit=15,file='AVcontour_'//trim(fname)//'.dat',status='unknown')
  nthet=int(180.d0/0.25d0)
  do iit=1,nthet-1
     thet=pi/180.d0*dble(iit)*0.25d0
     call locate(thetarr,ntg,thet,it)
     if (thetarr(it).gt.thet) it=it-1
     taud=0.d0
     ip=1
     ir=1
     do while (taud.lt.Av/1.086d0.and.ir.lt.nrg-1)
        ir=ir+1
        dr=rarr(ir+1)-rarr(ir)
        rad=0.5d0*(rarr(ir+1)+rarr(ir))
        do id=1,4
           if (densarr(ir,it,ip,id).gt.0.d0) then
              taud=taud+kapd(findopac(ir,it,ip,id))*densarr(ir,it,ip,id)*dr
           end if
        end do
     end do
     write(15,*) rad/autors,thet
  end do
  close(15)

  return

end subroutine AVcontour
