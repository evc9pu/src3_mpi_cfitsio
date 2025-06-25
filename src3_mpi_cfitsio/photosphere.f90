subroutine photosphere(fname)

  use grid_mod
  use dust_mod
  use tts_mod
  use opacin_mod

  implicit none

  real(8) :: thet,rad,taud,dr

  integer :: ir,it,ip,id
  character(len=30) :: fname

  open(unit=15,file='photosphere_'//trim(fname)//'.dat',status='unknown')
  do it=1,ntg-1
     thet=thetavearr(it)
     taud=0.d0
     ip=1
     ir=nrg-1
     do while (taud.lt.1.d0.and.ir.gt.1)
        ir=ir-1
        dr=rarr(ir+1)-rarr(ir)
        rad=0.5d0*(rarr(ir+1)+rarr(ir))
        do id=1,ndg+2
           if (.not.is_sg(id).and.densarr(ir,it,ip,id).gt.0.d0) then
              taud=taud+kapd(findopac(ir,it,ip,id))*densarr(ir,it,ip,id)*dr
!              if (it.eq.(ntg-1)/2) print*,ir,it,kapd(findopac(ir,it,ip,id))*densarr(ir,it,ip,id),dr
           end if
        end do
     end do
     write(15,*) rad/autors,thet
  end do
  close(15)

  return

end subroutine photosphere
