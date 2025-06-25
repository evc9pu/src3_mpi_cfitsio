subroutine photospherepeel(fname) 

  use grid_mod
  use dust_mod
  use tauint_mod
  use taunum_mod
  use constants
  use tts_mod
  use opacin_mod

  implicit none

  real(8) :: peelangle(4),peelthet,sint,cost,thet,tsum

  integer :: ndirec,idirec,ipeelthet,iphot
  integer :: ir,it,ip

  character(len=30) :: fname
  character(len=4) :: suffix

  peelangle=[0.001d0, 30.d0, 60.d0, 90.d0] ! in degree

  ndirec=10000
 
  do ipeelthet=1,4

     peelthet=peelangle(ipeelthet)*pi/180.d0

     write(suffix,'("_",I2.2,"_")') int(peelangle(ipeelthet))

     open(unit=99,file='photospherepeel'//suffix//trim(fname)//'.dat',status='unknown')

     do idirec=1,ndirec

        iphot=ndirec*(ipeelthet-1)+idirec
        iphot=1

        xp=0.d0
        thet=pi/dble(ndirec)*dble(idirec)-(pi/2.d0-peelthet)
        sint=sin(thet)
        cost=cos(thet)
        yp=rarr(nrg)*(1.d0-1.d-4)*sint
        zp=rarr(nrg)*(1.d0-1.d-4)*cost
        
        ux=0.d0
        uy=-sin(peelthet)
        uz=-cos(peelthet)

        rsq=xp**2+yp**2+zp**2
        rtot=sqrt(rsq)

        ir=nrg-1
        call locate(thetarr,ntg,thet,it)
        ip=1

        tau=1.d0
        tsum=0.d0

        call tauint_peel(iphot,tsum,ir,it,ip)
        
        write(99,*) yp/autors,zp/autors

     end do

     close(99)

  end do

  return

end subroutine photospherepeel
