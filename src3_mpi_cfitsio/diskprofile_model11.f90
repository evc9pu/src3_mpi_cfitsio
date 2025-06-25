subroutine diskprofile(massc_in,rmincgs_in,mdotcgs_in,rmaxd_in,rinfall_in,fd_in,mdotcorecgs_in,sigmadotarr,alpha,ifout,Co_md,Co_Lw,Co_Ld,sigmaarr,zmrarr)

  use constants
  use grid_mod
  use tts_mod
  use dust_mod
  use opacin_mod

  implicit none

  real(8) :: massc_in,rmincgs_in,mdotcgs_in,rmaxd_in,rinfall_in
  real(8) :: fd_in,mdotcorecgs_in,alpha
  real(8) :: Co_mdotd,Co_mdotin,Co_mdotw,Co_Jdotd,Co_Jdotin,Co_Jdotw
  real(8) :: Co_Edotd,Co_md,Co_Ld,Co_Lw,surfden,surfden1,Co_min,Co_Lin
  real(8) :: Co_Omega,Co_Omega1,Co_mdot,Co1,Co2,Co3,Co_Lw1,f
  real(8) :: rad,dr,radcgs,drcgs,sigmadot,par1,par2,par3,Th,Tl,rho
  real(8) :: aKext,aKextl,aKexth,aKext1,zmr,T,rator0,thetemp,muc
  real(8) :: par4,rho1

  real(8) :: sigmadotarr(nrg-1),sigmaarr(nrg-1),zmrarr(nrg-1)

  integer :: ir,itmp,iflag

  logical :: ifout

  rator0=sqrt(30.d0)
  muc=2.3d0

  if (ifout) open(unit=99,file='diskprofile.dat',status='unknown')
  
  Co_mdotd=0.d0
  Co_mdotin=0.d0
  Co_mdotw=0.d0
  Co_Jdotd=0.d0
  Co_Jdotin=0.d0
  Co_Jdotw=0.d0
  Co_Edotd=0.d0
  Co_md=0.d0
  Co_Ld=0.d0
  Co_Lw=0.d0
  surfden1=0.d0
  
  Co_min=msol*massc_in
  Co_Omega=sqrt(gn*Co_min/rmincgs_in**3)
  Co_Omega1=Co_Omega
  Co_mdot=mdotcgs_in+Co_mdotd+Co_mdotw-Co_mdotin
  
  Ldint=0.d0
  hdisk=0.d0
  midrho=0.d0

  sigmaarr=0.d0
  zmrarr=0.d0

  do ir=1,nrg-1

     rad=ravearr(ir)
     dr=rarr(ir+1)-rarr(ir)
     radcgs=rad*rmincgs_in
     drcgs=dr*rmincgs_in
     sigmadot=sigmadotarr(ir)

     if (rad.gt.rddust.and.rad.lt.rmaxd_in) then

        Co1=3.d0*pi-2.d0*pi**2*surfden1*radcgs**2/Co_min
        Co2=Co_mdot-mdotcgs_in*sqrt(msol*massc_in/Co_min)/sqrt(rad)
        Co3=(Co_Jdotd+Co_Jdotw-Co_Jdotin)/Co_Omega/radcgs**2
        f=(Co2-Co3)/Co1
        
        par1=(8.d0*sigt*alpha*k/3.d0/muc/mH/sqrt(2.d0*pi))/f**2 &
             & /(9.d0/4.d0*Co_Omega**3 &
             & +pi**2*gn**2*surfden1**2/Co_Omega/radcgs**2 &
             & -3.d0*pi*gn*surfden1*Co_Omega/radcgs)
        par2=f/alpha*Co_Omega**2/(k/muc/mH)**1.5d0/sqrt(2.d0*pi)
        par3=sqrt(k/muc/mH)/Co_Omega
        par4=4.d0*sigt*muc*mH/3.d0/c/k

!!$        if (ifprint.and.alpha.eq.10.d0) then
!!$           print*,alpha,Co1,Co2,Co3,f,co_mdot,co_mdotd,sigmadot,co_md/msol/massc_in
!!$        end if
        
        itmp=nT-1
        
        Th=Tint(itmp)*11605.d0
        Tl=Tint(itmp-1)*11605.d0
        rho=par2*Th**(-1.5d0)
!        call findrho(rho1,par4*Th**3,rho,iflag)
!        rho=rho1
        call cop(rho,Th,aKext)
        aKexth=aKext/Th**5
!        aKexth=aKext/Th**5/(1.d0+par4/rho*Th**3)
        rho=par2*Tl**(-1.5d0)
!        call findrho(rho1,par4*Tl**3,rho,iflag)
!        rho=rho1
        call cop(rho,Tl,aKext)
        aKextl=aKext/Tl**5
!        aKextl=aKext/Tl**5/(1.d0+par4/rho*Tl**3)
!        print*,'1',Th,Tl,aKexth,aKextl

        do while ((aKexth-par1)*(aKextl-par1).gt.0.d0.and.Tl/11605.d0.gt.Tint(1))
           itmp=itmp-1
           Th=Tint(itmp)*11605.d0
           Tl=Tint(itmp-1)*11605.d0
           rho=par2*Th**(-1.5d0)
!           call findrho(rho1,par4*Th**3,rho,iflag)
!           rho=rho1
           call cop(rho,Th,aKext)
           aKexth=aKext/Th**5
!           aKexth=aKext/Th**5/(1.d0+par4/rho*Th**3)
           rho=par2*Tl**(-1.5d0)
!           call findrho(rho1,par4*Tl**3,rho,iflag)
!           rho=rho1
           call cop(rho,Tl,aKext)
           aKextl=aKext/Tl**5
!           aKextl=aKext/Tl**5/(1.d0+par4/rho*Tl**3)
!           print*,'2',Th,Tl,aKexth,aKextl
        end do

        if (Tl/11605.d0.le.Tint(1)) then
!           print*,'search all table'
!           print*,'R',rad
           Tl=Tint(1)*11605.d0
           rho=par2*Tl**(-1.5d0)
           call cop(rho,Tl,aKext)
           aKextl=aKext/Tl**5
           Th=Tint(nT-1)*11605.d0
           rho=par2*Th**(-1.5d0)
           call cop(rho,Th,aKext)
           aKexth=aKext/Th**5
           if (abs(log10(aKexth)-log10(par1)).lt.abs(log10(aKextl)-log10(par1))) then
              Tl=Th
              aKextl=aKexth
           else
              Th=Tl
              aKexth=aKextl
           end if           
!           print*,'T',Th
        end if
        
        if (aKexth.eq.aKextl) then
           T=Th
        else
           T=Th*10.d0**(log10(Th/Tl)/log10(aKexth/aKextl)*log10(par1/aKexth))
        end if
!        if (T.lt.Tl) print*,T,Th,Tl,aKexth,aKextl,par1
        rho=par2*T**(-1.5d0)
!        call findrho(rho1,par4*T**3,rho,iflag)
!        rho=rho1
        aKext=par1*T**5
!        call cop(rho,T,aKext)
!        print*,aKext,par1*T**5*(1.d0+par4/rho*T**3)
        zmr=par3*T**0.5d0
!        zmr=par3*T**0.5d0*(1.d0+par4/rho*T**3)**0.5d0
        if (rho.lt.0.d0) rho=0.d0
        surfden=rho*zmr*sqrt(2.d0*pi)
        zmr=zmr/rmincgs_in
        if (ifout) then
           hdisk(ir)=zmr
           midrho(ir)=rho
        end if
        sigmaarr(ir)=surfden
        zmrarr(ir)=zmr
!        print*,'3',T,rho,par2,par1,par3,zmr,surfden


        surfden1=surfden
        Co_md=Co_md+2.d0*pi*surfden*radcgs*drcgs
        Co_min=massc_in*msol+Co_md
        if (Co_md.lt.0.d0) Co_min=massc_in*msol
        Co_Omega=sqrt(gn*Co_min/radcgs**3)
        Co_Omega1=sqrt(gn*massc_in*msol/radcgs**3)
        if (surfden.gt.0.d0) then
           Co_Ld=Co_Ld+16.d0*pi*sigt/3.d0*T**4/aKext/surfden*radcgs*drcgs/sqrt(2.d0*pi)
        end if
        if (ifout) Ldint(ir)=Co_Ld

!        if (ifout.and.ifprint) print*,Co_Ld,T,aKext,surfden,radcgs,drcgs
        
        Co_mdotd=Co_mdotd+sigmadot*2.d0*pi*radcgs*drcgs
        Co_Jdotd=Co_Jdotd+sigmadot*3.d0*pi*radcgs**3*drcgs*Co_Omega
        Co_Edotd=Co_Edotd+sigmadot*2.d0*pi*radcgs**3*drcgs*Co_Omega**2
        
        if (rad.le.rinfall_in) then
           Co_mdotw=fw*fwesc*mdotcgs_in*log(rad)/log(rinfall_in)
           Co_Jdotw=Co_Jdotw+fw*fwesc*mdotcgs_in/log(rinfall_in)*radcgs*drcgs*Co_Omega*rator0**2
           Co_Lw=Co_Lw+fw*fwesc*mdotcgs_in/log(rinfall_in)*radcgs*drcgs*Co_Omega**2*(rator0**2-1.5d0)
           Co_mdotin=Co_mdotin+0.d0
           Co_Jdotin=Co_Jdotin+0.d0
        else
           Co_mdotw=Co_mdotw+0.d0
           Co_Jdotw=Co_Jdotw+0.d0
           Co_Lw=Co_Lw+0.d0
           Co_mdotin=mdotcorecgs_in*(cos(thet1*deg2rad)-sqrt(1.d0-rad/rmaxd_in))
           Co_Jdotin=mdotcorecgs_in/3.d0 &
                & *sqrt(gn*msol*massc_in*(1.d0+fd_in)*rmaxd_in*rmincgs_in) &
                & *(((sin(thet1*deg2rad))**2+2.d0)*cos(thet1*deg2rad)- &
                & (rad/rmaxd_in+2.d0)*sqrt(1.d0-rad/rmaxd_in))
           thetemp=asin(sqrt(rad/rmaxd_in))
           if (ifout) then
              Ldint(ir)=Ldint(ir)+ &
                   & gn*mdotcorecgs_in*msol*massc_in*(1.d0+fd_in)/rmaxd_in/rmincgs_in* &
                   & (log(tan(thetemp/2.d0))-log(tan(thet1/2.d0*deg2rad)) &
                   & +cos(thetemp)-cos(thet1*deg2rad))
           end if
        end if

        Co_mdot=mdotcgs_in+Co_mdotd+Co_mdotw-Co_mdotin

        if (ifout) write(99,'(8(E15.8))') rad/autors,zmr/autors,T,rho,surfden,sigmadot,Co_mdot,aKext
        
     else if (rad.ge.rmaxd_in) then
        if (ifout) Ldint(ir)=Ldint(ir-1)
     end if

  end do

  if (ifout) close(99)

  Co_Lin=-(gn*mdotcorecgs_in*msol*massc_in*(1.d0+fd_in)/rmaxd_in/rmincgs_in) &
       & *(log(tan(thet1/2.d0*deg2rad))+cos(thet1*deg2rad) )

  Co_lw1=gn*massc_in*msol*mdotcgs_in/2.d0/rmincgs_in-Co_Lin-Co_Ld+Co_Edotd

  Co_Ld=Co_Ld+Co_Lin
  if (ifout) then
     do ir=1,nrg-1
        Ldint(ir)=Ldint(ir)/Ldint(nrg-1)
     end do
     open(unit=98,file='Ldint.dat',status='unknown')
     do ir=1,nrg-1
        rad=ravearr(ir)
        write(98,*) rad,rad/autors,Ldint(ir)
     end do
     close(98)
  end if

  Co_md=Co_md/msol/massc_in
  
  return

end subroutine diskprofile

!!$subroutine findrho(rho,co_a,co_b,iflag)
!!$
!!$  implicit none
!!$  
!!$  real(8) :: rho,co_a,co_b,rho1
!!$
!!$  integer :: iiter,iflag
!!$
!!$  rho=co_b
!!$
!!$  iiter=0
!!$
!!$10 rho1=rho
!!$
!!$  rho=co_b/(1.d0+co_a/rho)**1.5d0
!!$
!!$  print*,rho
!!$
!!$  iiter=iiter+1
!!$  
!!$  if (iiter.gt.40) goto 20
!!$
!!$  if (abs(rho1-rho).gt.rho*1.d-2) goto 10
!!$
!!$20 if (iiter.gt.40) print*,'iteration of rho',iiter,rho1,rho,co_a,co_b
!!$
!!$  return
!!$
!!$end subroutine findrho

subroutine findrho(xmu0,co_a,co_b,iflag)

  use constants

  implicit none

  real(8) :: a1,a2,a3,p,q,r,theta,xmu0
  real(8) :: xmu01,xmu02,xmu03,factor,t1,t2

  real(8) :: co_a,co_b

  integer :: iflag

  a1 = 3.d0*co_a
  a2 = 3.d0*co_a**2-co_b**2
  a3 = co_a**3

  p=a2-a1**2/3.d0
  q=2.d0/27.d0*a1**3-a1*a2/3.d0+a3

  factor=(q/2.d0)**2+(p/3.d0)**3

  if (factor.le.0.d0) then

     if (p.gt.0.d0) then
        print*,'p wrong'
        stop
     end if

     r=sqrt(-(p/3.d0)**3)
     theta=acos(-q/2.d0/r)/3.d0

     xmu01=2.d0*r**(1.d0/3.d0)*cos(theta)
     xmu02=2.d0*r**(1.d0/3.d0)*cos(theta+2.d0*pi/3.d0)
     xmu03=2.d0*r**(1.d0/3.d0)*cos(theta+4.d0*pi/3.d0)

     if (-q.le.0.d0) then

        if (xmu01.gt.0.d0.and.xmu02.gt.0.d0) then
           xmu0=max(xmu01,xmu02)
        else if (xmu01.gt.0.d0.and.xmu03.gt.0.d0) then
           xmu0=max(xmu01,xmu03)
        else if (xmu02.gt.0.d0.and.xmu03.gt.0.d0) then
           xmu0=max(xmu02,xmu03)
        end if
        iflag=1

     else
        if (xmu01.gt.0.d0) xmu0=xmu01
        if (xmu02.gt.0.d0) xmu0=xmu02
        if (xmu03.gt.0.d0) xmu0=xmu03
        iflag=2
     end if

  else

     t1=-q/2.d0+sqrt(factor)
     t2=-q/2.d0-sqrt(factor)
     xmu01=dsign((abs(t1))**(1.d0/3.d0),t1)+dsign((abs(t2))**(1.d0/3.d0),t2)
     xmu0=xmu01
     iflag=3

  end if

  xmu0=xmu0-a1/3.d0

  if (xmu0.le.0.d0) then
     xmu0=co_a**3/co_b**2
  end if

!  print*,iflag,xmu0*(1.d0+co_a/xmu0)**1.5d0-co_b

  return

end subroutine findrho
