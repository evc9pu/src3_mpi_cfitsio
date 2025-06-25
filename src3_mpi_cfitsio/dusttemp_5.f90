! ********************************************************************
! Function to calculate dust temperature from condition of radiative
! equilibrium
! ********************************************************************

! with expansion cooling

real(8) function dusttemp(ir,it,ip,na,T0,rstar,tstar,nphot,idust2,mcell,gradtr,gradtt)

  use grid_mod
  use dust_mod
  use constants
  implicit none

  integer :: ir,it,ip             ! grid cell index
  integer :: nphot,idust2,numtot,idust3,i
  integer :: mini,maxi
  real(8) :: T0,Td,Told
  real(8) :: rstar,tstar,na
  real(8) :: minkap,maxkap
  real(8) :: KapP
  real(8) :: kappaT(nT)
  real(8) :: kappaTl,kappaTu,Tl,Tu
  real(8) :: mcell
  real(8) :: divv,vr,vthet,gradtr,gradtt,gradrhor,gradrhot
  real(8) :: adi,adv,emi,abs,absl,absu

  real(8),parameter :: eps = 1.d-5

  real(8) :: Tmax

  logical :: ifout

  ifout=.false.

  numtot=0

  Tmax=199999.d0

  idust3=findopac(ir,it,ip,idust2)

  divv=divvarr(ir,it,1)

  vr=vrarr(ir,it,1)
  vthet=vthetarr(ir,it,1)

  do i=1,nT

     emi=Tint(i)**4*kappap(i,idust3)
     adi=(k/4.d0/sigt/2.3d0/mH)/kappav(idust3)/11605.d0**3* &
          & Tint(i)*divv/rstar/rsol
     adv=(k*3.d0/8.d0/sigt/2.3d0/mH)/kappav(idust3)/11605.d0**3* &
          & (vr*gradtr+vthet*gradtt)/rstar/rsol
     kappaT(i)=emi+adi+adv

  end do

  if (mcell.gt.0.d0) then

     abs=Tstar**4*(na*pi*Rstar*Rstar)/(nphot*kappav(idust3))
     abs=abs*rsol**2/mcell**4

     if (na.eq.0.d0) then
        Td=0.1d0/11605.d0
        goto 30
     end if

     minkap=kappaT(1)
     maxkap=kappaT(1)
     mini=1
     maxi=1
     do i=1,nT
        if (minkap.gt.kappaT(i)) then
           minkap=kappaT(i)
           mini=i
        end if
        if (maxkap.lt.kappaT(i)) then
           maxkap=kappaT(i)
           maxi=i
        end if
     end do

     if (maxkap.lt.abs) then
        Td=Tint(maxi)
!        print*,'abs is too high',idust3,maxkap,abs,maxi
        Td=min(Td,Tmax/11605.d0)
        Td=max(Td,0.1d0/11605.d0)
        absu=maxkap
        absl=maxkap
     end if

     if (minkap.gt.abs) then
        Td=Tint(mini)
!        print*,'abs is too low', idust3,minkap,abs,mini
        Td=min(Td,Tmax/11605.d0)
        Td=max(Td,0.1d0/11605.d0)
        absu=minkap
        absl=minkap
     end if

     if (mini.lt.maxi) then
        do i=mini,maxi-1
           if (kappaT(i).lt.abs.and.kappaT(i+1).ge.abs) then
              kappaTl=kappaT(i)
              kappaTu=kappaT(i+1)
              Tl=Tint(i)
              Tu=Tint(i+1)
              Td=Tl+(Tu-Tl)/(kappaTu-kappaTl)*(abs-kappaTl)
              Td=min(Td,Tmax/11605.d0)
              Td=max(Td,0.1d0/11605.d0)
              absu=kappaT(i+1)
              absl=kappaT(i)
              goto 30
           end if
        end do
     else
        do i=maxi,mini-1
           if (kappaT(i).gt.abs.and.kappaT(i+1).le.abs) then
              kappaTl=kappaT(i)
              kappaTu=kappaT(i+1)
              Tl=Tint(i)
              Tu=Tint(i+1)
              Td=Tl+(Tu-Tl)/(kappaTu-kappaTl)*(abs-kappaTl)
              Td=min(Td,Tmax/11605.d0)
              Td=max(Td,0.1d0/11605.d0)
              absu=kappaT(i+1)
              absl=kappaT(i)
              goto 30
           end if
        end do
     end if

  else

     Td=0.1d0/11605.d0

  end if

30  dusttemp=Td

  emi=Td**4*kapP(Td,idust3)
  adi=(k/4.d0/sigt/2.3d0/mH)/kappav(idust3)/11605.d0**3* &
       & Td*divv/rstar/rsol
  adv=(k*3.d0/8.d0/sigt/2.3d0/mH)/kappav(idust3)/11605.d0**3* &
       & (vr*gradtr+vthet*gradtt)/rstar/rsol
  
  adiarr(ir,it,ip)=adi/emi
  advarr(ir,it,ip)=adv/emi

  return

end function dusttemp

