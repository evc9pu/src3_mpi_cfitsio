subroutine densenv_read(rad,thetin,costin,sina,phi,lor,upr,lot,upt,denll,denuu,denlu,denul,dens,radamb,id,inhole,x0,vz,vcyn,ir,it)

  ! calculates density in envelope. called during grid setup
  ! so doesn't need to be optimized.

  ! history:
  ! 00/09/06 (baw): write subroutine, modified from opacenv.f
  ! 01/02/17 (baw): combine opacitybub subroutine into here.

  use tts_mod
  use grid_mod
  use opacin_mod
  use constants

  implicit none

  real(8) :: rad,dens,cosa,phi,thet,thetin,densw
  real(8) :: costin,sina,xmu,rx,rp,zup,xmu0,factor
  real(8) :: rx2,xmu0new,zlo,zp,radamb,x0,vz,vcyn,vphi
  real(8) :: lor,upr,lot,upt,denll,denuu,denlu,denul,denstu,denstl
  real(8) :: llor,lupr,ldenll,ldenuu,ldenul,ldenlu

  integer :: iflag,id,inhole
  integer :: ir,it

  if(.not.is_envelope(id).and..not.is_cavity(id)) then
     write(*,*) "ERROR: dust index",id,"is not for an envelope or cavity"
     stop
  end if

  inhole=0

  ! for now, assuming +z = -z
  thet=thetin
  cosa=abs(costin)
  if (thet.gt.pihalf) thet=pi-thet
  rp=rad*sina
  zp=rad*cosa

  zup=c1e*rp**ex1+z01
  xmu=zp/rad
  if (rad.gt.rmine) then

     llor=log10(lor)
     lupr=log10(upr)
     ldenll=log10(denll)
     ldenuu=log10(denuu)
     ldenul=log10(denul)
     ldenlu=log10(denlu)

     denstu=ldenul+(ldenuu-ldenul)/(lupr-llor)*(log10(rad)-llor)
     denstl=ldenll+(ldenlu-ldenll)/(lupr-llor)*(log10(rad)-llor)
     dens=denstl+(denstu-denstl)/(upt-lot)*(thetin-lot)
     dens=10.d0**dens

!     print*,dens

     if(istream.eq.1) then
        rx2=rad/rchole
        call zerod(xmu,rx2,xmu0new,iflag)
     end if

  else
     dens=0.d0
  endif

!! constant density in gap--will be overwritten in hole and by ambient density (I want this)
!  if (igape.eq.1.and.igapedens.eq.0) then
!!! igapedens=1 means scaling density in gap; igapedens=0 means density = constant in gap
!  if (rad.lt.rgape2.and.rad.gt.rgape1) then
!    dens=rhogape
!  endif
!  endif

  ! compare to ambient density; if rho lt rhoamb, set rho to rhoamb.
  if(rad.gt.rmine.and.dens.lt.rhoamb) then
     dens=rhoamb
     if (rad.lt.radamb) radamb=rad
  end if

  x0=0.d0
  vz=0.d0
  vcyn=0.d0

  call densdiskwind(rp,zp,densw,x0,vz,vcyn,vphi,ir,it)

  if(ihole.eq.1) then
     if (istream.eq.1.and.xmu0new.gt.windmu0) then
        inhole=1
        dens=densw
     end if
     if (ipoly.eq.1) then
        zlo=c2e*rp**ex2+z02
        ! baw: 2/5/99 new.  doing this so holes can intersect and create
        ! different shapes!
        if (zp.gt.zlo) then
           inhole=1
           dens=rhoconst2*rad**(-exf)
           ! if(zp.lt.zflowmin) dens=0.d0
        end if
        if(zp.gt.zup) then
           inhole=1
           dens=rhoconst1*rad**(-exf)
           ! if(zp.lt.zflowmin) dens=0.d0
        end if
     end if
     if(ibub.eq.1) then
        if (rad.lt.zbub2*xmu**nbub) then
           ! idust2=4
           inhole=1
           if (rad.lt.zbub1*xmu**nbub) then
              dens=rhoconst1
           else
              dens=rhoconst2
           end if
        end if
        ! if(rp.lt.roa.and.rad.gt..8d0*zbub1) dens=rhoconst1
        if(xmu.gt.cosbuboa) dens=rhoconst1
     end if
  end if

!!  scaled density in gap
!  if (igape.eq.1.and.igapedens.eq.1) then
!  if (rad.lt.rgape2.and.rad.gt.rgape1) then
!     dens=dens*fracte 
!  endif
!  endif

  ! make density stop at r>rmax and r<rmine
  ! extend outflow to outside of the core
  if (ifclump.eq.'NO'.or.ifclump.eq.'no') then
     if (rad.ge.rcore*autors.and.inhole.eq.0) dens=0.d0
  else
     if (rad.ge.rcore*autors.and.inhole.eq.0) dens=clumpden*(rad/rcore/autors)**(-clumppow)
  end if
  if (rad.gt.rmax) dens=0.d0
  if (rad.le.rmine) dens=0.d0

  ! set cavity or envelope density based on dust type
  ! if the envelope is higher than the dust sublimation temperature,
  ! for now just consider it as id=10
  if (is_envelope(id)) then
     if (inhole.eq.1) then
        dens=0.d0
     else if (tdave(ir,it,1).ge.1400.d0/11605.d0) then
        dens=0.d0
     end if
  else
     if (inhole.eq.0.and.tdave(ir,it,1).lt.1400.d0/11605.d0) then
        dens=0.d0
     end if
  end if

!!$  if (is_cavity(id).and.inhole.eq.0) dens=0.d0
!!$  if (is_envelope(id).and.inhole.eq.1) dens=0.d0

  return

end subroutine densenv_read

