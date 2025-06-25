
!!revised version of the wind solution 052612

subroutine densdiskwind(rpin,zpin,dens,x0,vz,vcyn,vphi,ir,it)

  use tts_mod
  use opacin_mod
  use constants
  use grid_mod
  
  implicit none

  real(8) :: rpin,zpin,dens,x0,rp,zp,vz,vcyn,vphi
  real(8) :: cyn_c0,co_phi,cyn_max0
  real(8) :: cyn_c,rad_max,cyn_max,x_max,co_eta,co_delta
  real(8) :: x_t,dcyndz,dRmdZm,dRcdZc,rhoch,v_Kc,zd,co_z
  real(8) :: co_zmax,co_rmax,dens1,zd1,thet
  real(8) :: x0_old,co_r,tanwind2,temp1,temp2,co_rc,co_m
  real(8) :: Omega0,co_ra,co_a,co_b,co_c,co_zra,co_rmax0,co_rc0
  real(8) :: co_r1

  integer :: irr,irr_old,ii,irr_l,irr_u
  integer :: ir,it,count

  co_a=7.8d0
  co_b=1.d0
  co_c=0.023d0

  co_zra=(exp(1.d0/co_a)-co_b)/co_c

  cyn_c0=rstar*rsol
  
  rp=rpin*rstar*rsol/cyn_c0
  zp=zpin*rstar*rsol/cyn_c0

  rad_max=zp/windmu0+rchole*(1.d0-windmu0**2)
  cyn_max=sqrt(rad_max**2-zp**2)
  cyn_c=1.d0+14.d0*log(1.d0+0.07d0*zp)
  if (cyn_max.lt.cyn_c) then
     print*,'outflow cavity wall smaller than the inner cavity wall!!!!!!'
     stop
  end if
  if (rp.lt.cyn_c) then
     dens=0.d0
     x0=0.d0
     vz=0.d0
     vcyn=0.d0
     vphi=0.d0
     return
  end if
  if (zp/rp.lt.zdisk(ir)/cyndisk(ir)) then
     dens=0.d0
     x0=1.d0
     vz=0.d0
     vcyn=0.d0
     vphi=0.d0
     return
  end if
  if (ravearr(ir).gt.radmax.and.rp.gt.cyn_max) then
     dens=0.d0
     x0=1.d0
     vz=0.d0
     vcyn=0.d0
     vphi=0.d0
     return
  end if

  v_Kc=sqrt(gn*massc*msol/cyn_c0)
  rhoch=mdotdisk*msol/3600.d0/24.d0/365.25d0*fw*fwesc
  rhoch=rhoch/4.d0/pi/(cyn_c0)**2/v_Kc/log(x_max0)

  tanwind2=(1.d0-windmu0**2)/windmu0**2

!!$  x0=log(rp/cyn_c)/log(cyn_max/cyn_c)*log(x_max0)
!!$  x0=10.d0**x0 ! initial guess of x0
!!$  x0_old=1.
!!$  count=0
!!$  do while (abs(x0-x0_old)/x0.gt.1.d-3)
!!$     x0_old=x0
!!$     co_delta=log(x0)/log(x_max0)
!!$     call locate(cyndisk,nrg-1,x0,irr)
!!$     if (irr.gt.1) then
!!$        if (cyndisk(irr).gt.x0) then
!!$           irr_l=irr-1
!!$           irr_u=irr
!!$        else
!!$           irr_l=irr
!!$           irr_u=irr+1
!!$        end if
!!$        zd=zdisk(irr_l)+(zdisk(irr_u)-zdisk(irr_l))/(cyndisk(irr_u)-cyndisk(irr_l))*(x0-cyndisk(irr_l))
!!$     else
!!$        zd=0.d0
!!$     end if
!!$     co_z=(zp-zd)/x0
!!$     co_rmax=1.d0+co_z**2*tanwind2+2.d0*co_z*rchole/x_max0*(1.d0-windmu0**2)/windmu0+2.d0*co_z*zdmax/x_max0*tanwind2
!!$     co_rmax=sqrt(co_rmax)
!!$     co_rc=1.d0+14.d0*log(1.d0+0.07d0*co_z)
!!$     co_r=co_rc**(1.d0-co_delta)*co_rmax**co_delta
!!$     x0=rp/co_r
!!$     count=count+1
!!$     print*,rp,zp,x0,co_r,zd
!!$     if (count.gt.10000) then
!!$        print*,'x0 cannot converge',x0,x0_old
!!$        stop
!!$     end if
!!$  end do

  do irr=1,nrg-1
     zd=zdisk(irr)
     x0=cyndisk(irr)
     co_delta=log(x0)/log(x_max0)
     co_z=(zp-zd)/x0
     if (co_z.lt.0.d0) co_z=0.d0
     co_rmax=1.d0+co_z**2*tanwind2+2.d0*co_z*rchole/x_max0*(1.d0-windmu0**2)/windmu0+2.d0*co_z*zdmax/x_max0*tanwind2
     co_rmax=sqrt(co_rmax)
     co_rc=1.d0+14.d0*log(1.d0+0.07d0*co_z)
     co_r=co_rc**(1.d0-co_delta)*co_rmax**co_delta
     if (x0*co_r.gt.rp) exit
  end do

  if (it.lt.(ntg-1)/2) then
     thet=thetarr(it+1)
  else
     thet=thetarr(it-1)
  end if
  if (ravearr(ir)*abs(cos(thet)).le.zdisk(ir)) zd=zdisk(ir)

  co_z=(zp-zd)/x0
  if (co_z.lt.0.d0) co_z=0.d0
  co_r=rp/x0
  co_rmax=1.d0+co_z**2*tanwind2+2.d0*co_z*rchole/x_max0*(1.d0-windmu0**2)/windmu0+2.d0*co_z*zdmax/x_max0*tanwind2
  co_rmax=sqrt(co_rmax)
  co_rc=1.d0+14.d0*log(1.d0+0.07d0*co_z)
  x_max=co_rmax*x_max0/co_rc
  x_t=rp/co_rc

  co_rmax0=1.d0+co_zra**2*tanwind2+2.d0*co_zra*rchole/x_max0*(1.d0-windmu0**2)/windmu0+2.d0*co_zra*zdmax/x_max0*tanwind2
  co_rmax0=sqrt(co_rmax0)
  co_rc0=1.d0+14.d0*log(1.d0+0.07d0*co_zra)
  co_ra=co_rc0**(1.d0-co_delta)*co_rmax0**co_delta
  co_r1=co_rc**(1.d0-co_delta)*co_rmax**co_delta
  
  dRmdZm=(co_z*tanwind2+rchole/x_max0*(1.d0-windmu0**2)/windmu0+zdmax/x_max0*tanwind2)/co_rmax
  dRcdZc=14.d0*0.07d0/(1.d0+0.07d0*co_z)
  temp1=dRmdZm*co_z/co_rmax
  temp2=dRcdZc*co_z/co_rc
  co_eta=log(x_max)/log(x_max0)-co_delta*temp1-(1.d0-co_delta)*temp2
  co_eta=1.d0/co_eta
  vz=log(1.01d0+5.d0*co_z**0.8)
  co_m=co_a*log(co_b+co_c*co_z) ! fit the alfven mach number m from BP82 wind. m^(1/2) value in BP82
  if (x_t.lt.1.d0) then
     dens=0.d0
  else
     dens=rhoch*co_eta/vz/rp**2*x0**0.5d0
  end if

!!$  if (ir.eq.333.and.(it.eq.437.or.it.eq.438.or.it.eq.439)) then
!!$     print*,'!!!!',it,rp,x0,zp,zd,vz,co_eta,dens,rhoch/log(1.01d0)/rp**1.5d0
!!$  end if

  dcyndz=(1.d0-co_delta)*(co_rmax/co_rc)**co_delta*dRcdZc &
       & +co_delta*(co_rmax/co_rc)**(co_delta-1.d0)*dRmdZm
  if (x_t.gt.1.d0.and.x_t.lt.x_max.and.dcyndz.lt.0.d0) then
     print*,'densdiskwind setting wrong!!!!!'
     print*,rp,zp,dcyndz,x_t
     print*,co_z,tanwind2,rchole,x_max0,windmu0,zdmax,co_rmax
     stop
  end if

  vz=vz*v_Kc*x0**(-0.5d0)
  vcyn=vz*dcyndz

  Omega0=sqrt(gn*massc*msol/cyn_c0**3)*x0**(-1.5d0)
!  vphi=(1.-(co_r1/co_ra)**2)/(co_m**2-1.d0)+1.d0 ! eq. 12 and 13 in Konigl & Pudritz 2000 chapter
!  vphi=vphi*(co_ra/co_r1)*Omega0*co_ra*x0*cyn_c0
  co_ra=co_rc0! inconsitency here?
  co_r1=co_rc! inconsitency here?
  vphi=(1.-(co_r1/co_ra)**2)/(co_m**2-1.d0)+1.d0 ! eq. 12 and 13 in Konigl & Pudritz 2000 chapter
  vphi=vphi*(co_ra/co_r1)*Omega0*co_ra*x0*cyn_c0

  if (vphi.lt.0.d0) then
     print*,'negative vphi',co_r1/co_ra,co_m
  endif
  
  x0=log(x0)/log(x_max0)

  return

end subroutine densdiskwind
