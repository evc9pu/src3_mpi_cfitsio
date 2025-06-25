subroutine gridset(iter,i1)
  
  use output_mod
  use log_mod
  
  ! set up grid:  3-D spherical grid with variable spacing in r and
  ! theta (set by exponents, rexp, texp)
  ! make a bunch of arrays necessary for find_wall subroutine
  
  ! history:
  ! 00/09/06 (baw): write subroutine
  
  use tts_mod
  use grid_mod
  use opacin_mod
  use out_mod
  use dust_mod
  use spot_mod
  use constants
  use random
  use ttsre_mpi_mod

  implicit none

  real(8) :: dr,dt,dp,rad,phi,thet,densd,dense,vol
  real(8) :: cost,sint,eps,tau,rming1,tiny
  real(8) :: rmincgs,taud,taue,mu0,tauv,taufmin,mu,thetatm
  real(8) :: thetmin,tauzrat,sigmu,rsonr,Tdisk,taudisk
  real(8) :: mudisk,a1,tauw,rtop,rfac,tauRr,taumid,tauRth
  real(8) :: tauRos,res,thettop,maxdens,tauRmu,radamb,rbeg2
  real(8) :: tauave,lstar,mdot,lacc,lacc2,tave,Vc
  real(8) :: rminsub,fluxratio
  real(8) :: fplay,rexpnew,ltot2,diskscale,rmintest
  real(8) :: dr3,dcost,r,masssg,x,z,masstot,mu0tmp
  real(8) :: aave,tmp,taup,dphi,totvol,kappafave
  real(8) :: fl,const,fm1,fm2,fm3,fm4,densav,densav2,rntot,std
  real(8) :: maxdensvar,mindensvar,taug
  real(8) :: denuu,denll,denul,denlu,lot,upr,lor,upt,vu,vl
  real(8) :: thetatm1,thetatm2,thetmin1,thetmin2
  real(8) :: tdust_temp,x0,x0_fp,t_fp,vz,vcyn,vphi,dens1
  real(8) :: rad1,rad2,theta1,theta2,vr1,vr2,vthet1,vthet2
  real(8) :: vrarru(nrg,ntg,npg),vrarrw(nrg,ntg,npg),vrarre(nrg,ntg,npg)
  real(8) :: vthetarru(nrg,ntg,npg),vthetarrw(nrg,ntg,npg)
  real(8) :: vthetarre(nrg,ntg,npg),vphiarrw(nrg,ntg,npg)
  real(8) :: divvarru(nrg,ntg,npg),divvarrw(nrg,ntg,npg),divvarre(nrg,ntg,npg)
  real(8) :: vfall,cost0,vew(nrg)
  real(8) :: thetatm3,thetatm4,thetatm5,gspl
  real(8) :: cpusec
  real(8) :: rcav,thetcav
  real(8) :: r_foot

  real(8),pointer :: radtab(:),thetatab(:),dentab(:,:),vtab(:)

  ! real :: xyarr(200,200),x,y,xmax,dx,r,dens1,z
  real(8) :: chiR,ierfc,erfc
  integer :: ir,it,ip,nttmp,nr0,ntatm,nrwall,irbeg
  integer :: dcount,id,i,iter,igas,inhole,j,i1,icent,nf
  integer :: indr,indrl,indru,indt,indtl,indtu,indv
  integer :: ntatm1,ntatm2,ntatm3,ntatm4,ntatm5
  integer :: idmax,idmin
  integer :: ir_fp,it_fp
  integer :: iflag
  integer :: thetcase
  integer :: idust2

  integer :: peelid

  character(len=4) :: suffix

  call diskdat_section('gridset')

!  call findopacid()
  
  tiny=1.d-15
  rmincgs=rstar*rsol
  radamb=rmax

  do id=1,nopac
     kapd(id)=kappav(id)*rmincgs
  end do
  ! set up the grid only in the first iteration
  if (iter.lt.iterstart+1) then

  ! I have some gridding below in case there's a gas disk, but that's a future update
  ! for now, set gas disk to no.
  igas=0

  if (ifprint) print*,'idiskacc',idiskacc

  if (idiskacc.eq.1) then

     ! set disk accretion fraction based on alpha disk paramters
     ! see eqns 5 and 6 in Whitney et al. 2003, apj, 591, 1049
     lstar=4.d0*pi*sigt*tstar**4*rmincgs**2
     ! 2011 july 6. modification, only counting non-hotspot photons
   !  if (fspot.lt.1.d0) then
     !  lstar=4.d0*pi*sigt*tstar**4*rmincgs**2*(1-fspot)
   !  else
  
     if (ialpha.eq.1) then
        ! if specifying alphad in input instead of mdotdisk
        if (ifprint) print*,'calculating disk accretion rate from alphad)'
        if (ifprint) print*,'ialpha,alphad',ialpha,alphad
        Vc=dsqrt(gn*msol*massc/rmincgs)
        ! I think this will work...20090731
        mdot=0.d0
        do id=1,ndg
           if(is_disk(id)) then
              mdot=mdot+dsqrt(18.d0*pi**3.d0)*alphad*Vc*rho0(id)*(z1(id)*rmincgs)**3.d0/rmincgs
           end if
        end do
     else
        ! if specifying mdotdisk in input
        mdot=mdotdisk/3600.d0/24.d0/365.25d0*msol
     end if

     ! accretion in disk
     lacc=gn*(msol*massc)*mdot/2.d0/(rddust*rmincgs)

     ! accretion in shock, the optically thick part from the heated
     ! atmosphere  (calvet & gullbring 1998).
     ! we will split it into half blackbody at Tshock, and half x-rays.
     ! assumes co-Rotation radius is 5 Rstar.
     if (massd.gt.0.d0) then
        lacc2=gn*(msol*massc)*mdot*(1.d0/rmincgs-1.d0/(rmincgs*rtrunc))
     else
        lacc2=0.d0
     end if
          
     ! total luminosity
     ltot=lacc+lacc2+lstar+l_isrf
     if (ifprint) print*,'lacc,lacc2,lstar,l_isrf',lacc/lsun,lacc2/lsun,lstar/lsun,l_isrf/lsun
!     print*,'ltot',ltot/lsun

     accfrac=lacc/ltot
     accfrac2=lacc2/ltot
     isrf_frac=l_isrf/ltot
     mdot=mdot*3600.d0*24.d0*365.25d0/msol
     if (ifprint) print*, 'including disk accretion luminosity '
     if (ifprint) print*, 'disk accretion rate (msol/yr) ',mdot
     if (ifprint) print*, 'fraction of acc luminosity from disk ',accfrac
     if (ifprint) print*, 'fraction of lum. on stellar hotspot ',accfrac2
     if (ifprint) print*, 'total system luminosity (Lsun)',ltot/lsun

     call diskdat_write('isrf_frac',isrf_frac, 'fractional ISRF lum')
     call diskdat_write('diskacc',.true.,'whether disk luminosity is included')
     call diskdat_write('mdot',mdot,'disk accretion rate (msol/yr)')
     call diskdat_write('accfrac',accfrac,'fraction of lum. from disk')
     call diskdat_write('accfrac2',accfrac2,'fraction of lum. on stellar hotspot')
     call diskdat_write('ltot',ltot/lsun,'total system lum. (Lsun)')


     ! f=calvet's filling factor, ranges from 0.01 to 0.001
     ! median is 0.007
!     if (fspot.gt.0.d0.and.lacc2.gt.0.d0) then
     if (fspot.gt.0.d0.and.lacc2.ge.0.d0) then
       fluxratio=0.5d0*lacc2/lstar/fspot
       tshock=tstar*(1.d0+fluxratio)**0.25d0
       tspot=tshock
       if (nspot.eq.1) then
          thspot=acos(1.d0-2.d0*fspot)*rad2deg
       else if (nspot.eq.2) then
          thspot=acos(1.d0-fspot)*rad2deg
       else
          print*,'error!  nspot not set,stopping program'
          stop
       end if
       
       if (ifprint) print*,'fspot',fspot
       if (ifprint) print*,'number and size of spots ',nspot,thspot
       if (ifprint) print*,'calculated thermal shock T ',tshock
       !       stop
       
       call diskdat_write('tshock',tshock,'calculated thermal shock T')
       
       fplay=0.5d0*lacc2/fspot/4.d0/pi/rmincgs**2
       if (ifprint) print*,'log flux shock ',log10(fplay)
    else
       thspot = 0.d0
    end if
     
    if(partial_peeloff) call output_accretion(lacc, lacc2)
    
 else
    
     lstar=4.d0*pi*sigt*tstar**4*rmincgs**2
     ltot=lstar+l_isrf
     isrf_frac=l_isrf/ltot
     mdot=0.d0
     accfrac=0.d0
     lacc=0.d0
     lacc2=0.d0

     call diskdat_write('diskacc',.false., 'whether disk luminosity is included')
     call diskdat_write('ltot',lstar/lsun, 'total stellar lum. (Lsun)')
     call diskdat_write('ltot',ltot/lsun, 'total system lum. (Lsun)')
     call diskdat_write('isrf_frac',isrf_frac,'fractional ISRF lum')

     if (ifprint) print*, 'disk accretion rate (msol/yr) ',mdot
     if (ifprint) print*, 'fraction of luminosity from disk ',accfrac
     if (ifprint) print*,'lstar/lsun,l_isrf/lsun',lstar/lsun,l_isrf/lsun
     if (ifprint) print*, 'total stellar luminosity (Lsun)',lstar/lsun
     if (ifprint) print*, 'total system luminosity (Lsun)',ltot/lsun
     if (ifprint) print*, 'isrf_frac',l_isrf/ltot

  end if
  
  !BAW TESTING!!!!   2010 apr 28.   delete all these lines after testing
!  lacc2=0.5*lstar
!  accfrac2=lacc2/ltot
!  tspot=tstar*2.5
!  tshock=tspot
!  ! total luminosity
!  ltot=ltot+lacc2
!  print*,'lacc,lacc2,lstar,l_isrf',lacc,lacc2,lstar,l_isrf
!  print*,'ltot',ltot
!  if (nspot.eq.1) then
!     thspot=acos(1.d0-2.d0*fspot)*rad2deg
!  else if (nspot.eq.2) then
!     thspot=acos(1.d0-fspot)*rad2deg
!  else
!     print*,'error!  nspot not set,stopping program'
!     stop
!  end if
!  print*,'number and size of spots ',nspot,thspot
!  print*,'calculated thermal shock T ',tshock
!!  stop 

  ! calculate an average stellar temp based on luminosity
  ! of star+accretion
  tave=tstar*(1+0.5*lacc2/lstar)**0.25
  if (ifprint) print*,'average stellar temp is now',tave

  if (ifprint) print*,'tshock,tave',tshock,tave
  if (ifprint) print*,'note: should be same if fspot=1'
 ! if (fspot.eq.1) tshock=tave

  if (fspot.eq.1.0.and.iplanckst.eq.0) then
    print*,'WARNING: due to accretion, emitting from tstar = ',tave
    print*, 'you may want to change stellar atmosphere file to agree with this!'
  endif

  rminsub=(tsub/tave)**(-2.1d0)
  call diskdat_write('Tave',Tave,'averave Tstar including hotspot')

  ! reset tstar to tave, since it is used for temperature
  ! calculation (sets the luminosity scale).
  tstarave=tave
  ! tstar=tave

  if (rminsub.gt.rmaxd.and.massd.gt.0.0d0) then
     print*,'WARNING!!!, dust sublimation radius larger than RMAXD'
!!$     print*,''
!!$     print*,'****************'
!!$     print*,'ERROR.  dust sublimation radius larger than RMAXD'
!!$     print*,'probably because you chose a hot star and'
!!$     print*,'dust sublimation radius is large.'
!!$     print*,'Increase RMAXD.'
!!$     print*,'RMIND, RMAXD in AU',rminsub/autors,rmaxd/autors
!!$     print*,'Stopping program'
!!$     stop
  end if

  if (irminsub.eq.1) then
     rmine=rminsub*rmine_in
     rmind=rminsub*rmind_in
     rddust=rminsub*rmind_in
     rmind2=rmind2*rminsub
     rmin_sg=rminsub*rmin_sg
     if (ifprint) print*,'resetting Rsub to ',rminsub
     if (ifprint) print*,'inner envelope radius ',rmine
     if (ifprint) print*,'inner disk radius ',rmind
     call diskdat_write('rdust',rmine,'updated dust destruction radius')
     call diskdat_write('rmine',rmine,'updated inner envelope radius')
     call diskdat_write('rmind',rmind,'updated inner disk radius')
  end if

  if (rddust.gt.rmaxd) then
     print*,'ERROR, rddust > rmaxd'
     print*,'rddust,rmaxd (AU)',rddust/autors,rmaxd/autors
     print*,'stopping program'
     stop
  end if

  ! calculate scale height based on Tsub at Rsub
  ! can't do this because we already set disk parameters, would
  ! have to iterate.  save this for hseq code.
  ! testing shows it only increases scale height slightly.
  ! k=1.38e-16
  ! muc=2.3
  ! mH=1.67e-24
  ! rcgs=rminsub*rstar*rsol
  ! mcgs=massc*msol
  ! honr=sqrt(k*Tsub/Gn/mcgs*rcgs/muc/mH)
  ! h=honr*rminsub
  ! print*,'h on r at Rsub',honr
  ! print*,'h (rsub) ',h
  ! zmin=honr*rminsub/rminsub**b
  ! print*,'h_0 based on h(r_sub)',zmin
  ! print*,'assumes Tsub=',Tsub,'    beta =',b

  ! convert opacity to units of cm^2/gm*rstar(cgs) because distance
  ! is in units if 1/rstar, and dtau=kapd*rho*ds

  ! 20070119 BAW, only recalcalate rminsub after first iteration
  ! note we aren't doing this since we aren't calling this subroutine after the first iteration
!  if (iter.gt.1) then
!     it=(ntg+1)/2
!     ir=3
!     do while (tdust2(ir,it,1,1).ge.1600.d0/11605.d0)
!        ir=ir+1
!        print*,tdust2(ir,it,1,1)*11605.d0
!     end do
!     rmintest=rarr(ir-1)

     ! ip=1
     ! ir=1
     ! rmintest=0.d0
     ! do it=1,ntg
     ! do while (densarr(ir,it,ip).eq.0.d0)
     ! ir=ir+1
     ! if(rarr(ir).gt.rmintest) rmintest=rarr(ir)
     ! end do
     ! end do

     ! new rminsub
!     rminsub=rmintest
!     if (rminsub.gt.rmaxd.and.massd.gt.0.0d0) then
!        print*,''
!        print*,'****************'
!        print*,'ERROR.  dust sublimation radius larger than RMAXD'
!        print*,'probably because you chose a hot star and'
!        print*,'dust sublimation radius is large.'
!        print*,'Increase RMAXD.'
!        print*,'RMIND, RMAXD in AU',rminsub/autors,rmaxd/autors
!        print*,'Stopping program'
!        stop
!     end if

!     if (irminsub.eq.1) then
!        rmine=rminsub*rmine_in
!        rmind=rminsub*rmind_in
!        rddust=rminsub*rmind_in
!        print*,'resetting Rsub to ',rminsub
!        print*,'inner envelope radius ',rmine
!        print*,'inner disk radius ',rmind
!        write(12,*) rmine, ' updated dust destruction radius'
!        write(12,*) rmine, ' updated inner envelope radius'
!        write(12,*) rmind, ' updated inner disk radius'
!     end if
!  end if

  if (rddust.gt.rmaxd) then
     print*,'ERROR, rddust > rmaxd'
     print*,'rddust,rmaxd (AU)',rddust/autors,rmaxd/autors
     print*,'stopping program'
     stop
  end if

  ! calculate some constants for the TSC envelope density
  call envset()

  ! make grid include minimum and maximum values.

  if (ifprint) print*, ' '
  if (ifprint) print*, 'grid setup'

  if (ifprint) print*,'rddust,rmine,rmind,rmin',rddust,rmine,rmind,rmin
  ! rgrid
  ! 20050628, change rexp for disks with large inner holes
  if (rddust.gt.0.2d0*rmaxd) then
     rexpnew=2.d0
     if (ifprint) print*,'disk has large inner hole'
     if (ifprint) print*,'changing radial density exponent to rexp=2'
  else
     rexpnew=rexp
  end if
  ! gflag=0   !commented out 20090731
  rarr(1)=rmin

  if (massd.eq.0.d0) then
     ! if there's empty space between rmine and rddust....
     rarr(2)=rddust
     dr=(rmax-rmine)/(dble(nrg-2))**rexpnew
     if (ifprint) print*,'rexp,dr,rmin',rexpnew,dr/autors,rmin/autors
     do ir=3,nrg
        rarr(ir)=rmine+dr*(dble(ir-2))**rexpnew
        ! print*,'rarr ',ir,rarr(ir)/au
     end do

  else

     irbeg=2
     rarr(irbeg)=rddust

     ! lots of fine spacing in inner region of disk

     if (igas.eq.1.and.rtrunc.lt.rddust) then
        ir=2
        if (ifprint) print*,'rtrunc,rddust',rtrunc,rddust
        if (rtrunc.lt.rddust) then
           rarr(ir)=rtrunc
           dr=(rmax-rtrunc)/(dble(nrg/3-1))**rexpnew
           ! print*,'dr',dr
           do while (rarr(ir).lt.rddust)
              ir=ir+1
              rarr(ir)=rtrunc+dr*(dble(ir))**rexpnew
              ! print*,'rarr',rarr(ir)
           end do
           ! reset irbeg and rarr(ir)
           rarr(ir)=rddust
           irbeg=ir
        end if
     end if

     ! a1=1.d0-a

     ! 20090731, do averages over different dust types
     mu0=0.d0
     tauv=0.d0
     tauzrat=0.d0
     tauRos=0.d0
     do id=1,ndg
        if (is_disk(id)) then
           mu0tmp=fmass(id)*z1(id)*Rddust**b(id)/Rddust
           mu0=mu0+mu0tmp
           tauv=tauv+taur(id)
           a1=1.0d0-a(id)
           tauzrat=tauzrat+sqrt(pihalf)*a1*mu0tmp/ &
                & ((rmaxd/rddust)**a1-1.d0)
           tauRos=tauRos+taur(id)*chiR(1400.d0/11605.d0,id)
        end if
     end do
     tauRth=tauzrat*tauRos
     ! taudisk=10.0
     taudisk=min(10.d0,tauRth)
     if (ifprint) print*,'taudisk',taudisk

     if (ifprint) print*,'tauRos,tauv,tauRth',tauRos,tauv,tauRth
     tauw=0.5d0
     rming1=rddust
     aave=0.5d0*(a(1)+a(2))   !average rasdial density exponent for disk thermal grains FIX THISxs
     a1=1.0d0-aave
     dr=(rmax-rming1)/(dble(nrg-irbeg))**rexpnew
     rbeg2=rddust+dr
     if (ifprint) print*,'rbeg2',rbeg2
     if (tauRos.gt.0.51d0) then
        ! okay, just use a average for dust grains
        rfac=(1.d0-(tauw*(1.d0-(rmaxd/rddust)**a1))/tauRos)**(1.d0/a1)
        if (rfac.lt.1.0000001d0) then
           if (ifprint) print*,''
           if (ifprint) print*,'WARNING: not enough precision to resolve inner disk'
           if (ifprint) print*,'rfac = ',rfac
           ! rfac=rfac*2.d0
           rfac=1.000001d0
           if (ifprint) print*,'resetting rfac = ',rfac
           if (ifprint) print*,''
        end if

        mu=sqrt(2.0d0)*ierfc(dble(taudisk/tauRth))
        rtop=Rddust*(1.d0-(taudisk*(1.d0-(rmaxd/rddust)**a1))/ &
             & (tauRos*exp(-0.5d0*mu*mu))      )**(1.d0/a1)
        ! +        (tauRos)      )**(1./a1)

        if (ifprint) print*,'Rddust,taudisk,rmaxd/rddust,tauRos,a1', &
             & Rddust,taudisk,rmaxd/rddust,tauRos,a1
        if (ifprint) print*,'rfac,rtop,mu',rfac,rtop,mu
        ir=irbeg
        rarr(ir)=rddust
        if (rfac*rddust.lt.rbeg2) then
           if (ifprint) print*,'making gaussian spacing in r, rtop',rtop
           do while ((rarr(ir).lt.rtop).and.(ir.lt.(nrg-1)))
              rarr(ir+1)=rfac*rarr(ir)
              rfac=rfac**1.1d0
              ir=ir+1
              ! print*,'rarr',rarr(ir),rtop,ir,rfac
           end do
           if(ir.ge.nrg) then
              write(*,*) 'too many points in rgrid'
              stop
           end if
        end if
     else
        ir=irbeg
     end if

     nrwall=ir
     if (ifprint) print*,'nrwall,rarr(nrwall)',nrwall,rarr(nrwall)

     !     rming1=rarr(nrwall)
     !     dr=(rmax-rming1)/(dble(nrg-nrwall))**rexpnew
     !     print*,'rexp,dr,rmin',rexpnew,dr/autors,rming1
     !     do ir=nrwall+1,nrg
     !        rarr(ir)=rming1+dr*(dble(ir-nrwall))**rexpnew
     !     end do


     rming1=rarr(nrwall)
     dr=(log10(rmax)-log10(rming1))/(dble(nrg-nrwall))
     do ir=nrwall+1,nrg
        rarr(ir)=10.**(log10(rming1)+dr*(dble(ir-nrwall)))
     end do

  end if

  ! set rmine to a grid location
  ! if (gflag.eq.1) then
  ! changed condition, 20090731
  if (rddust.lt.rmine) then
     call locate(rarr,nrg,rmine,ir)
     rarr(ir+1)=rmine
  end if

  ! r-squared array
  do ir=1,nrg
     r2arr(ir)=rarr(ir)**2.d0
     ! print first few points of rarr
     if (ir.lt.10) then
        if (ifprint) print*,'rarr/rmin,ir ',rarr(ir)/rmin,ir
     end if
  end do

  ! r-ave array
  do ir=1,nrg-1
     ravearr(ir)=0.5d0*(rarr(ir)+rarr(ir+1))
  end do

  do ir=1,nrg-1
     if (rarr(ir+1)-rarr(ir).le.0.d0) then
        print*,'rarr wrong!!!',ir,rarr(ir+1),rarr(ir)
        stop
     end if
  end do

  if (ifprint) then
     open(unit=15,file='rarr.dat',status='unknown')
     write(15,*) nrg,' = number of grid points in r'
     write(15,*) &
          & '       index      r/rstar         r(rsub)         r(au)'
     do ir=1,nrg
        write(15,900) ir,rarr(ir), &
             & rarr(ir)/rminsub, &
             & rarr(ir)/autors
     end do
     close(15)
  end if
900 format(i10,3(1x,f15.5))

  if (ifprint) print*,'rarr(nrg),r2arr(nrg)',rarr(nrg),r2arr(nrg)

  ! if (ntg.gt.2) then
  ! set up a tmp array going from 0 - 90 from equ. to pole
  ! make theta=0 bin 5 degrees wide (otherwise, too much noise,
  ! nothing happens there anyway)
!  thettop=5.d0*deg2rad
  thettop=5.d0*deg2rad

  ! taufmin=min(0.1,0.1*kappaf*tauv)
  ! use idust2=2, for kappaf, for disk dust properties
  kappafave=sum(fmass*kappaf(1:8),mask=is_disk(1:8))
  taufmin=min(.001d0,0.001d0*kappafave*tauv)
  
  !     print*,'kappaf',kappaf
  if (ifprint) print*,'kappafave',kappafave
  if (ifprint) print*,'kappaf*tauv',kappafave*tauv
  
  mu=min(mu0*sqrt(2.d0*log(tauv*kappafave/taufmin)),0.5d0)
  if (ifprint) print*,'tauv,kappaf,taufmin',tauv,kappafave,taufmin
  
  ! ************
  ! for disk-only model, let mu=0.5 to sample entire disk height
  ! at high-res
  ! mu=0.5
  ! ***********
  
  thetatm=pihalf-acos(mu)
  ! test
  !     print*,'thetatm',thetatm*rad2deg
  nttmp=(ntg+1)/2
  
  ! *********************
  ! change this for disk or envelope runs
  ! for envelope with 1 degree resolution in polar region
  ! except for first bin (theta=5).
  !     res=1.d0
  res=.25d0
  ! for disk with nothing in the polar region
  ! res = size of angle bin in poles, approximately
  ! res=10.
  ! *********************

  ntatm=nttmp-int((pihalf-thetatm-thettop)*90.d0/pihalf/res)
  thetmin=min(0.1d0*asin(mu0),thetatm/dble(ntatm))
  thetmin1=thetmin/50.d0
  thetmin2=1.d-5

  if (ifprint) print*,'thetmin,thetmin1,thetmin2',thetmin/deg2rad,thetmin1/deg2rad,thetmin2/deg2rad

  !only consider a case with cavity

  if (ifclump.eq.'YES'.or.ifclump.eq.'yes') then
     rcav=rmax
  else
     rcav=rcore*autors
  end if

  thetcav=cos(thet1*deg2rad)*(1.d0-sin(thet1*deg2rad)*sin(thet1*deg2rad)*rchole/rcav)
  thetcav=acos(thetcav)

  thetatm5=pihalf-thettop
  thetatm3=pihalf-thetcav
  thetatm4=thetatm3+10.d0*deg2rad
  thetatm2=thetatm3-10.d0*deg2rad
  thetatm1=thetatm3/4.d0

  if (thetatm2.lt.thetatm1) then
     thetatm2=thetatm1
  end if

  if (thetatm4.ge.thetatm5) then
     thetatm5=thetatm3+(pihalf-thetatm3)/3.d0
     thetatm4=thetatm5
  endif

  ntatm5=nttmp-100
  ntatm4=ntatm5-int((thetatm5-thetatm4)/deg2rad/res)
  ntatm3=ntatm4-50
  ntatm1=(ntatm3-int((thetatm2-thetatm1)/deg2rad/res))/2
  ntatm2=ntatm1+int((thetatm2-thetatm1)/deg2rad/res)

  if (ifprint) then
     print*,''
     print*,'there are', nttmp-ntatm5 ,'cells from thet=', 0.d0 ,'deg to', (pihalf-thetatm5)/deg2rad, 'deg.'
     print*,'there are', ntatm5-ntatm4 ,'cells from thet=', (pihalf-thetatm5)/deg2rad,'deg to', (pihalf-thetatm4)/deg2rad, 'deg.'
     print*,'there are', ntatm4-ntatm3 ,'cells from thet=', (pihalf-thetatm4)/deg2rad,'deg to', (pihalf-thetatm3)/deg2rad, 'deg.'
     print*,'there are', ntatm3-ntatm2 ,'cells from thet=', (pihalf-thetatm3)/deg2rad,'deg to', (pihalf-thetatm2)/deg2rad, 'deg.'
     print*,'there are', ntatm2-ntatm1 ,'cells from thet=', (pihalf-thetatm2)/deg2rad,'deg to', (pihalf-thetatm1)/deg2rad, 'deg.'
     print*,'there are', ntatm1-1 ,'cells from thet=', (pihalf-thetatm1)/deg2rad,'deg to', 90.d0, 'deg.'
     print*,'there are totally',nttmp-1,'cells from thet=0 to 90 deg.'
     print*,''
  end if

  tmptharr(1)=0.d0
  
  gspl=2.d0
  do it=2,ntatm1-1
     tmptharr(it)=thetatm1*(dble(it-1)/dble(ntatm1-1))**gspl
  end do

  tmptharr(ntatm1)=thetatm1
  
  gspl=1.d0
  do it=ntatm1+1,ntatm2-1
     tmptharr(it)=thetatm1+(thetatm2-thetatm1)*(dble(it-ntatm1)/dble(ntatm2-ntatm1))**gspl
  end do

  tmptharr(ntatm2)=thetatm2

  gspl=3.5d0
  do it=ntatm2+1,ntatm3-1
     tmptharr(it)=asin(sin(thetatm3)-(sin(thetatm3)-sin(thetatm2))*(dble(ntatm3-it)/dble(ntatm3-ntatm2))**gspl)
  end do
  
  tmptharr(ntatm3)=thetatm3
  
  gspl=2.d0
  do it=ntatm3+1,ntatm4-1
     tmptharr(it)=thetatm3+(thetatm4-thetatm3)*(dble(it-ntatm3)/dble(ntatm4-ntatm3))**gspl
  end do
  
  tmptharr(ntatm4)=thetatm4
  
  gspl=1.d0
  do it=ntatm4+1,ntatm5-1
     tmptharr(it)=thetatm4+(thetatm5-thetatm4)*(dble(it-ntatm4)/dble(ntatm5-ntatm4))**gspl
  end do

  tmptharr(ntatm5)=thetatm5

  gspl=2.d0
  do it=ntatm5+1,nttmp-1
     tmptharr(it)=pihalf-(pihalf-thetatm5)*(dble(nttmp-it)/dble(nttmp-ntatm5))**gspl
  end do
  
  tmptharr(nttmp)=pihalf
  
  do it=2,nttmp-1
     thetarr(nttmp+1-it)=pihalf-tmptharr(it)
     thetarr(nttmp-1+it)=pihalf+tmptharr(it)
  end do
  thetarr(1)=0.d0
  thetarr(nttmp)=pihalf
  thetarr(ntg)=pi
  
  do it=1,size(thete_arr)
     if(minval(abs(thete_arr(it) - thetarr))==0.) then
        print *,'ERROR: peeloff angle is exactly the same as one of the wall angles'
        print *,'       and this is not recommended. The angle is ',thete_arr(it) * rad2deg
        stop
     end if
  end do
  
  ! this is not the eps used in the rest of the code!  see
  ! newdisktherm for that
  eps=1.d-8
  do it=1,ntg
     if (it.eq.1) then
        costarr(it)=1.d0
        sintarr(it)=0.d0
        tan2arr(it)=0.d0
     else if (it.eq.ntg) then
        costarr(it)=-1.d0
        sintarr(it)=0.d0
        tan2arr(it)=0.d0
     else if (it.eq.nttmp) then
        costarr(it)=0.d0
        sintarr(it)=1.d0
        tan2arr(it)=-1.d0
     else
        costarr(it)=cos(thetarr(it))
        sintarr(it)=sin(thetarr(it))
        tan2arr(it)=tan(thetarr(it))**2.d0
     end if
     ! print*,'thetarr,costarr,tan2arr '
     ! print*,'thetarr',thetarr(it)*rad2deg
  end do
  ! else
  ! thetarr(1)=0.d0
  ! costarr(1)=1.d0
  ! sintarr(1)=0.d0
  ! tan2arr(1)=0.d0
  ! end if

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  ! theta-ave array
  do it=1,ntg-1
     
     ! The following check is needed because
     ! ifort -m32 produces a discontinuity
     if(thetarr(it+1)<thetarr(it)) then
        stop "Error - thetarr is not monotonically increasing"
     end if
     
     thetavearr(it)=0.5d0*(thetarr(it)+thetarr(it+1))
     
  end do

  do it=1,ntg-1
     if (thetarr(it+1)-thetarr(it).le.0.d0) then
        print*,'thetarr wrong!!!',it,thetarr(it+1),thetarr(it)
        stop
     end if
  end do

  if (ifprint) then
     open(unit=15,file='tharr.dat',status='unknown')
     write(15,'(I6," / number of grid points in theta")') ntg
     write(15,'(" index    thet(rad)    thet(deg)      '// &
          & 'cost      tan**2(thet)")')
     do it=1,ntg
        write(15,'(i5,4(1x,f12.5))') it,thetarr(it), &
             & thetarr(it)*rad2deg,costarr(it),tan2arr(it)
     end do
     close(15)
  end if
  
  dp=r2p/dble(npg-1)
  do ip=1,npg
     phiarr(ip)=dp*(dble(ip-1))
     aarr(ip)=dsin(phiarr(ip))
     barr(ip)=-dcos(phiarr(ip))
     ! print*,'phiarr, aarr, barr ',phiarr(ip),aarr(ip),barr(ip)
     ! carr(ip)=
     ! darr(ip)=
  end do
  
  ! phi-ave array
  do ip=1,npg-1
     phiavearr(ip)=0.5d0*(phiarr(ip)+phiarr(ip+1))
  end do
  
  if (ifprint) then
     open(unit=15,file='phiarr.dat',status='unknown')
     write(15,'(I6," / number of grid points in phi")') npg
     write(15,'(" index    phi(rad)    phi(deg)")')
     do ip=1,npg
        write(15,'(i5,2(1x,f12.5))') ip,phiarr(ip), &
             & phiarr(ip)*rad2deg
     end do
     close(15)
  end if
  
  
  if(diffusion) then
     ! set up diffusion grid
     dcount=0
     diffdir=0
     if (massd.eq.0.d0) then
        diffus = .false.
        diffdir = 0
     else
        
        if (rddust.eq.rmin) then
           nr0=1
        else
           nr0=3
           diffus(1,:,:) = .false.
           diffdir(1,:,:) = 0
        end if
        
        do ir=nr0,nrg-1
           
           r=ravearr(ir)/rddust
           Rsonr=1.d0/(rddust*r)
           Tdisk=min(1400.d0/11605.d0,max(3.d0/11605.d0, &
                & (Tstar/11605.d0) &
                & *(max(2.d0/3.d0*(Rsonr)**3, &
                & (asin(Rsonr)-(Rsonr)*sqrt(1.d0-(Rsonr)**2)))/pi &
                & )**0.25d0))
           
           sigmu=0.d0
           taumid=0.d0
           
           do id=1,ndg
              if(is_disk(id)) then
                 sigmu = sigmu + fmass(id)*mu0*r**(b(id)-1.d0)
                 taumid = taumid + fmass(id)*tauzrat*tauv*chiR(Tdisk,id)/r**(a(id)-b(id))
              end if
           end do
           
           if (taumid.gt.taudisk) then
              mudisk=sigmu*sqrt(2.d0)*ierfc(taudisk/taumid)
           else
              mudisk=0.d0
           end if

           do it=1,ntg-1

              mu=abs(cos(0.5d0*(thetarr(it)+thetarr(it+1))))
              tauRmu=taumid*erfc(mu/(sigmu*sqrt(2.0d0)))
              tauRr=tauRos*exp(-0.5d0*mu*mu/(sigmu*sigmu))* &
                   & (1.d0-r**a1)/(1.d0-(rmaxd/rddust)**a1)

              do ip=1,npg-1

                 if ((mu.lt.mudisk).and.(tauRr.gt.taudisk)) then
                    diffus(ir,it,ip) = .true.
                    dcount=dcount+1
                 else
                    diffus(ir,it,ip) = .false.
                 end if

                 if (tauRr.lt.tauRmu) then
                    diffdir(ir,it,ip) = -1
                 else
                    if (thetarr(it).lt.pihalf) then
                       diffdir(ir,it,ip) = -2
                    else
                       diffdir(ir,it,ip) = 2
                    end if
                 end if

              end do
           end do
        end do
     end if
  else
     diffus = .false.
     diffdir = 0
     dcount = 0
  end if

  if (ifprint) print*,'number of diffusion cells in grid',dcount
  call diskdat_write('dcount',dcount,'number of diffusion cells in grid')

!  call output_grid('diffuse',diffus)
!  call output_grid('diffdir',diffdir)

  end if
! rest part will be run every time


  !calculate fractal density variations
    if (ifractal.eq.1) then
      call diskdat_write('ifractal',ifractal,'adding fractal density variations')
      if (ifprint) print*,''
      if (ifprint) print*,'computing fractal density variations'
	  call set_random_seed(ifseed)
      nf=32
      fl=fractl
  !   fl=3.792
      icent=0
      const=0.
      fm1=1.
      fm2=1.
      fm3=1.
      fm4=1.
 !   for more info, see ../fractal/fractal_sph.f
      call fractal_sph(nf,fl,const,fm1,fm2,fm3,fm4,icent)
  !modify densvararr to have average value of 1, std dev = densratio
	  densav=0.d0
	  densav2=0.d0
	  rntot=(nrg-1)*(ntg-1)*(npg-1)
      do ir=1,nrg-1
	  do it=1,ntg-1
	  do ip=1,npg-1
		 densav=densvararr(ir,it,ip)+densav
		 densav2=densvararr(ir,it,ip)**2+densav2
	  enddo
      enddo
      enddo
      if (ifprint) print*,'rntot',rntot
      densav=densav/rntot
      densav2=densav2/rntot
      if (ifprint) print*,'densav2,densav**2',densav2,densav**2
      if (ifprint) print*,'std^2',densav2-densav**2
      std=sqrt(densav2-densav**2)
      if (ifprint) print*,'densav,std',densav,std
      if (ifprint) print*,'std/densav',std/densav
      if (ifprint) print*,'std',std
  !    call diskdat_write('std/densav',std/densav,'std dev of dens. vars normalized to ave')
      do ir=1,nrg-1
         do it=1,ntg-1
            do ip=1,npg-1
               !      densvararr(ir,it,ip)=densvararr(ir,it,ip)/densav
               densvararr(ir,it,ip)=densvararr(ir,it,ip)/densav*densratio+1-(densratio)
               !   densvararr(ir,it,ip)=(densvararr(ir,it,ip)-densav)*densratio/std+1.
               !   densvararr(ir,it,ip)=densvararr(ir,it,ip)*densratio/std
            enddo
         enddo
      enddo
      ! test; can comment this out when it works
      densav=0.d0
      densav2=0.d0
      mindensvar=0.d0
      maxdensvar=0.d0
   do ir=1,nrg-1
      do it=1,ntg-1
         do ip=1,npg-1
            densav=densvararr(ir,it,ip)+densav
            densav2=densvararr(ir,it,ip)**2+densav2
            if (densvararr(ir,it,ip).lt.0) print*,'densvararr lt 0',densvararr(ir,it,ip)
            if (densvararr(ir,it,ip).gt.maxdensvar) maxdensvar=densvararr(ir,it,ip)
            if (densvararr(ir,it,ip).lt.mindensvar) mindensvar=densvararr(ir,it,ip)
         enddo
      enddo
   enddo
   densav=densav/rntot
   densav2=densav2/rntot
   std=sqrt(densav2-densav**2)
   if (ifprint) print*,'fractal densav,std, min,max',densav,std,mindensvar,maxdensvar
   call diskdat_write('fractal dens ave',densav,'should be 1')
   call diskdat_write('fractal std',std,'standard deviation')
   call set_random_seed(i1)
   if (ifprint) print*,''
   endif

   compo=-1

!!$   call cpu_time(cpusec)
!!$   print*,''
!!$   print*,'CPU TIME',cpusec
!!$   print*,''

   if (ifprint) print*,'now opening rad file'
   close(1)
   open(1,file=tabname(1),status='old')
   read(1,*) indr
   if (ifprint) print*,'rad ind',indr
   allocate(radtab(indr))
   do i=1,indr
      read(1,*) radtab(i)
   end do
   close(1)

   if (ifprint) print*,'now opening theta file'
   close(1)
   open(1,file=tabname(2),status='old')
   read(1,*) indt
   if (ifprint) print*,'theta ind',indt
   allocate(thetatab(indt))
   do i=1,indt
      read(1,*) thetatab(i)
   end do
   close(1)

   if (ifprint) print*,'now opening density file'
   close(1)
   open(1,file=tabname(3),status='old')
   allocate(dentab(indt,indr))
   do i=1,indr
      do j=1,indt
         read(1,*) dentab(j,i)
      end do
   end do
   close(1)
!   print*,'(1,2),(indt,1)',dentab(1,2),dentab(indt,1)

   if (ifprint) print*,'now opening velocity file'
   close(1)
   open(1,file=tabname(4),status='old')
   read(1,*) indv
   if (ifprint) print*,'v ind',indv
   allocate(vtab(indv))
   do i=1,indv
      read(1,*) vtab(i)
   end do
   close(1)

  ! set scale factor for disk mass.  involves looping through and
  ! calculating mass twice.  :)
  ! 1 June 2009 - updated cell volume calculation to exact solution
!  if (massd.gt.0.d0) then
!     massdisk=0.d0
!     do id=1,ndg+2
!        if(is_disk(id)) then
!           do ir=1,nrg-1
!              rad=ravearr(ir)
!              dr=rarr(ir+1)-rarr(ir)
!              dr3=rarr(ir+1)**3.d0-rarr(ir)**3.d0
!              do it=1,ntg-1
!                 thet=0.5d0*(thetarr(it)+thetarr(it+1))
!                 cost=cos(thet)
!                 sint=sin(thet)
!                 dcost=cos(thetarr(it+1))-cos(thetarr(it))
!                 dt=thetarr(it+1)-thetarr(it)
!                 do ip=1,npg-1 !3-D atmosphere
!                    phi=0.5d0*(phiarr(ip)+phiarr(ip+1))
!                    dp=phiarr(ip+1)-phiarr(ip)
!                    vol = - dr3 * dcost * dp / 3.d0 * rmincgs**3.d0
!                    ! vol=rad**2*sint*dt*dp*dr*rmincgs**3
!                    call densdisk(rad,sint,cost,phi,densd,ir,it,ip,id)
!                    massdisk=massdisk+densd*vol
!                 end do
!              end do
!           end do
!        end if
!     end do
!     massdisk=massdisk/msun
!     diskscale=massd/massdisk
!  else
!     diskscale=1.d0
!  end if
!
!  print*,'massdisk before scaling (all grains)',massdisk
!  print*,'diskscale',diskscale
   
   call densdiskg_new(rmincgs,iter,lstar,lacc,lacc2)

!   call output_grid('darr_temp',densarr)

  ! calculate density in grid
  massenv=0.d0
!  massdisk=0.d0
  maxdens=0.d0
  masssg=0.d0
  totvol=0.d0

  x0arr=0.d0

!!$  ir_fp=0
!!$  do ir=1,nrg
!!$     t_fp=tdave(ir,(ntg-1)/2,1)*11605.d0
!!$!     print*,t_fp
!!$     if (t_fp.gt.2000.d0) then
!!$        ir_fp=ir
!!$     end if
!!$  end do
!!$
!!$  ir_fp=ir_fp+1
!!$  print*,'rad_fp',ravearr(ir_fp)/autors
!!$  it=(ntg-1)/2
!!$  do while (densarr(ir_fp,it,1,1)+densarr(ir_fp,it,1,2)+densarr(ir_fp,it,1,9).gt.0.d0)
!!$     it=it-1
!!$  end do
!!$  it_fp=it+1
!!$
!!$  rad=rarr(ir_fp)
!!$  thet=thetarr(it_fp)
!!$  cost=cos(thet)
!!$  sint=sin(thet)
!!$  call densdiskwind(rad*sint,rad*cost,dens1,x0,vz,vcyn)
!!$  x0_fp=x0

  ir_fp=1
  do ir=1,nrg
     thet=acos(zdisk(ir)/ravearr(ir))
     call locate(thetarr,ntg-1,thet,it)
     if (ravearr(ir)*cos(thetarr(it)).gt.zdisk(ir)) it=it+1
     t_fp=tdave(ir,it,1)*11605.d0
     if (t_fp.gt.1600.d0) then
        ir_fp=ir
     end if
  end do
  x0_fp=ravearr(ir_fp)
  x0_fp=log(x0_fp)/log(x_max0)

  if (ifprint) print*,'ir_fp,r_fp,x0_fp',ir_fp,ravearr(ir_fp)/autors,x0_fp

!!$  call cpu_time(cpusec)
!!$  print*,''
!!$  print*,'CPU TIME',cpusec
!!$  print*,''

  do ir=1,nrg-1

     rad=ravearr(ir)
     dr=rarr(ir+1)-rarr(ir)
     dr3=rarr(ir+1)**3.d0-rarr(ir)**3.

     !radtab/thetatab is sortted decreasely
     call locate(radtab,indr-1,rad,indru)
     if (indru.eq.0) indru=1
     indrl=indru+1
     upr=radtab(indru)
     lor=radtab(indrl)
!    if (ir.eq.1) print*,upr,lor
     vu=vtab(indru)
     vl=vtab(indrl)
     vew(ir)=(vu-vl)/(upr-lor)*(rad-lor)+vl

     do it=1,ntg-1

        thet=0.5d0*(thetarr(it)+thetarr(it+1))
        cost=cos(thet)
        sint=sin(thet)
        dt=thetarr(it+1)-thetarr(it)
        dcost=cos(thetarr(it+1))-cos(thetarr(it))

        call locate(thetatab,indt-1,thet,indtu)
        if (indtu.eq.0) indtu=1
        indtl=indtu+1
        upt=thetatab(indtu)
        lot=thetatab(indtl)
        
        denll=dentab(indtl,indrl)
        denlu=dentab(indtl,indru)
        denul=dentab(indtu,indrl)
        denuu=dentab(indtu,indru)

!        if (ir.eq.1) print*,denll,denlu,denul,denuu
!        print*,lor,upr,lot,upt

        do ip=1,npg-1 ! 3-D atmosphere

           phi=0.5d0*(phiarr(ip)+phiarr(ip+1))
           dp=phiarr(ip+1)-phiarr(ip)

           vol = - dr3 * dcost * dp / 3.d0 * rmincgs**3.d0
           totvol=totvol+vol
           vcell(ir,it,ip)=vol

           tdust_temp=tdave(ir,it,ip)*11605.d0

           ! Disk grains
           do id=1,ndg+2

!              if(is_disk(id)) then
!                 call densdisk(rad,sint,cost,phi,densd,ir,it,ip,id)
!                 densarr(ir,it,ip,id)=densarr(ir,it,ip,id)+densd*diskscale ! densd already includes fmassd scale factor in rho0
!                 densarr(ir,it,ip,id)=densd*diskscale
!                 massdisk=massdisk+densd*vol*diskscale
!              end if

              if(is_envelope(id) .or. is_cavity(id)) then
!                 if (ienvtype.eq.1) then
!                    call densenv(rad,thet,cost,sint,phi,dense,radamb,id,inhole)
!                 else
!                    call densenvpw(rad,thet,cost,sint,phi,dense,radamb,id,inhole,rmincgs)
!                 endif

                 call densenv_read(rad,thet,cost,sint,phi,lor,upr,lot,upt,denll,denuu,denlu,denul,dense,radamb,id,inhole,x0,vz,vcyn,ir,it)
                 x0arr(ir,it,ip)=x0
!                 if (ir.eq.1) then
!                 print*,dense
!                 end if

!                 if (id.eq.10) then
!                    densarr(ir,it,ip,id)=dense
!                 else
!                    densarr(ir,it,ip,id)=dense*fmass(id)
!                 end if
!                 massenv=massenv+dense*vol*fmass(id)

                 if (is_envelope(id)) then
                    densarr(ir,it,ip,id)=dense*fmass(id)
                    massenv=massenv+dense*vol*fmass(id)
                 end if

                 if (is_cavity(id)) then

                    if (id.eq.4.or.id.eq.8) then

                       if (x0.lt.x0_fp) then
                          densarr(ir,it,ip,id)=0.d0
                       else
                          if (tdust_temp.gt.1400.d0) then
                             densarr(ir,it,ip,id)=0.d0
                          else
                             densarr(ir,it,ip,id)=dense*fmass(id)
                             massenv=massenv+dense*vol*fmass(id)
                          end if
                       end if
                       
                    else

                       if (x0.lt.x0_fp) then
                          densarr(ir,it,ip,id)=dense
                          massenv=massenv+dense*vol
                       else
                          if (tdust_temp.ge.1400.d0) then
                             densarr(ir,it,ip,id)=dense
                             massenv=massenv+dense*vol
                          else
                             densarr(ir,it,ip,id)=0.d0
                          end if
                       end if

                    end if

                 end if

              end if
              
              if (ifractal.eq.1.and.fractmask(id).eq.1) then
                densarr(ir,it,ip,id)=densarr(ir,it,ip,id)*(densvararr(ir,it,ip))
     !       densarr(ir,it,ip,id)=densarr(ir,it,ip,id)*(densvararr(ir,it,ip))
              endif

           end do

           if (densarr(ir,it,ip,1)+densarr(ir,it,ip,2)+densarr(ir,it,ip,5)+densarr(ir,it,ip,6)+densarr(ir,it,ip,9).gt.0.d0) then
              densarr(ir,it,ip,3)=0.d0
              densarr(ir,it,ip,7)=0.d0
              densarr(ir,it,ip,4)=0.d0
              densarr(ir,it,ip,8)=0.d0
              densarr(ir,it,ip,10)=0.d0
              compo(ir,it,ip)=0 ! for disk
           else
              if (inhole.eq.0) then
                 compo(ir,it,ip)=1 ! for envelope
              else 
                 compo(ir,it,ip)=2 ! for outflow
              end if
           end if           

           do id=1,ndg+2
!              massarr(ir,it,ip,id)=(densarr(ir,it,ip,id)*vol)**0.25d0 ! bad programming:  massarr is really mass**0.25
              if(is_sg(id)) masssg=masssg+densarr(ir,it,ip,id)*vol
           end do

           if (sum(densarr(ir,it,ip,:)).eq.0.d0) then
              diffus(ir,it,ip)=.false.
              diffdir(ir,it,ip)=0
           end if

           tmp = 0.d0
           do id=1,ndg+2
              tmp=tmp+densarr(ir,it,ip,id)
           end do

           if (tmp.gt.maxdens) maxdens=tmp

!!$           tdust_temp=0.d0
!!$           dcount=0
!!$           do id=1,ndg+2
!!$              tdust_temp=max(tdust_temp,tdust2(ir,it,ip,id))
!!$              if (tdust2(ir,it,ip,id).gt.0.1d0/11605.d0) dcount=dcount+1
!!$           end do
!!$           if (dcount.gt.1) then
!!$              print*,'more than one real temperature at',ir,it,ip
!!$              stop
!!$           end if
!!$
!!$           dcount=0
!!$           do id=1,ndg+2
!!$              if (densarr(ir,it,ip,id).eq.0.d0) then
!!$                 tdust2(ir,it,ip,id)=0.1d0/11605.d0
!!$              else
!!$                 dcount=dcount+1
!!$                 tdust2(ir,it,ip,id)=tdust_temp
!!$              end if
!!$           end do
!!$           if (dcount.gt.1) then
!!$              print*,'more than one non-zero density at',ir,it,ip
!!$              stop
!!$           end if
              
!           if (.not.diffus(ir,it,ip)) then
!              do id=1,ndg+2
!                 tdust(ir,it,ip,id)=tdust2(ir,it,ip,id)
!              end do
!           end if

        end do
     end do
  end do

  if (ifprint) print*,'density grid setup finished in process'

  vzarr=0.d0
  vcynarr=0.d0
  vphiarr=0.d0

  vrarr=0.d0
  vthetarr=0.d0
  divvarr=0.d0

  vrarru=0.d0
  vthetarru=0.d0
  divvarru=0.d0

  vrarre=0.d0
  vthetarre=0.d0
  divvarre=0.d0

  vrarrw=0.d0
  vthetarrw=0.d0
  vphiarrw=0.d0
  divvarrw=0.d0

  do ir=1,nrg
     rad=rarr(ir)
     do it=1,ntg
        thet=thetarr(it)
        cost=abs(cos(thet))
        sint=sin(thet)
        call densdiskwind(rad*sint,rad*cost,dens1,x0,vz,vcyn,vphi,ir,it)
        if (cos(thet).gt.0.d0) then
           vzarr(ir,it,:)=vz 
        else
           vzarr(ir,it,:)=-vz
        end if
        vcynarr(ir,it,:)=vcyn
        vrarrw(ir,it,:)=vz*cost+vcyn*sint
        if (cos(thet).gt.0.d0) then
           vthetarrw(ir,it,:)=vcyn*cost-vz*sint
        else
           vthetarrw(ir,it,:)=-vcyn*cost+vz*sint
        end if
        vphiarrw(ir,it,:)=vphi
     end do
  end do

  !div in spherical coordinates
  do ir=1,nrg-1
     rad=ravearr(ir)
     rad1=rarr(ir+1)
     rad2=rarr(ir)
     dr=rad1-rad2
     do it=1,(ntg-1)/2
        theta1=thetarr(it+1)
        theta2=thetarr(it)
        thet=0.5d0*(theta1+theta2)
        dt=theta1-theta2
        vr1=0.5d0*(vrarrw(ir+1,it+1,1)+vrarrw(ir+1,it,1))
        vr2=0.5d0*(vrarrw(ir,it+1,1)+vrarrw(ir,it,1))
        vthet1=0.5d0*(vthetarrw(ir+1,it+1,1)+vthetarrw(ir,it+1,1))
        vthet2=0.5d0*(vthetarrw(ir+1,it,1)+vthetarrw(ir,it,1))
        divvarrw(ir,it,:)=(rad1**2*vr1-rad2**2*vr2)/dr/rad**2 &
             & +(sin(theta1)*vthet1-sin(theta2)*vthet2)/dt/rad/sin(thet)
     end do
  end do

  do ir=1,nrg-1
     do it=(ntg+1)/2,ntg-1
        divvarrw(ir,it,:)=divvarrw(ir,ntg-it,1)
     end do
  end do

  do ir=1,nrg-1
     rad=ravearr(ir)
     rad1=rarr(ir+1)
     rad2=rarr(ir)
     dr=rad1-rad2
     vfall=sqrt(gn*(massc+massd)*msol/rad/rstar/rsol)
     do it=1,ntg
        thet=thetarr(it)
        sint=sin(thet)
        cost=abs(cos(thet))
        if (it.eq.(ntg+1)/2.and.rad.gt.rchole) then
           vrarru(ir,it,:)=-vfall*sqrt(2.d0-rchole/rad)
           vthetarru(ir,it,:)=0.d0
        else if (it.eq.(ntg+1)/2.and.rad.le.rchole) then
           cost0=sqrt(1.d0-rad/rchole)
           vrarru(ir,it,:)=-vfall
           vthetarru(ir,it,:)=vfall*cost0
        else
           call zerod(cost,rad/rchole,cost0,iflag)
           vrarru(ir,it,:)=-vfall*sqrt(1.d0+cost/cost0)
           vthetarru(ir,it,:)=vfall*(cost0-cost)/sint*sqrt(1.d0+cost/cost0)
        end if
     end do
  end do

  do ir=1,nrg-1
     vrarre(ir,:,:)=-vew(ir)
     vthetarre(ir,:,:)=0.d0
  end do

  do ir=1,nrg-1
     rad=ravearr(ir)
     rad1=rarr(ir+1)
     rad2=rarr(ir)
     dr=rad1-rad2
     do it=1,(ntg-1)/2
        theta1=thetarr(it+1)
        theta2=thetarr(it)
        thet=0.5d0*(theta1+theta2)
        dt=theta1-theta2
        vr1=0.5d0*(vrarru(ir+1,it+1,1)+vrarru(ir+1,it,1))
        vr2=0.5d0*(vrarru(ir,it+1,1)+vrarru(ir,it,1))
        vthet1=0.5d0*(vthetarru(ir+1,it+1,1)+vthetarru(ir,it+1,1))
        vthet2=0.5d0*(vthetarru(ir+1,it,1)+vthetarru(ir,it,1))
        divvarru(ir,it,:)=(rad1**2*vr1-rad2**2*vr2)/dr/rad**2 &
             & +(sin(theta1)*vthet1-sin(theta2)*vthet2)/dt/rad/sin(thet)
     end do
  end do

  do ir=1,nrg-1
     rad=ravearr(ir)
     rad1=rarr(ir+1)
     rad2=rarr(ir)
     dr=rad1-rad2
     do it=1,(ntg-1)/2
        theta1=thetarr(it+1)
        theta2=thetarr(it)
        thet=0.5d0*(theta1+theta2)
        dt=theta1-theta2
        vr1=0.5d0*(vrarre(ir+1,it+1,1)+vrarre(ir+1,it,1))
        vr2=0.5d0*(vrarre(ir,it+1,1)+vrarre(ir,it,1))
        vthet1=0.5d0*(vthetarre(ir+1,it+1,1)+vthetarre(ir,it+1,1))
        vthet2=0.5d0*(vthetarre(ir+1,it,1)+vthetarre(ir,it,1))
        divvarre(ir,it,:)=(rad1**2*vr1-rad2**2*vr2)/dr/rad**2 &
             & +(sin(theta1)*vthet1-sin(theta2)*vthet2)/dt/rad/sin(thet)
     end do
  end do

  do ir=1,nrg-1
     do it=(ntg+1)/2,ntg-1
        divvarru(ir,it,:)=divvarru(ir,ntg-it,1)
        divvarre(ir,it,:)=divvarre(ir,ntg-it,1)
     end do
  end do

  do ir=1,nrg-1
     do it=1,ntg-1
        theta1=thetarr(it+1)
        theta2=thetarr(it)
        thet=0.5d0*(theta1+theta2)
        sint=sin(thet)
        cost=abs(cos(thet))
        if (densarr(ir,it,1,4)+densarr(ir,it,1,8)+densarr(ir,it,1,10).gt.0.d0) then
           vrarr(ir,it,:)=vrarrw(ir,it,1)
           vthetarr(ir,it,:)=vthetarrw(ir,it,1)
           divvarr(ir,it,:)=divvarrw(ir,it,1)
           r_foot=(rmaxd*(1.d0-windmu0**2))**x0arr(ir,it,1)
!           vphiarr(ir,it,:)=sqrt(gn*massc*msol/r_foot/rstar/rsol)*r_foot/ravearr(ir)/sin(thetavearr(ir))
           vphiarr(ir,it,:)=vphiarrw(ir,it,1)
        else if (densarr(ir,it,1,3)+densarr(ir,it,1,7).gt.0.d0) then
           if (ravearr(ir).lt.sonicp*autors) then
              vrarr(ir,it,:)=vrarru(ir,it,1)
              vthetarr(ir,it,:)=vthetarru(ir,it,1)
              divvarr(ir,it,:)=divvarru(ir,it,1)
              vphiarr(ir,it,:)=sqrt(gn*massc*4.d0/3.d0*msol/ravearr(ir)/rstar/rsol)*sqrt(1.d0-windmu0**2) &
                   & /sint*sqrt(1.d0-cost/windmu0)
           else
              vrarr(ir,it,:)=vrarre(ir,it,1)
              vthetarr(ir,it,:)=vthetarre(ir,it,1)
              divvarr(ir,it,:)=divvarre(ir,it,1)
              vphiarr(ir,it,:)=0.d0
           end if
        else
           vrarr(ir,it,:)=0.d0
           vthetarr(ir,it,:)=0.d0
           divvarr(ir,it,:)=0.d0
           vphiarr(ir,it,:)=sqrt(gn*massc*msol/ravearr(ir)/rstar/rsol)
!           vphiarr(ir,it,:)=0.d0
        end if
     end do
  end do

!  tdust2=tdust2*11605.d0

  call findopacid()

  write(suffix,'("_",I3.3)') iter

  if (ifprint) then
     call output_grid('x0arr',x0arr)
     !  call output_grid('vzarr',vzarr)
     !  call output_grid('vcynarr',vcynarr)
     call output_grid('vrarr',vrarr)
     call output_grid('vthetarr',vthetarr)
     call output_grid('vphiarr',vphiarr)
     !  call output_grid('vrarru',vrarru)
     !  call output_grid('vthetarru',vthetarru)
     !  call output_grid('vrarrw',vrarrw)
     !  call output_grid('vthetarrw',vthetarrw)
     !  call output_grid('vrarre',vrarre)
     !  call output_grid('vthetarre',vthetarre)
     !  call output_grid('divvarr',divvarr)
  end if
     
!  tdust2=tdust2/11605.d0

!  print*,densarr(275,499,1,1)
  maxdens = maxval(sum(densarr,dim=4))

  tdust_temp=0.d0
  do ir=1,nrg
     do it=1,ntg
        do ip=1,npg
           tdust_temp=max(tdust_temp,tdave(ir,it,ip))
        end do
     end do
  end do
  if (ifprint) print*,'highest temperature of the gas in disk', tdust_temp*11605.d0

!!$  call cpu_time(cpusec)
!!$  print*,''
!!$  print*,'CPU TIME',cpusec
!!$  print*,''

  if (ifprint) print*,'z1,b,rho0',z1,b,rho0
  massenv=massenv/msun
  massdisk=massdisk/msun
  masssg=masssg/msun
  ! important note:  densarr(ir,it,ip) is the density halfway between
  ! ir and ir+1, it and it+1, ip and ip+1; it is NOT the density at
  ! ir,it,ip.
  if (ifprint) print*, 'massenv (all grains) ',massenv
  if (ifprint) print*, 'massdisk from grid, compared to input ',massdisk, &
       & massd
  ! masssg is actually mass sg+20-200 A grains
  if (ifprint) then
     print*, 'mass very small grains (non-thermal) (solar) ', masssg
     print*,'massvsgs/masstot',masssg/(massdisk+massenv)
     print*,'max density in disk/envelope',maxdens
     print*,'min radius where rho=rhoamb ',radamb/autors
     print*, ' '
     print*,'totvol',totvol
  end if

  ! write out rarr,thetarr,densarr at phi=0 so you can read them into
  ! IDL.

  ! Av arr
  do ip=1,npg-1
     do it=1,ntg-1
        tauave=0.d0
        do ir=1,nrg-1
           dr=rarr(ir+1)-rarr(ir)
           rad=ravearr(ir)
           do id=1,ndg+2
              if (densarr(ir,it,ip,id).gt.0.d0) then
                 idust2=findopac(ir,it,ip,id)
                 tauave=tauave+kapd(idust2)*densarr(ir,it,ip,id)*dr
!                 if (ifprint) print*,tauave,idust2,ir,it,ip,id,findopac(ir,it,ip,id)
              end if
           end do
           avarr(ir,it,ip)=tauave*1.086d0
        end do
     end do
  end do
  if (ifprint) call output_grid('avarr',avarr)

  ! Av arr along the theta direction (estimate Av from outflow cavity wall)
  do ip=1,npg-1
     do ir=1,nrg-1
        rad=ravearr(ir)
        tauave=0.d0
        do it=1,ntg-1
           dr=rad*(thetarr(it+1)-thetarr(it))
           do id=1,ndg+2
              if (densarr(ir,it,ip,id).gt.0.d0) then
                 idust2=findopac(ir,it,ip,id)
                 tauave=tauave+kapd(idust2)*densarr(ir,it,ip,id)*dr
              end if
           end do
           avthetarr(ir,it,ip)=tauave*1.086d0
        end do
     end do
  end do
  if (ifprint) call output_grid('avthetarr',avthetarr)

  if (ifprint) then
  open(unit=15,file='Av_view.dat',status='unknown')
  ! calculate average optical depth (A_v) along all theta & phi directions
  write(15,*) &
       & 'Av along viewing directions'
  write(15,*) &
       & 'theta       phi      Av_env    Av_disk      Av_sg     Av_tot'
  tauave=0.d0
  dcount=0
  dphi=r2p/dble(nph)
  do i=1,nmu
     thet=acos(u(i))
     if (thet.eq.0.d0) thet=1.d-5
     sint=sin(thet)
     call locate(thetarr,ntg,thet,it)
     if (ntg.eq.1) it=1
     ! print*,'it',it
     do j=1,nph
        phi=dble(j)*dphi-0.5d0*dphi
        ! print*,'phi',phi
        if (phi.gt.r2p) print*,'oops!  phi',phi
        call locate(phiarr,npg,phi,ip)
        tau=0.d0
        taud=0.d0
        taue=0.d0
        taup=0.d0
        cost=cos(thet)
        do ir=1,nrg-1
           dr=rarr(ir+1)-rarr(ir)
           rad=ravearr(ir)
           ! print*,'ir,it,ip',ir,it,ip
           do id=1,2
              taud=taud+kapd(id)*densarr(ir,it,ip,id)*dr
           end do
           do id=3,4
              taue=taue+kapd(id)*densarr(ir,it,ip,id)*dr
           end do
           do id=5,8
              taup=taup+kapd(id)*densarr(ir,it,ip,id)*dr
           end do
        end do
        write(15,901) thet*rad2deg,phi*rad2deg &
             & ,taue*1.086d0,taud*1.086d0,taup*1.086d0 &
             & ,(taue+taud+taup)*1.086d0
901     format((f10.5,1x,f10.5),5(1x,1pe12.5))
        tauave=tauave+taud+taue+taup
        dcount=dcount+1
     end do
  end do
  tauave=tauave/dble(dcount)
  print*,'dcount',dcount
  print*,'A_V average over all directions ',tauave*1.086d0
  write(15,*) 'A_V average over all directions ',tauave*1.086d0
  call diskdat_write('Av_ave',tauave*1.086d0, &
       & 'A_V ave over all directions')
  close(15)
  end if

  ! stop

  ! integrate optical depth along all the grid angle directions
  tauave=0.d0
  dcount=0
  if (ifprint) then
  open(unit=15,file='Av_grid.dat',status='unknown')
  write(15,*) 'Av through each grid theta bin '
  write(15,*) &
       & 'theta    phi    AV_env       AV_disk    AV_sg      AV_tot'
  do ip=1,npg-1
     phi=phiavearr(ip)
     do it=1,ntg-1
        tau=0.d0
        taud=0.d0
        taue=0.d0
        taup=0.d0
        thet=thetavearr(it)
        if (thet.gt.pi/2.d0-eps.and.thet.lt.pi/2.d0+eps) then
           ! TSC blows up at thet=90
           print*,'thet',thet*rad2deg
           thet=89.999d0*deg2rad
        end if
        cost=cos(thet)
        sint=sin(thet)
        phi=0.d0
        do ir=1,nrg-1
           dr=rarr(ir+1)-rarr(ir)
           ! rad=0.5d0*(rarr(ir)+rarr(ir+1))
           rad=ravearr(ir)
           do id=1,2
              taud=taud+kapd(id)*densarr(ir,it,ip,id)*dr*1.086d0
           end do
           do id=3,4
              taue=taue+kapd(id)*densarr(ir,it,ip,id)*dr*1.086d0
           end do
           do id=5,8
              taup=taup+kapd(id)*densarr(ir,it,ip,id)*dr*1.086d0
           end do
        end do
        write(15,901) thet*rad2deg,phi*rad2deg &
             & ,taue &
             & ,taud,taup,taue+taud+taup
        tauave=tauave+taud+taue+taup
        dcount=dcount+1
     end do
  end do
  tauave=tauave/dble(dcount)
  print*,'dcount',dcount
  print*,'A_V average over all directions ',tauave
  call diskdat_write('Av_ave',tauave, &
       & 'A_V ave over all directions')
  write(15,*) 'A_V average over all directions ',tauave
  close(15)
  end if
  ! 902     format(f10.5,4(1x,f15.5))

  ! calculate Av along ethet,ephi directions (peeled image)
  if (ifprint) then
  open(unit=15,file='Av_peel.dat',status='unknown')
  write(15,*) 'Av along peeled direction ethet'
  write(15,*) &
       & 'theta     phi      Av_env     Av_disk    AV_sg     Av_tot'
  do peelid=1,npeel
     thet=thete_arr(peelid)
     cost=coste_arr(peelid)
     sint=sinte_arr(peelid)
     call locate(thetarr,ntg,thet,it)
     phi=phie
     call locate(phiarr,npg,phi,ip)
     taud=0.d0
     taue=0.d0
     taup=0.d0
     do ir=1,nrg-1
        dr=rarr(ir+1)-rarr(ir)
        ! rad=0.5d0*(rarr(ir)+rarr(ir+1))
        rad=ravearr(ir)
        do id=1,2
           taud=taud+kapd(id)*densarr(ir,it,ip,id)*dr*1.086d0
        end do
        do id=3,4
           taue=taue+kapd(id)*densarr(ir,it,ip,id)*dr*1.086d0
        end do
        do id=5,8
           taup=taup+kapd(id)*densarr(ir,it,ip,id)*dr*1.086d0
        end do
     end do
     write(15,901) thet*rad2deg,phi*rad2deg &
          & ,taue &
          & ,taud,taup,taue+taud+taup
     write(15,*) 'average A_v along *all* directions ', &
          & tauave
  end do
  close(15)
  end if
  
  ! calculate Av along the disk midplane
  if (ifprint) then
  open(unit=15,file='Av_90.dat',status='unknown')
  write(15,*) 'Av in the disk midplane'
  write(15,*) &
       & 'theta     phi      Av_env     Av_disk    AV_sg     Av_tot'
     thet=pi/2.d0
     cost=0.d0
     sint=1.d0
     call locate(thetarr,ntg,thet,it)
     phi=0.d0
     call locate(phiarr,npg,phi,ip)
     taud=0.d0
     taue=0.d0
     taup=0.d0
     taug=0.d0
     do ir=1,nrg-1
        dr=rarr(ir+1)-rarr(ir)
        ! rad=0.5d0*(rarr(ir)+rarr(ir+1))
        rad=ravearr(ir)
        do id=1,2
           taud=taud+kapd(id)*densarr(ir,it,ip,id)*dr*1.086d0
        end do
        do id=3,4
           taue=taue+kapd(id)*densarr(ir,it,ip,id)*dr*1.086d0
        end do
        do id=5,8
           taup=taup+kapd(id)*densarr(ir,it,ip,id)*dr*1.086d0
        end do
        do id=9,9
           taug=taug+kapd(id)*densarr(ir,it,ip,id)*dr*1.086d0
        end do
     end do
     write(15,901) thet*rad2deg,phi*rad2deg &
          & ,taue &
          & ,taud,taup,taue+taud+taup
     write(15,*) 'average A_v along *all* directions ', &
          & tauave
  close(15)
  end if

  ! calculate tau_V=1 surf from origin (not quite the stellar surface!).
  if (ifprint) then
  open(unit=15,file='taudsurfin.dat',status='unknown')
  write(15,*) 'cyl r, z values for tau=1 surface'
  do it=1,ntg-1
     ! thet=(90.-(45./2000.*(it-1)))*deg2rad
     ! call locate(thetarr,ntg,thet,it)
     taud=0.d0
     cost=cos(thet)
     sint=sin(thet)
     phi=0.d0
     ip=1
     ir=0
     do while (taud.lt.1.d0.and.ir.lt.(nrg-1))
        ir=ir+1
        dr=rarr(ir+1)-rarr(ir)
        rad=0.5d0*(rarr(ir)+rarr(ir+1))
        do id=1,2
           taud=taud+kapd(id)*densarr(ir,it,ip,id)*dr
        end do
     end do
     x=rad*sint
     z=rad*cost
     write(15,*) x/autors,z/autors
  end do
  close(15)
  end if

  ! r=1.d0*autors
  ! call locate(rarr,nrg,dble(r),ir)
  ! dens1=densarr(ir,1,1)
  ! print*,'r,dens1',r,dens1

  idmax=0
  do ir=1,nrg
     do it=1,ntg
        do ip=1,npg
           do id=1,ndg+2
              idmax=max(idmax,findopac(ir,it,ip,id))
           end do
        end do
     end do
  end do

  idmin=idmax
  do ir=1,nrg
     do it=1,ntg
        do ip=1,npg
           do id=1,ndg+2
              idmin=min(idmin,findopac(ir,it,ip,id))
           end do
        end do
     end do
  end do

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  print*,'max and min in findopac array',myid,idmax,idmin
!  stop
  
  if (diffusion.and.iter.ge.1) call setdiffus(iter)

!!$  call cpu_time(cpusec)
!!$  print*,''
!!$  print*,'CPU TIME',cpusec
!!$  print*,''
  
!  do id=1,ndg+2
!     where(diffus)
!        tdust(:,:,:,id)=0.1d0/11605.d0
!     elsewhere
!        tdust(:,:,:,id)=tdust2(:,:,:,id)
!     end where
!  end do

  do ir=1,nrg
     do it=1,ntg
        do ip=1,npg
           if (diffus(ir,it,ip)) then
              tdust(ir,it,ip)=0.1d0/11605.d0
           else
              tdust(ir,it,ip)=tdave(ir,it,ip)
           end if
           if (sum(densarr(ir,it,ip,:)).eq.0.d0) tdust(ir,it,ip)=0.1d0/11605.d0
        end do
     end do
  end do

!  call output_grid('opacid'//suffix,findopac)

  ! write out massarr
!  call output_grid('marr'//suffix,massarr**4)

  ! write out densarr
!  call output_grid('darr'//suffix,densarr)
  if (ifprint) call output_grid('darr',densarr)

  ! oops, not f77 compatible
  ! masstot = sum(densarr) / msun
  ! massdisk = sum(densarr(:,:,:,1)*vcell + densarr(:,:,:,2)*vcell) / msun
  ! massenv = sum(densarr(:,:,:,3)*vcell) / msun
  ! masssgreg = sum(densarr(:,:,:,4)*vcell) / msun

  masstot=massenv+massdisk+masssg

  call diskdat_write('masstot',masstot,'total circumstellar mass')
  call diskdat_write('massenv',massenv,'envelope mass without sgs')
  call diskdat_write('massdisk',massdisk,'disk mass without sgs')
  call diskdat_write('masssgs',masssg,'total  mass, small grains')

  call cpu_time(cpusec)
  if (ifprint) then
     print*,''
     print*,'CPU TIME',cpusec
     print*,''
  end if

  ! stop

end subroutine gridset










