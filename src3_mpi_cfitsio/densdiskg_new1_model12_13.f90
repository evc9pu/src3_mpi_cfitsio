subroutine densdiskg_new(rmincgs,iter,lstar,lacc,lacc2)

  use tts_mod
  use grid_mod
  use dust_mod
  use opacin_mod
  use constants
  use out_mod
  use log_mod
  use ttsre_mpi_mod

  implicit none

  integer :: ir,it,ip,id,iter,dcount
  integer :: ial

  real(8) :: mdot,mdotcgs,rmincgs
  real(8) :: zwall,Co_md,Co_Lw,Co_Ld
  real(8) :: rad,dr,dr3,alpha,rho
  real(8) :: rho1,rinfall,mdotcorecgs,fd,rad1
  real(8) :: thet,cost,dcost,zp,gexp,zmr,zmr1,dens,temp
  real(8) :: frac1,frac2,rhothr1,rhothr2,dp,vol,diskscale
  real(8) :: lacc,lstar,lacc2,Tmin,thetemp,rp
  real(8) :: rho2,v_Kc,alphal,alphah,alpha_temp
  real(8) :: massc_in,rmincgs_in,mdotcgs_in,rmaxd_in,rinfall_in,fd_in
  real(8) :: mdotcorecgs_in,dt,diffsigma,sumsigma,tnow

  real(8) :: sigmaarr(nrg-1),zmrarr(nrg-1),sigmadotarr(nrg-1)
  real(8) :: sigmaarr1(nrg-1),sigmadotarr_new(nrg-1)

  integer :: ir_zdmax,iiter

  logical :: ifout

  Tmin=1400.d0

  rhothr1=3.34d-14

  if (ifprint) print*,'set up density for disk'

  massdisk=0.d0

  if (massd.gt.0.d0) then

     if (fmass(5)+fmass(1).eq.0.d0.or.fmass(6)+fmass(2).eq.0.d0) then
        print*,'both two disks need to be set!'
        stop
     end if

     if (fmass(5).eq.0.d0) then
        frac1=0.d0
     else
        frac1=fmass(5)/(fmass(1)+fmass(5))
     end if

     if (fmass(6).eq.0.d0) then
        frac2=0.d0
     else
        frac2=fmass(6)/(fmass(2)+fmass(6))
     end if

     if (ifprint) print*,'frac1,frac2',frac1,frac2

     mdot=mdotdisk
     if (ifprint) print*,'mdot',mdot
     mdotcgs=mdot*msol/3600.d0/24.d0/365.25d0
     rinfall=rchole*(sin(thet1*deg2rad))**2
     if (ifprint) print*,'rinfall in AU',rinfall/autors
     fd=massd/massc
     if (ifprint) print*,'fd',fd,'fw',fw,'fwesc',fwesc
     mdotcorecgs=(1.d0+fd+fw*fwesc)*mdotcgs/cos(thet1*deg2rad)

     if (iter.eq.iterstart) then

!!$     alphal=.5d0 !!! adjust this
!!$     alphah=100.d0  !!! adjust this

        alphal=alphamin
        alphah=alphamax
        
        sigmadotarr=0.d0
        
        iiter=0
        
10      massc_in=massc
        rmincgs_in=rmincgs
        mdotcgs_in=mdotcgs
        rmaxd_in=rmaxd
        rinfall_in=rinfall
        fd_in=fd
        mdotcorecgs_in=mdotcorecgs
        
        call alphadisk(massc_in,rmincgs_in,mdotcgs_in,rmaxd_in,rinfall_in,fd_in &
             & ,mdotcorecgs_in,sigmadotarr,alphah,alphal,alpha,sigmaarr)
        
        massc_in=massc*1.1d0
        rmincgs_in=rmincgs
        mdotcgs_in=mdotcgs*1.1d0**0.5d0
        rmaxd_in=rmaxd*1.1d0**(2.d0/3.d0)
        rinfall_in=rmaxd_in*(sin(thet1*deg2rad))**2
        fd_in=fd
        mdotcorecgs_in=(1.d0+fd+fw*fwesc)*mdotcgs_in/cos(thet1*deg2rad)
        
        alpha_temp=alpha

        call alphadisk(massc_in,rmincgs_in,mdotcgs_in,rmaxd_in,rinfall_in,fd_in &
             & ,mdotcorecgs_in,sigmadotarr,alphah,alphal,alpha,sigmaarr1)
        
        tnow=1.29d5*sqrt(massc*(1.d0+fd)/mcore)
        dt=(sqrt(1.1d0)-1.d0)*tnow*365.25d0*3600.d0*24.d0
        sigmadotarr_new=(sigmaarr1-sigmaarr)/dt
        
        diffsigma=0.d0
        sumsigma=0.d0
        do ir=1,nrg-1
           diffsigma=diffsigma+(sigmadotarr_new(ir)-sigmadotarr(ir))**2
           sumsigma=sumsigma+(sigmadotarr_new(ir))**2
        end do
        diffsigma=sqrt(diffsigma/sumsigma)
        
        iiter=iiter+1
        
        if (ifprint) print*,'diffsigma',iiter,diffsigma
        if (ifprint) print*,''
        
        sigmadotarr=sigmadotarr_new
        
        if (diffsigma.gt.1.d-3.and.iiter.le.10) goto 10
        
        if (ifprint) print*,'iteration for sigmadot',iiter

        alphafinal=alpha_temp
        sigmadotarrf=sigmadotarr

     end if
        
     massc_in=massc
     rmincgs_in=rmincgs
     mdotcgs_in=mdotcgs
     rmaxd_in=rmaxd
     rinfall_in=rinfall
     fd_in=fd
     mdotcorecgs_in=mdotcorecgs
     ifout=.true.

     sigmaarr=0.d0
     zmrarr=0.d0

     call diskprofile(massc_in,rmincgs_in,mdotcgs_in,rmaxd_in,rinfall_in, &
          & fd_in,mdotcorecgs_in,sigmadotarrf,alphafinal,ifout,Co_md,Co_Lw, &
          & Co_Ld,sigmaarr,zmrarr)

     print*,'myid,alpha,Co_md/fd,Co_Lw/Co_Ld',myid,alphafinal,Co_md/fd,Co_Lw/Co_Ld

     if (iter.eq.iterstart) then
        call diskdat_write('alpha',alphafinal, 'for final disk solution')
     end if

     rho1=0.d0
     do ir=1,nrg-1
        rad=ravearr(ir)
        zmr=zmrarr(ir)
        rho=0.d0
        if (zmr.gt.0.d0) then
           rho=sigmaarr(ir)/sqrt(2.d0*pi)/(zmr*rmincgs)
           rho1=rho
           zmr1=zmr
        end if
        do it=1,ntg-1
           thet=0.5d0*(thetarr(it)+thetarr(it+1))
           cost=cos(thet)
!           zp=rad*cost
!           gexp=(zp/zmr)**2
           gexp=(rad*(pi/2.d0-thet)/zmr)**2
           dens=rho*exp(-0.5d0*gexp)
           do ip=1,npg-1
              temp=tdave(ir,it,ip)
              densarr(ir,it,ip,:)=0.d0
              if (temp.gt.Tmin/11605.d0) then
                 densarr(ir,it,ip,9)=dens
              else if (dens.lt.rhothr1) then
                 densarr(ir,it,ip,2)=dens*(1.d0-frac2)
                 densarr(ir,it,ip,6)=dens*frac2
              else
                 densarr(ir,it,ip,1)=dens*(1.d0-frac1)
                 densarr(ir,it,ip,5)=dens*frac1
              end if
           end do
        end do
     end do

     v_Kc=sqrt(gn*massc*msol/rstar/rsol)
     zdisk=0.d0
     rhod=0.d0
     radmax=0.d0
     do ir=nrg-1,1,-1
        rad=ravearr(ir)
        if (rad.gt.rchole*(1.d0-windmu0**2).and.rad.lt.rmaxd) then
           zwall=windmu0-(rchole/rad)*windmu0*(1.d0-windmu0**2)
           if (zwall.gt.1.d0.or.zwall.lt.-1.d0) then
              print*,'zwall wrong for stream line shape'
              print*,'zwall,windmu0,rchole,rad',zwall,windmu0,rchole/autors,rad/autors,autors
              stop
           end if
           zwall=asin(zwall)*rad !spherical surface
           gexp=(zwall/hdisk(ir))**2
           rho1=midrho(ir)*exp(-0.5d0*gexp)
           rp=rad*cos(zwall/rad)
           rho2=mdotdisk*msol/3600.d0/24.d0/365.25d0*fw*fwesc
           rho2=rho2/4.d0/pi/(rstar*rsol)**2/v_Kc/log(rp)
           rho2=rho2*rp**(-1.5d0)/log(1.01d0)
!           print*,rad/autors,rho1,rho2
           if (rho2.gt.rho1) then
              radmax=ravearr(ir-1)
              ir_zdmax=ir-1
           end if
        end if
     end do

     if (ifprint) print*,'radmax',radmax,radmax/autors

     zwall=windmu0-(rchole/rmaxd)*windmu0*(1.d0-windmu0**2)
     if (zwall.gt.1.d0.or.zwall.lt.-1.d0) then
        print*,'zwall wrong for stream line shape'
        print*,'zwall,windmu0,rchole,rad',zwall,windmu0,rchole/autors,rad/autors,autors
        stop
     end if
     zwall=asin(zwall)*rmaxd !spherical surface
     gexp=(zwall/zmr1)**2
     rhothr2=rho1*exp(-0.5d0*gexp)

     if (radmax.eq.0.d0) then

!!$        print*,'radmax wrong'
!!$        stop

        if (ifprint) print*,'Warning!!! radmax set to be largest'
        radmax=rmaxd
        do ir=1,nrg-1
           if (ravearr(ir).lt.rmaxd) ir_zdmax=ir
        end do
        zdmax=radmax*windmu0-rchole*windmu0*(1.d0-windmu0**2)
        x_max0=sqrt(radmax**2-zdmax**2)
        do ir=1,nrg-1
           rad=ravearr(ir)
           if (rad.lt.radmax) then
              do it=1,(ntg-1)/2
                 thet=thetarr(it)
                 gexp=((pi/2.d0-thet)*rad/hdisk(ir))**2
                 rho1=midrho(ir)*exp(-0.5d0*gexp)
                 if (rho1.gt.rhothr2) exit
              end do
              zdisk(ir)=rad*cos(thet)
              cyndisk(ir)=rad*sin(thet)
              rhod(ir)=rho1
              do it=1,(ntg-1)/2
                 thet=thetarr(it)
                 gexp=((pi/2.d0-thet)*rad/hdisk(ir))**2
                 rho1=midrho(ir)*exp(-0.5d0*gexp)
                 rp=rad*sin(thet)
                 rho2=mdotdisk*msol/3600.d0/24.d0/365.25d0*fw*fwesc
                 rho2=rho2/4.d0/pi/(rstar*rsol)**2/v_Kc/log(x_max0)
                 rho2=rho2*rp**(-1.5d0)/log(1.01d0)
                 if (rho1.gt.rho2) exit
              end do
!              print*,rad/autors,rhod(ir),rho1
              if (rad*cos(thet).lt.zdisk(ir)) then
                 zdisk(ir)=rad*cos(thet)
                 cyndisk(ir)=rad*sin(thet)
                 rhod(ir)=rho1
              end if
           else
              rhod(ir)=rhod(ir-1)
              cyndisk(ir)=cyndisk(ir-1)
           end if
        end do

     else

        zdmax=radmax*windmu0-rchole*windmu0*(1.d0-windmu0**2)
        x_max0=sqrt(radmax**2-zdmax**2)
        do ir=1,nrg-1
           rad=ravearr(ir)
           if (rad.lt.radmax) then
              do it=1,(ntg-1)/2
                 thet=thetarr(it)
                 gexp=((pi/2.d0-thet)*rad/hdisk(ir))**2
                 rho1=midrho(ir)*exp(-0.5d0*gexp)
                 rp=rad*sin(thet)
                 rho2=mdotdisk*msol/3600.d0/24.d0/365.25d0*fw*fwesc
                 rho2=rho2/4.d0/pi/(rstar*rsol)**2/v_Kc/log(x_max0)
                 rho2=rho2*rp**(-1.5d0)/log(1.01d0)
                 if (rho1.gt.rho2) exit
              end do
              zdisk(ir)=rad*cos(thet)
              cyndisk(ir)=rad*sin(thet)
              rhod(ir)=rho1
           else
              rhod(ir)=rhod(ir-1)
              cyndisk(ir)=cyndisk(ir-1)
           end if
        end do

     end if

     if (ifprint) print*,'zdmax',zdmax,zdmax/autors
     if (ifprint) print*,'x_max0',ir_zdmax,x_max0,x_max0/autors

     if (iter.eq.iterstart.and.ifprint) then
        open(unit=98,file='zdisk.dat',status='unknown')
        do ir=1,nrg-1
           rad=ravearr(ir)
           write(98,*) rad/autors,zdisk(ir)/autors,cyndisk(ir)/autors
        end do
        close(98)
     end if

     do ir=1,nrg-1
        dr3=rarr(ir+1)**3.d0-rarr(ir)**3.d0
        do it=1,ntg-1
           dcost=cos(thetarr(it+1))-cos(thetarr(it))
           do ip=1,npg-1
              dp=phiarr(ip+1)-phiarr(ip)
              vol=-dr3*dcost*dp/3.d0*rmincgs**3.d0
              do id=1,ndg+2
                 if (is_disk(id)) then
                    dens=densarr(ir,it,ip,id)
                    if (dens.le.rhod(ir)) dens=0.d0
                    densarr(ir,it,ip,id)=dens
                    massdisk=massdisk+dens*vol
                 end if
              end do
           end do
        end do
     end do

     massdisk=massdisk/msol
     diskscale=massd/massdisk
     if (ifprint) print*,'massdisk before scaling', massdisk
     if (ifprint) print*,'diskscale',diskscale

     massdisk=0.d0
     do ir=1,nrg-1
        dr3=rarr(ir+1)**3.d0-rarr(ir)**3.d0
        do it=1,ntg-1
           dcost=cos(thetarr(it+1))-cos(thetarr(it))
           do ip=1,npg-1
              dp=phiarr(ip+1)-phiarr(ip)
              vol=-dr3*dcost*dp/3.d0*rmincgs**3.d0
              do id=1,ndg+2
                 if (is_disk(id)) then
                    densarr(ir,it,ip,id)=densarr(ir,it,ip,id)*diskscale
                    massdisk=massdisk+densarr(ir,it,ip,id)*vol
                 end if
              end do
           end do
        end do
     end do

     do ir=1,nrg
        do it=1,ntg
           do ip=1,npg
              dcount=0
              do id=1,ndg+2
                 if (.not.is_sg(id).and.densarr(ir,it,ip,id).gt.0.d0) then
                    dcount=dcount+1
                 end if
              end do
              if (dcount.gt.1) then
                 print*,'more than one non-zero density',ir,it,ip,densarr(ir,it,ip,:)
                 stop
              end if
           end do
        end do
     end do

     if (iter.eq.iterstart.and.idiskacc.eq.1) then

        if (ifprint) print*,'adjust accretion luminosity from',lacc,'to',Co_Ld
        lacc=Co_Ld
        ltot=lacc+lacc2+lstar+l_isrf
        accfrac=lacc/ltot
        if (ifprint) print*,'lacc,lacc2,lstar,l_isrf',lacc/lsun,lacc2/lsun,lstar/lsun,l_isrf/lsun

        accfrac=lacc/ltot
        accfrac2=lacc2/ltot
        isrf_frac=l_isrf/ltot
        if (ifprint) print*, 'new fraction of acc luminosity from disk ',accfrac
        if (ifprint) print*, 'new fraction of lum. on stellar hotspot ',accfrac2
        if (ifprint) print*, 'new total system luminosity (Lsun)',ltot/lsun

        call diskdat_write('isrf_frac',isrf_frac, 'new fractional ISRF lum')
        call diskdat_write('diskacc',.true.,'whether disk luminosity is included')
        call diskdat_write('accfrac',accfrac,'new fraction of lum. from disk')
        call diskdat_write('accfrac2',accfrac2,'new fraction of lum. on stellar hotspot')
        call diskdat_write('ltot',ltot/lsun,'new total system lum. (Lsun)')

     end if
     
  else

     if (ifprint) print*,'disk mass = 0'
     massdisk=0.d0

  end if

!  stop

  return

end subroutine densdiskg_new
