subroutine disk(converge,iter,i1)

  ! monte carlo radiative transfer for disk surrounding star.
  ! uses cartesian coordinates, arbitrary density distribution,
  ! arbitrary disk structure.
  ! This program works for extremely (geometrically) thin disks of large
  ! radial extent plus tenuous envelope.  In order to do the
  ! radiative transfer in the thin disk, the opacity is calculated
  ! at each step, and the photon path integration is variable.
  ! In the disk, the step size is small for directions perpendicular
  ! to the disk, say zmindisk/10.  Outside the
  ! disk, z gt. zmaxdisk, the step size is much larger, zmax/100
  ! or zmax/200, something on that order.
  ! Note that zmindisk is probably about .1 stellar radius.
  ! zmax is about 10000 stellar radii, that is, about 100 AU.
  ! The envelope extend is even larger for protostars--10^4 AU,
  ! so the step size increases with distance from the source.

  ! calls these subroutines:
  ! stokes
  ! initp
  ! opacdisk
  ! opacinfall
  ! dels
  !
  ! history:
  ! 99/04/12 (mjw): add cputime for photon loop update output
  ! (this breaks the overall CPU clock in newtts)
  ! 99/11/20 (baw):  big changes...
  !
  ! *************************************************************

  use tts_mod
  use grid_mod
  use stokes_mod
  use filt_mod
  use taunum_mod
  use log_mod
  use opacin_mod
  use out_mod
  use peeloff_mod
  use dust_mod
  use random
  use spot_mod
  use constants
  use output_mod
  use ttsre_mpi_mod

  implicit none

  real(8) :: sini(nmu),xmax
  real(8) :: thetb
  real(8) :: tsum,coslp,sinlp
  real(8) :: cosnorm,sipnew
  real(8) :: rmincgs,tcount,tau_therm,rn,lstar

  integer :: ii(3),ns,icount,ithet,iphot,ir
  integer :: it,ip,npacc,ips,idust2,id,npaccstar
  integer :: ispotph,iter,nsum,i1

  logical :: is_fuv

  character(len=4) :: suffix
  character(len=4) :: filename

  integer :: peelid

  ! cpu time variables
  real(8) :: cpusec,rjunk
  character(len=70) :: cmsgnm
  character(len=3) :: converge

  peeloff_time=0.d0
  peeloff1_time=0.d0
  peeloff2_time=0.d0
  peeloff3_time=0.d0
  peeloff4_time=0.d0
  peeliter = 0
  peeliter1 = 0
  peeliter2 = 0
  peeliter3 = 0
  peeliter4 = 0

  if (ifprint) write(6,*) 'check, pi',pi
  icount=0

  ! set up density grid
  if (ifprint) print*,'iteration #',iter,'calling gridset'
  call gridset(iter,i1)

!  call tempmap()
!  call linemap()
!  call depmap()
!  close(12)
!  call MPI_FINALIZE(irc)
!  stop

  write(suffix,'("_",I3.3)') iter

!  call output_grid('diffus' //suffix,diffus)
!  call output_grid('diffdir' //suffix,diffdir)

  if (ifprint) print*,''
  if (converge.eq.'yes'.and.ifprint) print*, 'final iteration'
  if (ifprint) print*,'iter',iter

  ! set up filter functions if making peeled images
  if (ipeel.eq.1.and.output_imfilt) call filt

  ! The flux is summed to image(ix,iy,it)
  ! polarization to imagei(ixp,iyp,it),imageq...,imageu...

  ! inclination arrays.
  do ithet=1,nmu
     sini(ithet)=sqrt(1.d0-u(ithet)**2)
  end do

  ! if the inner hole is large, star is effectively point source
  ! and theta,phi grid cells are really close together at rmin,
  ! so emit photons into radial direction so grid can handle it.
  if (rddust.gt.100.d0) then
     if (ifprint) print*,'emitting from star as point source'
     ips=1
  else
     ips=0
  end if

!  ips=1

  ! call this after gridset since some variables are set there.
  if (ispot.eq.1.and.spotflag.eq.1) then
     if (ifprint) print*,'calling spotset'
     call spotset()
  end if
  ! only call spotset once, it will set spotflag=0 after first call

  nscat=0.d0
  tot=0.d0

  ! origin of photons, set to 0 for now...
  i_orig=0

  ! x,y arrays for peeled-off image
  fractxh=2.d0*rmaxi/dble(nxhst)
  xmax=rmaxi

  rmincgs=rstar*rsol

  ! np=npsav
  if (ifprint) print*,'np',np
  flux=0.d0
  sflux=0.d0
  aflux=0.d0
  abflux=0.d0
  scount=0.d0
  dscount=0.d0

  npout=int(isrf_frac*dble(np))
  rn=isrf_frac*dble(np)
  if (ifprint) print*,'rn',rn

  if ((rn-npout).gt.0.5d0) npout=npout+1
! number of stellar photons
  ns=int(dble(np)*(1.d0-accfrac-isrf_frac))
  rn=real(np)*(1.d0-accfrac-isrf_frac)
  if (ifprint) print*,'rn',rn

  if ((rn-ns).gt.0.5d0) ns=ns+1
  npacc=int(dble(np)*accfrac)
  rn=int(dble(np)*accfrac)
  if (ifprint) print*,'rn',rn

  if ((rn-npacc).gt.0.5d0) npacc=npacc+1
!  number of hotspot photons
  npaccstar=int(accfrac2*dble(np))
  rn=int(dble(ns)*accfrac2)
  if (ifprint) print*,'rn',rn

  if ((rn-npaccstar).gt.0.5d0) npaccstar=npaccstar+1
  nsum=ns+npacc+npout

  if (ifprint) print*,'ns (star),npaccstar, npacc (disk),npout (outside),sum,should=np'
  if (ifprint) print*, ns,npaccstar,npacc,npout,nsum,np

  if (nsum.ne.np) then
     if (ifprint) print*,'resetting ns'
     ns=ns+np-nsum
     if (ifprint) print*,'ns (star),npacc (disk),npout (outside),sum,should=np'
     nsum=ns+npacc+npout
     if (ifprint) print*, ns,npacc,npout,nsum,np
  end if
  if (ifprint) print*,'np_disk,np_shock',npacc,npaccstar
  if (ifprint) print*,'np_spot = np_shock = npaccstar'
  
  lstar=4.d0*pi*sigt*tstar**4*rmincgs**2*(1.d0-fspot)
  if (ifprint) print*,'ltot/np*(ns-npaccstar),lstar',ltot/lsun/dble(np)*dble(ns-npaccstar),lstar/lsun

  ! convert Tstar to Jon's units
  if (iter.eq.iterstart) then
     tstar=tstar/11605.d0
     tstarave=tstarave/11605.d0
     tshock=tshock/11605.d0
  end if

  ! energy of fuv photon package
  ifuv_ori=sigt*(Tstar*11605.d0)**4*4.d0*pi*(rstar*rsol)**2/dble(ns)
  if (ifprint) print*,'i_fuv_ori',ifuv_ori/lsol

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  ! ********  do loop over stellar photons  *******
  if (ifprint) print*,'ns',ns
  if (ifprint) print*,'npout',npout

  np_fuv=0
  nps_fuv=0
  npd_fuv=0
  nphs_fuv=0

  do iphot=1,np

     on_wall = .false.

     if(mod(iphot,iwrite).eq.0) then
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        ! cpu time
        call cpu_time(cpusec)
        ! use explicit format to keep g77 happy
        write(cmsgnm,'(i12,a20,f11.2,a)') iphot,' photons completed. ',cpusec,' = CPU time (sec)'
        call WRIMSG('MAIN',cmsgnm)
        if (myid.eq.0) then
           print*,'photon completed',iphot,'CPU time',cpusec
        end if
        if (myid.eq.0) then
           write(filename,'(I4.4)') iphot/iwrite
           open(unit=15,file='IWRITE_'//filename//'.dat',status='unknown')
           write(15,*) 'photon completed',iphot
           close(15)
        end if
     end if

     ! if (ran().lt.accfrac) then
     if (iphot.lt.(npacc)) then ! DISK ACCRETION PHOTONS

        if (iphot.eq.1.and.ifprint) print*,'doing disk accretion photons now'

        call initpacc(ii,ns,iter)
        call opacset(nub)

        if (is_fuv(nub)) npd_fuv=npd_fuv+1

        wave=1.2398d0/nub

        ir=ii(1)
        it=ii(2)
        ip=ii(3)

        do id=1,nopac
           kapd(id)=kappa(id)*rmincgs*kappav(id)
           ! rlam(id)=albedo(id)
        end do

        iflag=0

        ! origin of photon, 2 for disk
        i_orig=2
        ! else
        ! i_orig=3
        ! end if
        ! i_orig=2

!        print*,ipeel
        if (ipeel.eq.1) then

           if (iflag.eq.1) then
              print*,'iflag=1'
              stop
           end if

           do peelid=1,npeel
              ! 20090730, idust is not defined, should be idust2, but should be okay, recalculated in peeloff
              sipnew = sip
              call peeloff(on_wall,xp,yp,zp,sipnew,sqp,sup,svp &
                   & ,cost,sint,cosp,sinp,phi &
                   & ,hit,htot,rsq,rtot &
                   & ,tsum,ii(1),ii(2),ii(3),idust &
                   & ,idust2,iflag,iphot &
                   & ,peelid)
           end do
        end if

        ! ******  check peeling off of accretion photon!!!!!!!  071902
        ! set normalization!!!!!

     else if (iphot.le.(npacc+ns)) then ! STELLAR AND HOT SPOT PHOTON

        if (iphot.lt.(npacc+npaccstar)) then ! HOT SPOT PHOTON
           if (iphot.eq.npacc+1.and.ifprint) print*,'doing hot spot photons now'

           ! 20080828 BAW selecting region inside hotspot
           ! iplanck=1
           ispotph=1
           if (ips.eq.0) then
              call initpspot(ispotph)
           else
              call initp_ps(ispotph)
           end if

           ! emit half the photons as x-rays
           if (ran().lt.0.5d0) then
              nub=1.2398d0/0.06d0+ran()*(1.2398d0/0.015d0-1.2398d0/0.06d0)
           end if

           if (is_fuv(nub)) nphs_fuv=nphs_fuv+1

        else ! STELLAR PHOTONS

           if (iphot.eq.npacc+npaccstar+1.and.ifprint) print*,'now doing stellar photons'

           ! 20080828 BAW, selecting region outside hotspot
           ! iplanck=0
           ispotph=0
           if (ips.eq.0) then
              call initpspot(ispotph)
           else
              call initp_ps(ispotph)
           end if
        end if

        if (is_fuv(nub)) nps_fuv=nps_fuv+1

        ! origin of photon i_orig = 1 for star
        i_orig=1

        call opacset(nub)
        wave=1.2398d0/nub

        do id=1,nopac
           kapd(id)=kappa(id)*rmincgs*kappav(id)
        end do

        coslp=cos(lp)
        sinlp=sin(lp)

        zp=cosb*rmin
        rp=sinb*rmin
        xp=rp*coslp
        yp=rp*sinlp

        rsq=zp**2+rp**2
        rtot=sqrt(rsq)

        ux=sint*cosp
        uy=sint*sinp
        uz=cost

        iflag=0
        ii(1)=1 ! index of radial grid (stellar surface is 1)
        on_wall = .true.

        thetb=acos(cosb)
        call locate(thetarr,ntg,thetb,it)
        ii(2)=it

        call locate(phiarr,npg,lp,ip)
        ii(3)=ip

        ! first, peel off direct flux
        ! weight photon intensity by angle between normal and
        ! photon direction, since emitted from a surface.
        if (ipeel.eq.1.and..not.partial_peeloff) then

           if (iflag.eq.1) then
              print*,'iflag=1'
              stop
           end if

           do peelid=1,npeel

              if (ips.eq.0) then

                 cosnorm=cosb*coste_arr(peelid)+(sinb*sinte_arr(peelid)*(cospe_arr(peelid)*coslp+sinpe_arr(peelid)*sinlp))

                 if (limb.eq.0) then
                    ! intensity constant, energy/Sr proportional to mu
                    sipnew=4.d0*sip*cosnorm !normalization = 2
                 else
                    ! intensity goes as (1+mu), energy/Sr has another factor of mu.
                    sipnew=12.d0/5.d0*sip*(cosnorm+cosnorm*cosnorm)
                 end if

                 if (cosnorm.gt.0.d0) then
                    call peeloff(on_wall,xp,yp,zp,sipnew,sqp,sup,svp &
                         & ,cost,sint,cosp,sinp,phi &
                         & ,hit,htot,rsq,rtot &
                         & ,tsum,ii(1),ii(2),ii(3),idust &
                         & ,idust2,iflag &
                         & ,iphot,peelid)
                 end if

              else

                 sipnew=sip*2.d0
                 call peeloff(on_wall,xp,yp,zp,sipnew,sqp,sup,svp &
                      & ,cost,sint,cosp,sinp,phi &
                      & ,hit,htot,rsq,rtot &
                      & ,tsum,ii(1),ii(2),ii(3),idust &
                      & ,idust2,iflag,iphot &
                      & ,peelid)

              end if

           end do
        end if

     else ! EXTERNAL ILLUMINATION

        if (iphot.eq.npacc+ns+1.and.ifprint) print*,'now doing external photons'

        call initpout()
        on_wall=.true.
        i_orig=5 !!!! changed from 4 to 5 since we use 4 for outflow
        call opacset(nub)
        wave=1.2398d0/nub

        do id=1,nopac
           kapd(id)=kappa(id)*rmincgs*kappav(id)
        end do

        coslp=cos(lp)
        sinlp=sin(lp)
        zp=cosb*rmax*0.999999
        rp=sinb*rmax*0.999999
        xp=rp*coslp
        yp=rp*sinlp

        rsq=zp**2+rp**2
        rtot=sqrt(rsq)

        ux=sint*cosp
        uy=sint*sinp
        uz=cost

        iflag=0
        ii(1)=nrg-1          !index of radial grid at rmax

        thetb=acos(cosb)
        call locate(thetarr,ntg,thetb,it)
        ii(2)=it

        call locate(phiarr,npg,lp,ip)
        ii(3)=ip

        ! first, peel off direct flux
        ! I think we weight photon simply by 1/4pi, or just one?  no, 2 because
        ! we only count half the flux.
        ! however, we don't want to include these photons unless they
        ! interact. anlogous to observing, it's like subtracting off the sky
        ! background. so don't call peel here.
        ! if (ipeel.eq.1) then
        ! cosnorm=cosb*coste+(sinb*sinte*(cospe*coslp+
        ! 1                 sinpe*sinlp))
        ! sipnew=sip*2.d0
        ! if (cosnorm.lt.0.d0) then
        ! call peeloff(xp,yp,zp,sipnew,sqp,sup,svp
        ! 1                    ,cost,sint,cosp,sinp,phi
        ! 1                    ,pi,r2p,hit,htot,rsq,rtot
        ! 1                    ,tsum,ii,idust2,iflag,iphot,eps)
        ! end if
        ! if (cosnorm.gt.0.and.iveeg.eq.1) then
        ! call vger(xp,yp,zp,sipnew,sqp,sup,svp
        ! 1                    ,cost,sint,cosp,sinp,phi
        ! 1                    ,pi,r2p,hit,htot,rsq,rtot
        ! 1                    ,tsum,ii,idust,iflag,iphot,eps)
        ! end if
        ! end if

     end if

     if (is_fuv(nub)) then
        ifuv=ifuv_ori 
     else 
        ifuv=0.d0
     end if

!     if (i_orig.eq.2) ifuv=0.d0

!     print*,'propagate',iphot
     call propagate(iphot,sini,ii,ns,idust2,iter,converge)
     
5 end do
  
  ! ***** end of loop over stellar photons  ******

  ! np=ns

  np_fuv=nps_fuv+npd_fuv+nphs_fuv

  if (ifprint) print*,'photons killed: ', myid,sum(killed)
  if (ifprint) write(12,*) 'photons killed: ', sum(killed)

  if (ifprint) print*, 'fuv photon from the star', nps_fuv
  if (ifprint) print*, 'fuv photon from the disk', npd_fuv
  if (ifprint) print*, 'fuv photon from the hotspot', nphs_fuv
  if (ifprint) print*, 'total fuv photon, total stellar photon, total photon', np_fuv, ns, np
  if (ifprint) print*, 'total fuv luminosity (lsun)', dble(np_fuv)*ifuv_ori/lsol 

  print*,'process',myid,':  done with photons'

  if (ifprint) print*,'total peeloff time',peeloff_time
  if (ifprint) print*,'total peeloff1 time',peeloff1_time
  if (ifprint) print*,'total peeloff2 time',peeloff2_time
  if (ifprint) print*,'total peeloff3 time',peeloff3_time
  if (ifprint) print*,'total peeloff4 time',peeloff4_time
  if (ifprint) print*,'total peeloff iter',peeliter, peeliter1, peeliter2,peeliter3,peeliter4

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)


! MPI reduction***********************************************
! ************************************************************

  ! some scalars
  call MPI_ALLREDUCE(MPI_IN_PLACE,icount,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,np,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,dscount,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,scount,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  write(12,*) 'REDUCE: Number of photons: np',np
  if (ifprint) print*,'REDUCE: Number of photons: np',np

! quantities which may need to be summed, but are not explicitly in the routine
! KILLED (reset before each call to DISK)
  icountbuf = nrg*ntg*npg
  call MPI_ALLREDUCE(MPI_IN_PLACE,killed,icountbuf,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  write(12,*) 'REDUCE: photons killed: ', sum(killed)
  if (ifprint) print*,'REDUCE: photons killed: ', sum(killed)

! quantities which may need to be summed, but are not explicitly in the routine
! NABS,DBIG,DSG (reset before each call to DISK)  
! JMEAN not needed as it is calculated from DBIG and DSG
  icountbuf = nrg*ntg*npg
  call MPI_ALLREDUCE(MPI_IN_PLACE,nabs,icountbuf,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,dtauabs,icountbuf,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,dbig,icountbuf,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,dsg,icountbuf,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,jmean,icountbuf,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,fuvmean,icountbuf,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr) ! note fuvmean includes energy per photon
  call MPI_ALLREDUCE(MPI_IN_PLACE,tdust,icountbuf,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  tdust=tdust/dble(numprocs)

! SI,SQ,SU,SV,SI2,SQ2,SU2,SV2,NUMS (these are reset before each call to DISK)
  icountbuf = nfreq*nmu*nph*nap*no
  call MPI_ALLREDUCE(MPI_IN_PLACE,si,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sq,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,su,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sv,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,si2,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sq2,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,su2,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sv2,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,nums,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
  write(12,*) 'REDUCE: SED photons used: ', sum(nums(:,:,:,:,1))
  if (ifprint) print*,'REDUCE: SED photons used: ', sum(nums(:,:,:,:,1))

 ! do reduction for arrays available only when temperature is converged
  if (converge.eq.'yes') then
     write(12,*) 'Reducing image arrays...'

     icountbuf = nfreq*npeel*nap*no
     call MPI_ALLREDUCE(MPI_IN_PLACE,ti,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE,tq,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE,tu,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE,tv,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE,ti2,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE,tq2,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE,tu2,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE,tv2,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE,numt,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
     write(12,*) 'REDUCE: peeled SED photons used: ', sum(numt(:,:,:,1))
     if (ifprint) print*, 'REDUCE: peeled SED photons used: ', sum(numt(:,:,:,1))

! broadband peeloff IMAGES
     if(output_imfilt) then
        icountbuf = npeel*nxhst*nxhst*nbnd
        call MPI_ALLREDUCE(MPI_IN_PLACE,image_b_i,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,image_b_q,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,image_b_u,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,image_b_v,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,image_b_d,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,image_b_s,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,image_b_e,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,image_b_o,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
     endif
     
     if(output_imcube) then
        ! monochromatic peeloff IMAGES
        icountbuf = npeel*nxhst*nxhst*nfreq
        call MPI_ALLREDUCE(MPI_IN_PLACE,image_m_i,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
        icountbuf = npeel*nxhst*nxhst*nfreq
        call MPI_ALLREDUCE(MPI_IN_PLACE,image_m_q,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,image_m_u,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,image_m_v,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
     endif
  end if

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  write(cmsgnm,'(a,i2)') 'Done with reduction in process ',myid
  call wrimsg('DISK',cmsgnm)
  if (ifprint) print*,'Done with reduction in process'

!*******************************************************
!******************************************************

  write(suffix,'("_",I3.3)') iter

  fuvmean=fuvmean/1.6d-3/dble(numprocs)
  if (iffuv) call output_grid('fuv'//suffix,fuvmean)

  if (ifprint) print*, 'calling tfinal'
  call tfinal(ns*numprocs,iter,converge)
!  call tfinal(ns,iter,converge)

  call int_mean(ns+npout+npacc)

  write(12,*) 'fraction of phots scattered outside image'
  write(12,*) real(icount)/real(np)
  write(12,*) np,scount,dscount
  if (ifprint) print*,np,scount,dscount
  tcount=(dble(np)-scount-dscount)/dble(np)
  tau_therm=-log(1.d0-tcount)

  call diskdat_section('newdisktherm')
  call diskdat_write('tcount',tcount,'')
  call diskdat_write('tautherm',tau_therm,'')

  if (sum(fmass,mask=is_sg(1:8)).gt.1.d-6) then
     rjunk=real(imave)/real(imct)
  else
     rjunk=0.
  end if

  write(12,*) 'imave,imct',imave,imct
  write(12,*) 'average imean index for PAH emissivity',rjunk
  call diskdat_write('imave/imct',rjunk,'average imean index for PAH emissivity')

!  do ir=1,nrg-1
!     do it=1,ntg-1
!        do ip=1,npg-1
!           if(killed(ir,it,ip) > 1 .and. sum(densarr(ir,it,ip,:), mask=is_thermal) > 0.) then
!              diffus(ir,it,ip) = .true.
!           end if
!        end do
!     end do
!  end do


!  call output_grid('denstotal'//suffix,sum(densarr(:,:,:,1:4), dim=4))
!  call output_grid('killed' //suffix,killed)

!  call output_grid('imean'//suffix,imean)

  return
end subroutine disk

! *********************************************************
