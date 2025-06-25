subroutine propagate(iphot,sini,dum,nphot,idust2,iter,converge)

  use tts_mod
  use grid_mod
  use stokes_mod
  use taunum_mod
  use opacin_mod
  use tauint_mod
  use peeloff_mod
  use dust_mod
  use random
  use constants
  use ttsre_mpi_mod

  implicit none

  real(8) :: xran,xpold,ypold,rho2,zpold,x,y
  real(8) :: tsum,pol2,sipnew
  real(8) :: sini(nmu)
  integer :: ii(3),dum(3),iphot,iscat,ia
  integer :: ir,it,ip,inub,iint,iintmax,nphot
  integer :: iabs,idust2,id,iter,i

  integer :: indx(3),i_inter, ii_prev(3)

  integer :: peelid

  logical :: is_fuv

  character(len=3) :: converge

  iintmax = 10000
  
  xpold=xp
  ypold=yp
  zpold=zp

  do i=1,3
     ii(i)=dum(i)
  end do

  ! iflag is set to 1 if photon scatters
  iflag=0
  exitflag=0
  aflag=0
  iscat=0
  iabs=0
  iint=0

  ! sample tau
  xran=ran()
  tau=-log(xran)

  ! integrate over distance until the optical depth equals tau.
  call tauint(iphot,tsum,ii(1),ii(2),ii(3))
  ir=ii(1)
  it=ii(2)
  ip=ii(3)

  if(exitflag.eq.1) go to 300
  if(aflag.eq.1) then
     write(6,*) 'shouldnt be here, aflag=1'
     return
  end if

  !  randomly sample dust type.0
!  print*,'idust2',idust2
  call select_dust(ir,it,ip,is_any,idust2)
!  print*,'idust2',idust2
!  stop

  tot=tot+1.d0

  ! photon scatters/absorbs/emits until exit exits disk

  do

     if(ran().le.albedo(findopac(ir,it,ip,idust2))) then ! SCATTERED
        iflag=1
        iscat=iscat+1
        iint=iint+1
        sipnew=sip
        if (ipeel.eq.1) then
           do peelid=1,npeel
              call peeloff(on_wall,xp,yp,zp,sipnew,sqp,sup,svp &
                   &,cost,sint,cosp,sinp,phi&
                   &,hit,htot,rsq,rtot&
                   &,tsum,ii(1),ii(2),ii(3),idust,idust2,iflag,iphot,peelid)
           end do
        end if

        pol2=sqp**2+sup**2+svp**2
        if (pol2 .gt. sip**2) then
           print*,'error, P^2, I^2 ', pol2,sip**2
           continue
        end if

        if (idust.lt.2) then
           call stokes(findopac(ir,it,ip,idust2))
        else
           print*,'oops, error, idust is wrong'
           stop
        end if

        pol2=sqp**2+sup**2+svp**2
        if (pol2 .gt. sip**2) then
           print*,'error, P^2, I^2 ', pol2,sip**2
           continue
        end if

!        ifuv=0.d0

     else ! ABSORBED + RE-EMITTED

        ! tabulate which emission process is chosen as a function of radius
        if (is_sg(idust2)) then
           dsg(ir,it,ip)=dsg(ir,it,ip)+1
        else
           dbig(ir,it,ip)=dbig(ir,it,ip)+1
        end if

        ! set photon origin to disk or envelope
        if (idust2.eq.1.or.idust2.eq.2.or.idust2.eq.5.or.idust2.eq.6.or.idust2.eq.9) then
           i_orig=2
        else if (idust2.eq.3.or.idust2.eq.7) then
           i_orig=3
        else if (idust2.eq.4.or.idust2.eq.8.or.idust2.eq.10) then
           i_orig=4 ! use 4 for outflow
        else
           print*,'i_orig wrong'
           stop
        end if

        if (albedo(findopac(ir,it,ip,idust2)).eq.1.d0) then
           print*,'error! albedo=1, yet photon absorbed'
           stop
        end if

        if(.not.diffus(ir,it,ip)) nabs(ir,it,ip)=nabs(ir,it,ip)+1

!        print*,'OK1',ir,it,ip
        call emit_common(ir,it,ip,idust2,iter,ii,nphot)
!        print*,'OK2'

        iflag=0

        call opacset(nub)
        wave=1.2398d0/nub

        if (is_fuv(nub)) then
           ifuv=ifuv_ori 
!           ifuv=0.d0
        else 
           ifuv=0.d0
        end if

        do id=1,nopac
           kapd(id)=kappa(id)*rstar*rsol*kappav(id)
        end do

        sipnew=sip
        if (ipeel.eq.1.and..not.partial_peeloff) then
           ! don't need to know idust because it's not a scattered photon
           do peelid=1,npeel
              call peeloff(on_wall,xp,yp,zp,sipnew,sqp,sup,svp&
                   &,cost,sint,cosp,sinp,phi,&
                   &hit,htot,rsq,rtot,&
                   &tsum,ii(1),ii(2),ii(3),idust,idust2,iflag,iphot,peelid)
           end do
        end if

        iint=iint+1

     end if

     ! sample tau
     tau=-log(ran())

     ! integrate distance until optical depth equals tau.
     xpold=xp
     ypold=yp
     zpold=zp

     ! not sure if we need to do this here, might already have this
     ! rsq=xp**2+yp**2+zp**2
     ! rtot=sqrt(rsq)
     ! rpold=rp

     ux=sint*cosp
     uy=sint*sinp
     uz=cost
     ii_prev = ii
     call tauint(iphot,tsum,ii(1),ii(2),ii(3))
     ir=ii(1)
     it=ii(2)
     ip=ii(3)

     if(exitflag.eq.1) exit

     if(any(ii.ne.ii_prev)) iint = 0.
     
     if (iffuv) then
        if (.not.is_fuv(nub)) then
           go to 400
        end if
     end if

     if (iint.gt.iintmax) then
        killed(ii(1),ii(2),ii(3)) = killed(ii(1),ii(2),ii(3)) + 1
        abflux=abflux+1
        ! print*,'killing photon (too many interactions)'
        go to 400
     end if

     !  randomly sample dust type.
     call select_dust(ir,it,ip,is_any,idust2)

  end do

300 continue

!  if (iter.eq.3) print*,iint

  ! Bin angle
  ! NOTE:  assuming axisymmetric, and z=-z.
  ! k=int(real(nmu-1)*abs(cost)+1.5d0)
  ! for comparison to jon
  ! k=min(int(real(nmu)*abs(cost)+1),nmu)
  ! 20080826 new binning in cost, not combining about midplane

  it=int(dble(nmu)*(1.d0-cost)/2.d0)+1
  if(it.lt.0.or.it.gt.nmu) then
     print*,'cost binning error, it, nmu',it,nmu
     print*, 'cost,sint',cost,sint
     print*, 'iphot',iphot
  end if

  if (phi.lt.0.d0) phi=phi+r2p
  ip=int(dble(nph)*phi/r2p)+1
  if(ip.lt.0.or.ip.gt.nph) then
     print*, 'phi binning error, m, nph ',ip,nph
     print*, 'phi ', phi
     print*, 'iphot ',iphot
  end if

  ! Old imaging.  (not even used anymore but don't ever want to have
  ! to figure these out again)
  ! if (cost.lt.0.d0) then
  !    x = -ypold*cosp+xpold*sinp
  !    y = -zpold*sini(k)-ypold*u(k)*sinp
  !     &        -xpold*u(k)*cosp
  ! else
  !    x = ypold*cosp-xpold*sinp
  !    y = zpold*sini(k)-ypold*u(k)*sinp
  !     &        -xpold*u(k)*cosp
  ! end if

  ! 20080826 imaging with no mirror of cost (not used anymore)
  x = ypold*cosp-xpold*sinp
  y = zpold*sini(it)-ypold*u(it)*sinp-xpold*u(it)*cosp

  ! first, sum fluxes
  rho2=(x**2+y**2)
  if (rho2.lt.aperture2(nap)*1.0001d0) then
     flux=flux+1.d0
     if (iflag.eq.1) then
        sflux=sflux+1
        nscat=nscat+iscat
     end if
  end if

  ! Find frequency bin, and make sure it is inside the limits
  inub=min(max(int(dble(nfreq)*log(nub/numin)/lnurat)+1,1),nfreq)
  if (inub.lt.1.or.inub.gt.nfreq) then
     print*,'inub out of limits!,inub,nub',inub,nub
  end if

  if (i_orig.lt.1.or.i_orig.gt.5) print*,'error, IOR wrong',i_orig

  ! Loop over apertures, and bin photon into SEDs

  if(iflag==1) then
     i_inter=2
  else
     if(i_orig==1) then
        i_inter=1
     else
        i_inter=3
     end if
  end if

  indx = (/1,1+i_orig,5+i_inter/)

  do ia=1,nap
     if(rho2 < aperture2(ia)) then

        if (.not.(iflag.eq.0.and.i_orig.eq.5)) then

           !    add photon to flux arrays
           si(inub,it,ip,ia,indx) = si(inub,it,ip,ia,indx) + real(sip)
           sq(inub,it,ip,ia,indx) = sq(inub,it,ip,ia,indx) + real(sqp)
           su(inub,it,ip,ia,indx) = su(inub,it,ip,ia,indx) + real(sup)
           sv(inub,it,ip,ia,indx) = sv(inub,it,ip,ia,indx) + real(svp)

           !    add photon to flux error arrays
           si2(inub,it,ip,ia,indx) = si2(inub,it,ip,ia,indx) + real(sip*sip)
           sq2(inub,it,ip,ia,indx) = sq2(inub,it,ip,ia,indx) + real(sqp*sqp)
           su2(inub,it,ip,ia,indx) = su2(inub,it,ip,ia,indx) + real(sup*sup)
           sv2(inub,it,ip,ia,indx) = sv2(inub,it,ip,ia,indx) + real(svp*svp)

           nums(inub,it,ip,ia,indx) = nums(inub,it,ip,ia,indx) + 1.

           aveinc(inub,it,ip,ia,indx) = aveinc(inub,it,ip,ia,indx) + real(abs(cost))

        end if

     end if
  end do

  if (iflag.eq.1) scount=scount+1.d0

400 continue

  return

end subroutine propagate


