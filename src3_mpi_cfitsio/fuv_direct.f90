subroutine fuv_direct()

  use grid_mod
  use dust_mod
  use tts_mod
  use opacin_mod
  use constants
  use output_mod

  implicit none

  real(8) :: elow,ehigh,nulow,nuhigh,nutoev,dr,rad,dnufuv,cpusec,taudold,taudnew

  real(8),pointer :: nufuv(:),fuvarr(:,:,:),fuv(:),fuvstar(:),fuvunatt(:)
  real(8),pointer :: kapfuv(:,:),taud(:),fuvunattarr(:,:,:)

  integer :: nfuv,i,ir,it,ip,id

  call cpu_time(cpusec)
  print*,'cpu time', cpusec

  elow=6 !in ev
  ehigh=13.6 !in ev

  nutoev=6.58211928d-16*2.d0*pi

  nulow=elow/nutoev
  nuhigh=ehigh/nutoev

  print*,'nu low, high', nulow, nuhigh

  nfuv=101

  allocate(nufuv(nfuv),fuv(nfuv),fuvstar(nfuv),taud(nfuv),fuvunatt(nfuv))

  allocate(kapfuv(nopac,nfuv))

  do i=1,nfuv
     nufuv(i)=(nuhigh-nulow)*dble(i-1)/dble(nfuv-1)+nulow
  end do
  dnufuv=(nuhigh-nulow)/dble(nfuv-1)
  
!  print*,nufuv(1),nufuv(nfuv)

  do i=1,nfuv
     fuvstar(i)=2.d0*pi*h*(nufuv(i))**3/c**2/(exp(h*nufuv(i)/k/tstar)-1.d0)
!     print*,'fuvstar',i,fuvstar(i)
  end do

  print*,'fuv luminosity from star',sum(fuvstar*dnufuv)*4.d0*pi*(rstar*rsol)**2/lsol,sigt*tstar**4*4.d0*pi*(rstar*rsol)**2/lsol

  do i=1,nfuv
     call opacset(nufuv(i)*nutoev)
     do id=1,nopac
        kapfuv(id,i)=kappa(id)*kappav(id)*rstar*rsol
!        if (i.eq.1.and.id.le.4) print*,'kapfuv',id,kapfuv(id,i)/rstar/rsol
     end do
  end do
  
  allocate(fuvarr(nrg,ntg,npg),fuvunattarr(nrg,ntg,npg))

  fuvarr=0.d0

  do it=1,(ntg-1)/2

     if (mod(it,10).eq.0) print*,'it',it,(ntg-1)/2,thetavearr(it)/pi*180.d0
     
     do ip=1, npg-1

        taud=0.d0

        ir=1
        do while (ir.lt.nrg)

           dr=rarr(ir+1)-rarr(ir)
           rad=0.5d0*(rarr(ir+1)+rarr(ir))
        
           fuv=0.d0
           fuvunatt=0.d0
           taudold=taud(1)
           do i=1,nfuv
              do id=1,ndg+2
                 if (.not.is_sg(id).and.densarr(ir,it,ip,id).gt.0.d0) then
                    taud(i)=taud(i)+kapfuv(findopac(ir,it,ip,id),i)*densarr(ir,it,ip,id)*dr
                 end if
              end do
              fuvunatt(i)=fuvstar(i)/rad**2
              fuv(i)=fuvstar(i)/rad**2*exp(-taud(i))
           end do
           taudnew=taud(1)

           fuvarr(ir,it,ip)=sum(fuv*dnufuv)/5.29d-14/c
           fuvunattarr(ir,it,ip)=sum(fuvunatt*dnufuv)/5.29d-14/c

!!$           if (it.eq.400.and.ip.eq.1) then
!!$              print*,ir,rad/autors,sum(densarr(ir,it,ip,:)),dr/autors,taud(1),fuvarr(ir,it,ip),fuvunattarr(ir,it,ip)
!!$           end if

           if (fuvarr(ir,it,ip).lt.1.d-2) exit

           ir=ir+1
        end do

!        print*,ir,rad/autors,fuvarr(ir-1,it,ip),fuvarr(ir,it,ip),taudold,taudnew,sum(densarr(ir-1,it,ip,:)),dr/autors,fuvunattarr(ir-1,it,ip),fuvunattarr(ir,it,ip)

     end do
  end do

  do it=(ntg+1)/2,ntg-1
     do ip=1,npg-1
        do ir=1,nrg-1
           fuvarr(ir,it,ip)=fuvarr(ir,ntg-it,ip)
        end do
     end do
  end do

  if (ifprint) call output_grid('fuv_direct',fuvarr)

  call cpu_time(cpusec)
  print*,'cpu time', cpusec

end subroutine fuv_direct
