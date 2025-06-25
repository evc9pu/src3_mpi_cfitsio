subroutine diskflux()

  ! calculates flux from disk
  !
  ! history:
  ! 2000/07/28 (mjw):  change filter on rtmp to .gt.1.d0

  use tts_mod
  use stokes_mod
  use constants

  implicit none

  real(8) :: rnorm,rnp,image_norm

  if (ifprint) print*,'in diskflux'

  ! total flux from disk
  rnp=dble(np)
  ! comparison with jon's code
  rnorm=rnp/dble(nmu)/dble(nph)/(ltot/lsol)
  image_norm=rnp/(ltot/lsol)
  ! rnorm=rnp/dble(nmu)
  ! image_norm=rnp
  if (ifprint) print*,'ltot check',ltot/lsol

  where(nums > 1.)
     si2 = sqrt((si2-si**2/nums)/(nums-1.)*nums)/si+1./nums
     sq2 = sqrt((sq2-sq**2/nums)/(nums-1.)*nums)/si
     su2 = sqrt((su2-su**2/nums)/(nums-1.)*nums)/si
     sv2 = sqrt((sv2-sv**2/nums)/(nums-1.)*nums)/si
     aveinc = aveinc / nums
  elsewhere
     si2 = 0.
     sq2 = 0.
     su2 = 0.
     sv2 = 0.
  end where

  si = si / real(rnorm)
  sq = sq / real(rnorm)
  su = su / real(rnorm)
  sv = sv / real(rnorm)

  if(ipeel==1) then

     where(numt > 1.)
        ti2 = sqrt((ti2*numt/ti**2-1.)/(numt-1.))
        tq2 = sqrt((tq2*numt/tq**2-1.)/(numt-1.))
        tu2 = sqrt((tu2*numt/tu**2-1.)/(numt-1.))
        tv2 = sqrt((tv2*numt/tv**2-1.)/(numt-1.))
     elsewhere
        ti2 = 0.
        tq2 = 0.
        tu2 = 0.
        tv2 = 0.
     end where

     ti = ti / real(image_norm)
     tq = tq / real(image_norm)
     tu = tu / real(image_norm)
     tv = tv / real(image_norm)

     if(output_imfilt) then
        image_b_i = image_b_i / real(image_norm)
        image_b_q = image_b_q / real(image_norm)
        image_b_u = image_b_u / real(image_norm)
        image_b_v = image_b_v / real(image_norm)
        image_b_d = image_b_d / real(image_norm)
        image_b_s = image_b_s / real(image_norm)
        image_b_e = image_b_e / real(image_norm)
        image_b_o = image_b_o / real(image_norm)
     end if

     if(output_imcube) then
        image_m_i = image_m_i / real(image_norm)
        image_m_q = image_m_q / real(image_norm)
        image_m_u = image_m_u / real(image_norm)
        image_m_v = image_m_v / real(image_norm)
     end if

     ! need to add nested image normalization

  end if

  ! check constants
  ! rnorm has factor of 2*pi because you made all the photons
  ! come out at one phi.
  ! rnorm=r2p*delnt*real(np)
  ! that was the value for the correct rnorm.  using that
  ! value, the flux at a given angle from the star, if there
  ! is not absorbing disk, is 1/4pi.  We will multiply the
  ! flux by 4pi so that it is normalized to the stellar flux.
  ! rnorm=0.5d0*dmu*real(np)/normstar
  ! rnorm=0.5d0*dmu*dble(np)/(ltot/lsol)

  ! multiply rnorm by 2 if you are doubling photons by
  ! symmetry through z axis (in diskim)
  ! rnorm=rnorm*2.d0
  ! fimage=0.d0
  ! print*,'normalizing flux image (not pol images) by ',rnorm
  ! whenever we deal with bins 1 and nmu, multiply by 2 because
  ! bin is half size of others.

  ! icenter=nx/2+1
  ! do it=1,nmu
  ! do ix=1,nx
  ! do iy=1,nx
  ! standard deviation of imagei array
  ! rtmp=dble(numi(ix,iy,it))
  ! if (image2(ix,iy,it).gt.1.d0) then
  ! print*,'rtmp',rtmp
  ! image2(ix,iy,it)=(image2(ix,iy,it)-imagei(ix,iy,it)**2/rtmp)
  ! 1           /rtmp/(rtmp-1.d0)
  ! end if
  ! if(ix.eq.icenter.and.iy.eq.icenter) then
  ! write(6,*) 'image before norm'
  ! 1           ,image(icenter,icenter,it),icenter,it
  ! end if
  ! if(image(ix,iy,it).lt.0.d0) then
  ! write(6,*) 'before normalizing'
  ! write(6,*) 'image lt 0',image(ix,iy,it),ix,iy,it
  ! end if
  ! if(it.eq.1.or.it.eq.nmu) then
  ! image(ix,iy,it)=image(ix,iy,it)/(rnorm)*2.d0
  ! else
  ! image(ix,iy,it)=image(ix,iy,it)/(rnorm)
  ! end if
  ! fimage=fimage+(image(ix,iy,it))/dble(nmu)
  ! if(ix.eq.icenter.and.iy.eq.icenter) then
  ! write(6,*) 'image(icenter,icenter,it),icenter,it'
  ! 1           ,image(icenter,icenter,it),icenter,it
  ! end if
  ! if(image(ix,iy,it).lt.0.d0) then
  ! write(6,*) 'image lt 0',image(ix,iy,it),ix,iy,it
  ! end if
  ! end do
  ! end do
  ! end do

  return

end subroutine diskflux
