subroutine initarr_final(first)

  use tts_mod
  use grid_mod
  implicit none

  logical,intent(in) :: first

  integer :: i,id
  integer :: ir,it,ip
  real(8) :: dump1(nrg-1),dump2(ntg-1),dump3(npg-1)
  real(8) :: tdave_temp(nrg-1,ntg-1,npg-1)

  if(first) then

     call grid_allocate()

     ! cos(theta) array
     allocate(u(nmu))
     dmu=2.d0/dble(nmu)
     do i=1,nmu
        u(i)=1.d0-dmu*(dble(i)-0.5d0)
        if (ifprint) print '("mu_",I0," = ",F10.5)',i,u(i)
     end do

     ! SEDs
     allocate(si(nfreq,nmu,nph,nap,no))
     allocate(sq(nfreq,nmu,nph,nap,no))
     allocate(su(nfreq,nmu,nph,nap,no))
     allocate(sv(nfreq,nmu,nph,nap,no))

     ! SED uncertainties
     allocate(si2(nfreq,nmu,nph,nap,no))
     allocate(sq2(nfreq,nmu,nph,nap,no))
     allocate(su2(nfreq,nmu,nph,nap,no))
     allocate(sv2(nfreq,nmu,nph,nap,no))

     ! Number of photons in SEDs
     allocate(nums(nfreq,nmu,nph,nap,no))

     ! Other arrays
     allocate(aveinc(nfreq,nmu,nph,nap,no))
     allocate(tauenv(nfreq,nmu,nph))

     densarr=0.d0
!     massarr=0.d0
     !tdust2=3.d0/11605.d0
!     tdust2=0.1d0/11605.d0
!     tdave=0.1d0/11605.d0
     diffdir=0
     diffus=.false.
!     nodust=.false.
     vcell=0.d0
     imean=0

  end if
  
  x0arr=0.d0
  vrarr=0.d0
  vthetarr=0.d0
  vzarr=0.d0
  vcynarr=0.d0
  divvarr=0.d0
  adiarr=0.d0
  advarr=0.d0
  
  fuvmean=0.d0

  killed = 0
  
  dbig=0
  dsg=0
  jmean=0.d0

  si = 0.0
  sq = 0.0
  su = 0.0
  sv = 0.0

  si2 = 0.0
  sq2 = 0.0
  su2 = 0.0
  sv2 = 0.0

  nums = 0.0

  aveinc = 0.0

  tdave=0.d0
!  open(unit=14,file='tdave_f.unf',status='old',form='unformatted',convert='little_endian')
  open(unit=14,file='tdave_f.unf',status='old',form='unformatted')
  read(14) ir,it,ip,id
  print*,''
  print*,ir,it,ip,id
  print*,''
  read(14) dump1
  read(14) dump2
  read(14) dump3
  read(14) tdave_temp
  close(14)

  do ir=1,nrg-1
     do it=1,ntg-1
        do ip=1,npg-1
           tdave(ir,it,ip)=tdave_temp(ir,it,ip)
        end do
     end do
  end do

  tdave=tdave/11605.d0

!  do id=1,ndg+2
!     where(diffus)
!        tdust(:,:,:,id) = 0.1d0/11605.d0
!     elsewhere
!        tdust(:,:,:,id)=tdust2(:,:,:,id)
!     end where
!  end do

  nabs=0
  dtauabs=0.d0

end subroutine initarr_final

subroutine initimages()

  use type_nested_image
  use tts_mod
  implicit none

  integer :: peelid
  real :: freq_min,freq_max

  ! SEDs
  ! (n_freq,n_apmax,n_output)

  allocate(ti(nfreq,npeel,nap,no)) ; ti=0.0
  allocate(tq(nfreq,npeel,nap,no)) ; tq=0.0
  allocate(tu(nfreq,npeel,nap,no)) ; tu=0.0
  allocate(tv(nfreq,npeel,nap,no)) ; tv=0.0

  ! SED uncertainties
  ! (n_freq,n_apmax,n_output)

  allocate(ti2(nfreq,npeel,nap,no)) ; ti2=0.0
  allocate(tq2(nfreq,npeel,nap,no)) ; tq2=0.0
  allocate(tu2(nfreq,npeel,nap,no)) ; tu2=0.0
  allocate(tv2(nfreq,npeel,nap,no)) ; tv2=0.0

  allocate(numt(nfreq,npeel,nap,no)) ; numt=0.0

  ! Images - real(8), dimension(npeel,n_x,n_x,n_bands)

  if(output_imfilt) then
     allocate(image_b_i(npeel,nxhst,nxhst,nbnd)) ; image_b_i=0.0
     allocate(image_b_q(npeel,nxhst,nxhst,nbnd)) ; image_b_q=0.0
     allocate(image_b_u(npeel,nxhst,nxhst,nbnd)) ; image_b_u=0.0
     allocate(image_b_v(npeel,nxhst,nxhst,nbnd)) ; image_b_v=0.0
     allocate(image_b_d(npeel,nxhst,nxhst,nbnd)) ; image_b_d=0.0
     allocate(image_b_s(npeel,nxhst,nxhst,nbnd)) ; image_b_s=0.0
     allocate(image_b_e(npeel,nxhst,nxhst,nbnd)) ; image_b_e=0.0
     allocate(image_b_o(npeel,nxhst,nxhst,nbnd)) ; image_b_o=0.0
  end if

  ! Monochromatic:

  if(output_imcube) then
     allocate(image_m_i(npeel,nxhst,nxhst,nfreq)) ; image_m_i=0.0
     allocate(image_m_q(npeel,nxhst,nxhst,nfreq)) ; image_m_q=0.0
     allocate(image_m_u(npeel,nxhst,nxhst,nfreq)) ; image_m_u=0.0
     allocate(image_m_v(npeel,nxhst,nxhst,nfreq)) ; image_m_v=0.0
  end if

  ! Nested images

  if(output_imsparse) then

     allocate(image_ms_i(npeel))
     allocate(image_ms_q(npeel))
     allocate(image_ms_u(npeel))
     allocate(image_ms_v(npeel))

     do peelid=1,npeel
        freq_min = log10(real(numin)*2.9979e14/1.2398e0)
        freq_max = log10(real(numax)*2.9979e14/1.2398e0)
        call nested_image_setup(image_ms_i(peelid),sparse_dim,sparse_rmax*real(autors),&
        &sparse_rmin,nfreq,freq_min,freq_max,sparse_factor)
        call nested_image_setup(image_ms_q(peelid),sparse_dim,sparse_rmax*real(autors),&
        &sparse_rmin,nfreq,freq_min,freq_max,sparse_factor)
        call nested_image_setup(image_ms_u(peelid),sparse_dim,sparse_rmax*real(autors),&
        &sparse_rmin,nfreq,freq_min,freq_max,sparse_factor)
        call nested_image_setup(image_ms_v(peelid),sparse_dim,sparse_rmax*real(autors),&
        &sparse_rmin,nfreq,freq_min,freq_max,sparse_factor)
     end do

  end if

end subroutine initimages
