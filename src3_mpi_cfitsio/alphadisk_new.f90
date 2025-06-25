subroutine alphadisk(massc_in,rmincgs_in,mdotcgs_in,rmaxd_in,rinfall_in,fd_in,mdotcorecgs_in,sigmadotarr,alphah_in,alphal_in,alpha,sigmaarr)

  use constants
  use grid_mod
  use tts_mod, only : ifprint

  implicit none

  real(8) :: massc_in,rmincgs_in,mdotcgs_in,rmaxd_in,rinfall_in
  real(8) :: fd_in,mdotcorecgs_in,alpha,alphal_in,alphah_in
  real(8) :: Co_md,Co_Lw,Co_Ld,Co_md_l,Co_md_h,Co_md_t
  real(8) :: alphal,alphah,dalpha,alpha_t
  real(8) :: comdmax,comdmin

  real(8) :: sigmadotarr(nrg-1),sigmaarr(nrg-1),zmrarr(nrg-1)

  real(8),pointer :: alpharr(:),comdarr(:)

  integer :: nalpha,ial,ial_t,ialmax,ialmin,itera

  logical :: ifout,isnan,iffound

  ifout=.false.

  alphal=alphal_in
  alphah=alphah_in

  nalpha=1001

  itera=0

!!$5 call diskprofile(massc_in,rmincgs_in,mdotcgs_in,rmaxd_in,rinfall_in, &
!!$       & fd_in,mdotcorecgs_in,sigmadotarr,alphal,ifout,Co_md_l,Co_Lw, &
!!$       & Co_Ld,sigmaarr,zmrarr)
!!$
!!$  if (isnan(Co_md_l)) then 
!!$     alphal=alphal*1.01d0
!!$     goto 5
!!$  end if
!!$  if (alphal.ne.alphal_in) then
!!$     if (ifprint) print*,'adjust alphal from', alphal_in, 'to', alphal 
!!$  end if

  allocate(alpharr(nalpha),comdarr(nalpha))
  
5  dalpha=(log10(alphah)-log10(alphal))/dble(nalpha-1)

  do ial=1,nalpha
!     alpharr(ial)=10.d0**(log10(alphal)+dalpha*dble(ial-1))
     alpharr(ial)=10.d0**(log10(alphah)-dalpha*dble(ial-1))
  end do

  iffound=.false.
  comdmax=0.d0
  comdmin=1.d10
  ialmax=1
  ialmin=1
  do ial=1,nalpha
     
     alpha_t=alpharr(ial)

     sigmaarr=0.d0
     zmrarr=0.d0

     call diskprofile(massc_in,rmincgs_in,mdotcgs_in,rmaxd_in,rinfall_in, &
       & fd_in,mdotcorecgs_in,sigmadotarr,alpha_t,ifout,Co_md_t,Co_Lw, &
       & Co_Ld,sigmaarr,zmrarr)

     comdarr(ial)=Co_md_t

!     if (ifprint) print*,alpha_t,Co_md_t

     if (Co_md_t.gt.comdmax) then 
        comdmax=Co_md_t
        ialmax=ial
     end if

     if (Co_md_t.lt.comdmin) then
        comdmin=Co_md_t
        ialmin=ial
     end if

     if (ial.gt.1.and.(comdarr(ial)-fd_in)*(comdarr(ial-1)-fd_in).lt.0.d0.and.   &
          & abs((comdarr(ial)/comdarr(ial-1))-1).lt.0.1d0) then
!!$        alphal=alpharr(ial-1)
!!$        alphah=alpharr(ial)
        alphal=alpharr(ial)
        alphah=alpharr(ial-1)
        if (ifprint) print*,'alphal,alphah,co_md',alphal,alphah,comdarr(ial),comdarr(ial-1),fd_in
        iffound=.true.
        exit
     end if

!!$     if (isnan(Co_md_t)) then 
!!$        iffound=.true.
!!$     end if

  end do

  if (.not.iffound) then
     if (itera.le.0) then
        if (ifprint) print*,'check alpha range!!!'
        if (ifprint) print*,'max alpha, co_md',alpharr(ialmax),comdarr(ialmax),comdmax,fd_in
        if (ifprint) print*,'min alpha, co_md',alpharr(ialmin),comdarr(ialmin),comdmin,fd_in
        alphah=alphah_in*10.
        alphal=alphal_in/10.
        if (ifprint) print*,'changing alphal to', alphal
        if (ifprint) print*,'changing alphah to', alphah
        if (ifprint) print*,'search alpha one more time'
        itera=itera+1
        goto 5
     else
        print*,'check alpha range!!!'
        print*,'ial,nalpha',ial,nalpha
        print*,'max alpha, co_md',alpharr(ialmax),comdarr(ialmax),comdmax,fd_in
        print*,'min alpha, co_md',alpharr(ialmin),comdarr(ialmin),comdmin,fd_in
        stop
     end if
  end if

  alpha=sqrt(alphah*alphal)

  call diskprofile(massc_in,rmincgs_in,mdotcgs_in,rmaxd_in,rinfall_in, &
       & fd_in,mdotcorecgs_in,sigmadotarr,alpha,ifout,Co_md,Co_Lw, &
       & Co_Ld,sigmaarr,zmrarr)

  if (ifprint) print*,'alpha, Co_md/fd',alpha,Co_md/fd_in

  return

end subroutine alphadisk
