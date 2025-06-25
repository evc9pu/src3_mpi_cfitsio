subroutine alphadisk(massc_in,rmincgs_in,mdotcgs_in,rmaxd_in,rinfall_in,fd_in,mdotcorecgs_in,sigmadotarr,alphah_in,alphal_in,alpha,sigmaarr)

  use constants
  use grid_mod
  use tts_mod, only : ifprint

  implicit none

  real(8) :: massc_in,rmincgs_in,mdotcgs_in,rmaxd_in,rinfall_in
  real(8) :: fd_in,mdotcorecgs_in,alpha,alphal_in,alphah_in
  real(8) :: Co_md,Co_Lw,Co_Ld,Co_md_l,Co_md_h,Co_md_t
  real(8) :: alphal,alphah,dalpha,alpha_t
  real(8) :: comdmax,comdmin,comdl,comdh
  real(8) :: dis

  real(8) :: sigmadotarr(nrg-1),sigmaarr(nrg-1),zmrarr(nrg-1)

  real(8),pointer :: alpharr(:),comdarr(:)

  integer :: nalpha,ial,ial_t,ialmax,ialmin,itera

  logical :: ifout,isnan,iffound

  ifout=.false.

  alphal=alphal_in
  alphah=alphah_in

  nalpha=1001

5 call diskprofile(massc_in,rmincgs_in,mdotcgs_in,rmaxd_in,rinfall_in, &
       & fd_in,mdotcorecgs_in,sigmadotarr,alphal,ifout,Co_md_l,Co_Lw, &
       & Co_Ld,sigmaarr,zmrarr)

  if (isnan(Co_md_l)) then 
     alphal=alphal*1.01d0
     goto 5
  end if
  if (alphal.ne.alphal_in) then
     if (ifprint) print*,'adjust alphal from', alphal_in, 'to', alphal 
  end if

  allocate(alpharr(nalpha),comdarr(nalpha))

  itera=0
10  dalpha=(log10(alphah)-log10(alphal))/dble(nalpha-1)

  do ial=1,nalpha
     alpharr(ial)=10.d0**(log10(alphal)+dalpha*dble(ial-1))
  end do

  iffound=.false.
  comdmax=0.d0
  comdmin=1.d10
  ialmax=1
  ialmin=1

!!$  if (ifprint) then
!!$     print*,massc_in,rmincgs_in,mdotcgs_in
!!$     print*,rmaxd_in,rinfall_in,fd_in
!!$     print*,mdotcorecgs_in
!!$     print*,maxval(sigmadotarr)
!!$  end if

  do ial=1,nalpha
     
     alpha_t=alpharr(ial)

     sigmaarr=0.d0
     zmrarr=0.d0

     call diskprofile(massc_in,rmincgs_in,mdotcgs_in,rmaxd_in,rinfall_in, &
       & fd_in,mdotcorecgs_in,sigmadotarr,alpha_t,ifout,Co_md_t,Co_Lw, &
       & Co_Ld,sigmaarr,zmrarr)

     comdarr(ial)=Co_md_t

!!$     if (ial.gt.1.and.comdarr(ial)/comdarr(ial-1)-1.d0.gt.1.d-3) then
!!$        print*,''
!!$        print*,'alpha,co_md',ial-1,alpharr(ial-1),comdarr(ial-1)
!!$        print*,'alpha,co_md',ial,alpharr(ial),comdarr(ial)
!!$     end if

     if (Co_md_t.gt.comdmax) then 
        comdmax=Co_md_t
        ialmax=ial
     end if

     if (Co_md_t.lt.comdmin) then
        comdmin=Co_md_t
        ialmin=ial
     end if

     dis=log10(alphah)
!     if (ifprint.and.itera.gt.0) print*,comdarr(ial),comdarr(ial-1),fd_in
     if (ial.gt.1.and.(comdarr(ial)-fd_in)*(comdarr(ial-1)-fd_in).lt.0.d0) then
          if (abs((comdarr(ial)/comdarr(ial-1))-1).lt.0.1d0) then
             alphal=alpharr(ial-1)
             alphah=alpharr(ial)
             if (ifprint) print*,'alphal,alphah,co_md',alphal,alphah,comdarr(ial-1),comdarr(ial),fd_in
             iffound=.true.
             exit
!!$          else
!!$             if (abs(log10(alpharr(ial))-log10(alpha)).lt.dis) then
!!$                dis=abs(log10(alpharr(ial))-log10(alpha))
!!$                alphal=alpharr(ial-1)
!!$                alphah=alpharr(ial)
!!$                comdl=comdarr(ial-1)
!!$                comdh=comdarr(ial)
!!$             end if
          end if
     end if

  end do

  if (.not.iffound) then
     print*,'check alpha range!!!'
     print*,'ial,nalpha',ial,nalpha
     print*,'max alpha, co_md',alpharr(ialmax),comdarr(ialmax),comdmax,fd_in
     print*,'min alpha, co_md',alpharr(ialmin),comdarr(ialmin),comdmin,fd_in
     stop
!!$     if (ifprint) then
!!$        print*,'check alpha range!!!'
!!$        print*,'alphal,alphah,alpha',alphal,alphah,alpha,comdl,comdh
!!$     end if
!!$     alphal=alpha
!!$     alphah=alpha
!     itera=itera+1
!     goto 10
  end if

!!$  if (comdarr(1).lt.fd_in.or.comdarr(nalpha).gt.fd_in) then
!!$     print*,'alpha range wrong',alphal,alphah,comdarr(1),comdarr(nalpha),fd_in
!!$     stop
!!$  end if
!!$
!!$  call locate(comdarr,nalpha,fd_in,ial_t)
!!$  if (comdarr(ial_t).lt.fd_in) then
!!$     alphal=alpharr(ial_t-1)
!!$     alphah=alpharr(ial_t)
!!$  else
!!$     alphal=alpharr(ial_t)
!!$     alphah=alpharr(ial_t+1)
!!$  end if

  alpha=sqrt(alphah*alphal)

  call diskprofile(massc_in,rmincgs_in,mdotcgs_in,rmaxd_in,rinfall_in, &
       & fd_in,mdotcorecgs_in,sigmadotarr,alpha,ifout,Co_md,Co_Lw, &
       & Co_Ld,sigmaarr,zmrarr)

  if (ifprint) print*,'alpha, Co_md/fd',alpha,Co_md/fd_in

  return

end subroutine alphadisk
