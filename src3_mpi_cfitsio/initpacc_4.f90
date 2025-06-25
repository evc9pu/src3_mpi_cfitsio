! ********************************************************

subroutine initpacc(ii,nphot,iter)

  ! initialize variables for new accretion photon

  use tts_mod
  use grid_mod
  use stokes_mod
  use taunum_mod
  use dust_mod
  use random
  use opacin_mod
  use output_mod
  use constants

  implicit none

  real(8) :: winv,xi
  real(8) :: c2p,s2p,cp,sp
  integer :: ii(3),ir,it,ip,id,nphot,idust2,iter
  integer :: try

  ! find photon position

  ! only emitting photons in disk.  adjust accretion
  ! luminosity accordingly.
  xi=ran()
  call locate(Ldint,nrg-1,xi,ir)
  rp=ravearr(ir)
  
11  call rantheta(sp,cp,s2p,c2p)
  xp=rp*cp
  yp=rp*sp
  !20090731, weight zp by different dust-type masses in the disk
  zp = hdisk(ir)*gasdev()

  sip=1.d0
  sqp=0.d0
  sup=0.d0
  svp=0.d0

  rsq=xp*xp+yp*yp+zp*zp
  rtot=sqrt(rsq)
  call locate(r2arr,nrg,rsq,ir)

  ! if (ntg.gt.1) then
  cosb=zp/rtot
  ! sinb=sqrt(1.d0-cosb**2)
  call locate(costarr,ntg,cosb,it)
  ! else
  ! it=1
  ! end if
  ! if (npg.gt.1) then
  lp=atan2(yp,xp)
  if (lp.lt.0.d0) lp=lp+r2p
  call locate(phiarr,npg,lp,ip)
  ! if (ip.gt.1) then
  ! print*,'error in 2-D code; ip should be 1',ip
  ! end if
  ! else
  ! ip=1
  ! end if
  ii(1)=ir
  ii(2)=it
  ii(3)=ip

  ! what is id...
  ! randomly sample thermal dust type.
  if(sum(densarr(ir,it,ip,:), mask=is_thermal).gt.0.d0) then
!     print*,'1',sum(densarr(ir,it,ip,:),mask=is_thermal)
     call select_dust(ir,it,ip,is_thermal,idust2)
  else
!     print*,'2'
!     idust2 = -1 ! indicates emit from gas
     goto 11
  end if

  ! SHOULD WE ALLOW ENVELOPE TYPES TO BE SELECTED?

  !     tabulate which emission process is chosen as a function of radius
  if(idust2 > 0) then
     if (is_sg(idust2)) then
        dsg(ir,it,ip)=dsg(ir,it,ip)+1
     else
        dbig(ir,it,ip)=dbig(ir,it,ip)+1
     end if
  end if

  if (diffus(ir,it,ip)) then
     nabs(ir,it,ip)=nabs(ir,it,ip)+1
  else if(idust2 > 0) then
     if (.not.is_sg(idust2).and.densarr(ir,it,ip,idust2).gt.0.d0) &
          & nabs(ir,it,ip)=nabs(ir,it,ip)+1
  end if

  call emit_common(ir,it,ip,idust2,iter,ii,nphot)

  return
end subroutine initpacc


! **********************************************************
