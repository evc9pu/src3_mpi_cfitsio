subroutine testgrid(xp,yp,zp,rsq,rtot,ir,it,ip,on_wall,iphot,routine)

  ! ******************************************************

  ! test to see if photon position corresponds to grid
  ! position, ii(1),ii(2),ii(3)

  ! input:
  ! xp,yp,zp,rsq photon position
  ! ii(3), grid indices
  ! routine, name of subroutine called from
  ! iphot, photon number

  ! output:
  ! ii(3), updated if necessary.

  ! ******************************************************

  use grid_mod
  use messages
  use constants
  use random

  implicit none

  real(8) ::  xp,yp,zp,rsq,cosb,rtot
  integer :: ir,it,ip,iphot
  character(len=*) :: routine
  real(8) :: check

  real(8) :: rfrac,tfrac,pfrac

  real(8),parameter :: tol = 1.d-4
  real(8) :: phi,dphi1,dphi2

  logical(1) :: phierr
  logical(1),intent(inout) :: on_wall
  logical(1) :: reset

  cosb=zp/rtot
  cosb=check(cosb,routine)

  phi = atan2(yp,xp)
  if(phi < 0.d0) phi = phi + r2p

  reset = .false.

  if((rsq - r2arr(1))/(r2arr(2) - r2arr(1)) < -tol .or. (rsq - r2arr(nrg-1))/(r2arr(nrg) - r2arr(nrg-1)) > 1.d0+tol) then
     print *, "WARNING - photon is not inside grid - resetting r to random position"
     rtot = (ran() * (rarr(nrg)**3. - rarr(1)**3.) + rarr(1)**3.)**(1./3.)
     reset = .true.
  end if

  if((cosb - costarr(1))/(costarr(2) - costarr(1)) < -tol .or. &
       & (cosb - costarr(ntg-1))/(costarr(ntg) - costarr(ntg-1)) > 1.d0+tol) then
     print *, "WARNING - photon is not inside grid - resetting theta to random position"
     cosb = (ran() * (costarr(ntg) - costarr(1)) + costarr(1))
     reset = .true.
  end if

  if((phi - phiarr(1))/(phiarr(2) - phiarr(1)) < -tol .or. (phi - phiarr(npg-1))/(phiarr(npg) - phiarr(npg-1)) > 1.d0+tol) then
     print *, "WARNING - photon is not inside grid - resetting phi to random position"
     phi = ran() * (phiarr(npg) - phiarr(1)) + phiarr(1)
     reset = .true.
  end if

  if(reset) then
     xp = rtot * sqrt(1.-cosb**2.) * cos(phi)
     yp = rtot * sqrt(1.-cosb**2.) * sin(phi)
     zp = rtot * cosb
     rsq = rtot * rtot
     on_wall = .false.
  end if

  if(ir.lt.1.or.ir.gt.nrg) then
     print *, "WARNING - ir is out of bounds, resetting"
     call locate(r2arr,nrg,rsq,ir)
  end if

  if(it.lt.1.or.it.gt.ntg) then
     print *, "WARNING - it is out of bounds, resetting"
     call locate(costarr,ntg,cosb,it)
  end if

  if(ip.lt.1.or.ip.gt.npg) then
     print *, "WARNING - ip is out of bounds, resetting"
     call locate(phiarr,npg,phi,ip)
  end if

  phierr = .false.

  rfrac = (rsq-r2arr(ir))/(r2arr(ir+1)-r2arr(ir))
  tfrac = (cosb-costarr(it))/(costarr(it+1)-costarr(it))
  pfrac = (phi-phiarr(ip))/(phiarr(ip+1)-phiarr(ip))

  if(pfrac < -tol .or. pfrac > 1.d0+tol) then
     dphi1 = phi - phiarr(ip)
     if(dphi1 > pi) dphi1 = dphi1 - r2p
     if(dphi1 < -pi) dphi1 = dphi1 + r2p
     dphi2 = phi - phiarr(ip+1)
     if(dphi2 > pi) dphi2 = dphi2 - r2p
     if(dphi2 < -pi) dphi2 = dphi2 + r2p
     phierr = abs(dphi1) > tol .and. abs(dphi2) > tol
  end if

  if(rfrac < -tol .or. rfrac > 1.d0+tol) then
     print *,'WARNING - r does not match ir, resetting ir, rfrac=',rfrac, ir, r2arr(ir), r2arr(ir+1), rsq
     call locate(r2arr,nrg,rsq,ir)
     if(on_wall) print *,'EXTRA WARNING: photon is supposed to be on wall, resetting might not have worked'
  end if

  if(tfrac < -tol .or. tfrac > 1.d0+tol) then
     if (iwarn.eq.1) print *,'WARNING - theta does not match it, resetting it, tfrac=',tfrac, it, costarr(it), costarr(it+1), cosb
     call locate(costarr,ntg,cosb,it)
     if(iwarn.eq.1.and.on_wall) print *,'EXTRA WARNING: photon is supposed to be on wall, resetting might not have worked'
  end if

  if(phierr) then
     print *,'WARNING - phi does not match ip, resetting ip, pfrac=',pfrac
     print*,'ir,it,ip,rtot',ir,it,ip,rtot
     print*,'cosb,phi,iphot',cosb,phi,iphot
     call locate(phiarr,npg,phi,ip)
     if(on_wall) print *,'EXTRA WARNING: photon is supposed to be on wall, resetting might not have worked'
  end if

  return

end subroutine testgrid


