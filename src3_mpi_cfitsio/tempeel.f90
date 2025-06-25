subroutine tempeel(ir,it,ip,ncrit,ifout)

  use tts_mod
  use grid_mod
  use taunum_mod
  use opacin_mod
  use messages
  use dust_mod
  use closest_wall, only : reset_t,find_next_wall
  use constants
  use tempmap_mod

  implicit none
  
  character(len=20),parameter :: routine = 'tempeel'
  
  real(8) :: dt,dens,temp,ncrit,vr,vthet,vphi,vlos,vx,vy,vz,sint,cost

  integer :: iface,iphot,it,ir,ip,id,idust2

  logical :: ifout,isnan
  
  iphot=1

  on_wall=.False.

  do

     call testgrid(xp,yp,zp,rsq,rtot,ir,it,ip,on_wall,iphot,routine)

     call reset_t()
     
     ! find shortest non-negative distance to two constant radius surfaces
     call radface(ux,uy,uz,xp,yp,zp,r2arr(ir),r2arr(ir+1))
     
     ! find shortest positive distance to two constant theta surfaces
     if (ntg.gt.2) call thetaface(ux,uy,uz,xp,yp,zp,tan2arr(it),tan2arr(it+1))
     
     ! find shortest non-negative distance to two constant phi surfaces
     if (npg.gt.2.and.ntg.gt.2) call phiface(aarr(ip),barr(ip),0.d0,0.d0,aarr(ip+1),barr(ip+1),0.d0,0.d0,ux,uy,uz,xp,yp,zp)
     
     call find_next_wall(on_wall,dt,iface)
     
!     print*,ir,it,ip,xp,yp,zp,rtot,ux,uy,uz,dt

     if (dt.lt.0.d0) then
        print*,'problem with find_wall! did not find surface! '
        print*,'in tauint,ir,rtot2,r2arr(ir)',ir,rsq,r2arr(ir)
        print*,'ux,uy,uz,iphot',ux,uy,uz,iphot
        print*,'xtot,ytot,ztot',xp,yp,zp
        stop
     end if
     
     if (ifout) then

     dens=sum(densarr(ir,it,ip,:))/2.34d-24/2.

     temp=tdave(ir,it,ip)*11605.d0

     if (dens.lt.ncrit) dens=0.d0

     vr=vrarr(ir,it,ip)
     if (isnan(vr)) then
        vr=0.d0
     end if
     vthet=vthetarr(ir,it,ip)
     
     cost=zp/rtot
     sint=sqrt(1.d0-cost**2)
     
     vphi=0.d0
     if (densarr(ir,it,ip,1)+densarr(ir,it,ip,2)+densarr(ir,it,ip,9).gt.0.d0) then
        vphi=sqrt(gn*massc*msol/rtot/rstar/rsol)
     else if (densarr(ir,it,ip,3).gt.0.d0.and.rtot.lt.sonicp*autors) then
        vphi=sqrt(gn*massc*4.d0/3.d0*msol/rtot/rstar/rsol)*sqrt(1.d0-windmu0**2) &
             & /sint*sqrt(1.d0-cost/windmu0)
     else if (densarr(ir,it,ip,4).gt.0.d0) then
        vphi=10.d0*sqrt(gn*massc*msol/(rmaxd*(1.d0-windmu0**2))**x0arr(ir,it,ip)/rstar/rsol)
     end if
     
     vx=vr*xp/rtot+vthet*zp/rtot*xp/rtot-vphi*yp/rtot
     vy=vr*yp/rtot+vthet*zp/rtot*yp/rtot+vphi*xp/rtot
     vz=vr*zp/rtot-vthet*sqrt(xp**2+yp**2)/rtot
     
     vlos=vx*ux+vy*uy+vz*uz
     vlos=vlos/1.d5 !in km/s

     vrhor=vrhor+vlos*dens*dt
     Trhor=Trhor+temp*dens*dt
     rhor=rhor+dens*dt

     end if

     select case(iface)
     case(1)
        ir=ir-1
     case(2)
        ir=ir+1
     case(3)
        it=it-1
     case(4)
        it=it+1
     case(5)
        ip=ip-1
     case(6)
        ip=ip+1
     case default
        stop "ERROR: invalid iface"
     end select
     
     ! account for crossing the phi=0 plane:
     if (ip == 0)   ip = npg - 1
     if (ip == npg) ip = 1
     
     xp=xp+dt*ux
     yp=yp+dt*uy
     zp=zp+dt*uz
     on_wall = .true.
     rsq=xp*xp+yp*yp+zp*zp
     rtot=sqrt(rsq)
     
     !     print*,ir,it,ip,xp,yp,zp,rtot
     
  ! if photon exists grid, exit integration
     if (ir == nrg-1) then
        exitflag=1
        return
     end if
     
     ! if photon intersects star, re-emit
     if (ir == 0) then
        print*,'opthinpeel meets star!!'
        stop
     end if
     
  end do
  
end subroutine tempeel

