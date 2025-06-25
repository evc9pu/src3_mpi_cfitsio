subroutine opthinpeel(ir,it,ip,Aul,abund,b0,jlevel,Eul,gammaul,ifreverse)

  use tts_mod
  use grid_mod
  use taunum_mod
  use opacin_mod
  use messages
  use dust_mod
  use closest_wall, only : reset_t,find_next_wall
  use constants
  use lineint_mod

  implicit none
  
  character(len=20),parameter :: routine = 'opthinpeel'
  
  real(8) :: dt,Aul,temp,dens,abund,co_z,b0,jlevel,Eul,gammaul
  real(8) :: vr,vthet,vphi,vlos,dnu,phinu,vx,vy,vz

  integer :: iface,iphot,it,ir,ip,id,idust2,inv

  logical :: ifreverse
  
  iphot=1

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
     
     if (dt.lt.0.d0) then
        print*,'problem with find_wall! did not find surface! '
        print*,'in tauint,ir,rtot2,r2arr(ir)',ir,rsq,r2arr(ir)
        print*,'ux,uy,uz,iphot',ux,uy,uz,iphot
        print*,'xtot,ytot,ztot',xp,yp,zp
        stop
     end if
     
     dens=sum(densarr(ir,it,ip,:))/2.34d-24
     dens=dens*abund

     temp=tdave(ir,it,ip)*11605.d0

     if (temp.gt.2000.d0) dens=0.d0

     co_z=sqrt(1.d0+(k*temp/b0)**2)

     dens=dens/co_z*(2.d0*jlevel+1.d0)*exp(-b0*jlevel*(jlevel+1)/k/temp)

     vr=vrarr(ir,it,ip)
     vthet=vthetarr(ir,it,ip)
     vphi=0.d0

     vx=vr*xp/rtot+vthet*zp/rtot*xp/rtot-vphi*yp/rtot
     vy=vr*yp/rtot+vthet*zp/rtot*yp/rtot+vphi*xp/rtot
     vz=vr*zp/rtot-vthet*sqrt(xp**2+yp**2)/rtot

     vlos=vx*ux+vy*uy+vz*uz

     if (ifreverse) vlos=-vlos

     do inv=1,nv
        dnu=(vnu(inv)-vlos)/3.d5*Eul/h
        phinu=4.d0*gammaul/(16.d0*pi**2*dnu**2+gammaul)
        Iv(inv)=Iv(inv)+dens*Aul/4.d0/pi*Eul*phinu
     end do
          
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
     
     ! if photon exists grid, exit integration
     if (ir == nrg) then
        exitflag=1
        return
     end if
     
     ! if photon intersects star, re-emit
     if (ir == 0) then
        print*,'opthinpeel meets star!!'
        stop
     end if
       
  end do
  
end subroutine opthinpeel
  
