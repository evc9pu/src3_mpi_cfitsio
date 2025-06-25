subroutine linepeel(ir,it,ip,Aul,abund,b0,jlevel,Eul,gammaul,ifout,ifdep)

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
  
  real(8) :: dt,Aul(5),temp,dens,abund(5),co_z,b0(5)
  real(8) :: jlevel(5),Eul(5),gammaul(5)
  real(8) :: vr,vthet,vphi,vlos,dnu,phinu,vx,vy,vz,sint,cost
  real(8) :: ju,jl,gu,gl,nup,nlo,nuul,lambdaul,jnu,kappanu,co_b,cogas
  real(8) :: densd,dense,denso
  real(8) :: nupd,nupe,nupo
  real(8) :: jnud,jnue,jnuo

  integer :: iface,iphot,it,ir,ip,id,idust2,inv,isp

  logical :: ifout,isnan,ifdep
  
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

        do isp=1,nspec

           if (ifdep) then
              call depletion_co(ir,it,ip,cogas)
           else
              cogas=1.d0
           end if

           dens=sum(densarr(ir,it,ip,:))/2.34d-24
           densd=(densarr(ir,it,ip,1)+densarr(ir,it,ip,2)+densarr(ir,it,ip,9))/2.34d-24
           dense=(densarr(ir,it,ip,3))/2.34d-24
           denso=(densarr(ir,it,ip,4)+densarr(ir,it,ip,10))/2.34d-24

           dens=dens*abund(isp)*cogas
           densd=densd*abund(isp)*cogas
           dense=dense*abund(isp)*cogas
           denso=denso*abund(isp)*cogas

           temp=tdave(ir,it,ip)*11605.d0
        
           co_b=sqrt(2.d0*k*temp/2.3d0/mH)/1.d5 ! in km/s

!           if (temp.gt.2000.d0) dens=0.d0
!           if (densarr(ir,it,ip,10).gt.0.d0) dens=0.d0

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
!           print*,vlos
        
           ju=jlevel(isp)
           jl=jlevel(isp)-1.d0
        
           gu=2.d0*ju+1.d0
           gl=2.d0*jl+1.d0
           
           co_z=sqrt(1.d0+(k*temp/b0(isp))**2)

           nup=dens/co_z*gu*exp(-b0(isp)*ju*(ju+1)/k/temp)
           nupd=densd/co_z*gu*exp(-b0(isp)*ju*(ju+1)/k/temp)
           nupe=dense/co_z*gu*exp(-b0(isp)*ju*(ju+1)/k/temp)
           nupo=denso/co_z*gu*exp(-b0(isp)*ju*(ju+1)/k/temp)

           nlo=dens/co_z*gl*exp(-b0(isp)*jl*(jl+1)/k/temp)
     
           nuul=Eul(isp)/h
           lambdaul=c/nuul

           do inv=1,nv
              dnu=(vnu(inv)-vlos)/3.d5*nuul
!              phinu=4.d0*gammaul(isp)/(16.d0*pi**2*dnu**2+gammaul(isp)**2)
              phinu=1.d0/sqrt(pi)/nuul*3.d5/co_b*exp(-(vnu(inv)-vlos)**2/co_b**2)
              jnu=nup*Aul(isp)/4.d0/pi*Eul(isp)*phinu
              jnud=nupd*Aul(isp)/4.d0/pi*Eul(isp)*phinu
              jnue=nupe*Aul(isp)/4.d0/pi*Eul(isp)*phinu
              jnuo=nupo*Aul(isp)/4.d0/pi*Eul(isp)*phinu
              kappanu=nlo*gu/gl*Aul(isp)/8.d0/pi*lambdaul**2*phinu*(1.d0-exp(-Eul(isp)/k/temp))
              if (isnan(jnu)) then
                 print*,'jnu is NAN'
                 print*,nup,Aul(isp),Eul(isp),phinu,dnu,vlos,vr,vphi,vthet
                 print*,dens,temp
                 stop
              end if
              Ivthin(inv,isp)=Ivthin(inv,isp)+jnu*dt*rstar*rsol
              Ivthind(inv,isp)=Ivthind(inv,isp)+jnud*dt*rstar*rsol
              Ivthine(inv,isp)=Ivthine(inv,isp)+jnue*dt*rstar*rsol
              Ivthino(inv,isp)=Ivthino(inv,isp)+jnuo*dt*rstar*rsol
              Ivthick(inv,isp)=Ivthick(inv,isp)+jnu*dt*rstar*rsol-Ivthick(inv,isp)*kappanu*dt*rstar*rsol
              if (Ivthick(inv,isp).lt.0.d0) Ivthick(inv,isp)=0.d0
           end do

        end do

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
  
end subroutine linepeel
  
