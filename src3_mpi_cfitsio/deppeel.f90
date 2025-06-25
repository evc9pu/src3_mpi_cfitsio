subroutine deppeel(ir,it,ip,ifout)

  use tts_mod
  use grid_mod
  use taunum_mod
  use opacin_mod
  use messages
  use dust_mod
  use closest_wall, only : reset_t,find_next_wall
  use constants
  use depint_mod

  implicit none

  character(len=20),parameter :: routine = 'opthinpeel'

  real(8) :: S0,Rdg,crosssec,mu,Eco,abunddep,nu,Agr,Nb,crrate
  real(8) :: co_on,co_off,gas_crit,tau_on,tau_off1,tau_off2,tau_offcr
  real(8) :: dens,temp,gas_t,mass_t,yr,nh_t,t_t,dt

  integer :: iphot,iface,ir,it,ip,id,idust2,idep

  logical :: ifout,isnan

  iphot=1

  on_wall=.False.

  yr=365.25d0*24.d0*3600.d0

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

        do idep=1,ndep

           S0=S0_arr(idep)
           Rdg=Rdg_arr(idep)
           crosssec=crosssec_arr(idep)
           mu=mu_arr(idep)
           Eco=Eco_arr(idep)
           abunddep=abunddep_arr(idep)
           nu=nu_arr(idep)
           Agr=Agr_arr(idep)
           Nb=Nb_arr(idep)
           crrate=crrate_arr(idep)

           co_on=S0*Rdg*crosssec*sqrt(8.*k/pi/mH)
           co_off=Agr*Rdg
           gas_crit=1.-Nb*Rdg/abunddep

!           dens=densarr(ir,it,ip,3)/2.34d-24 ! density only in envelope
           dens=(densarr(ir,it,ip,1)+densarr(ir,it,ip,2)+densarr(ir,it,ip,3)+densarr(ir,it,ip,4))/2.34d-24

           temp=tdave(ir,it,ip)*11605.d0

           mass_t=dens*dt

           T_t=temp
           nH_t=dens

           if (dens.gt.0.d0) then

              tau_on=1./(co_on*nH_t*sqrt(T_t/mu))/yr
              tau_off1=1./(co_off*nu*exp(-Eco/T_t)/abunddep)/yr
              tau_off2=1./(co_off*nu*exp(-Eco/T_t)/Nb/Rdg)/yr
              tau_offcr=3.3e6/(crrate/1.e-17)*exp(Eco/70.)/yr
              gas_t=(tau_on/tau_off1)*(tau_offcr+tau_off1)/(tau_offcr+tau_on)
              if (gas_t.gt.gas_crit) then
                 gas_t=(tau_on*tau_off2+tau_on*tau_offcr)/(tau_on*tau_off2+tau_on*tau_offcr+tau_offcr*tau_off2)
                 if (gas_t.lt.gas_crit) then
                    print*,'could not find proper gas_t in eq.'
                    stop
                 end if
              end if
              gascol(idep)=gascol(idep)+gas_t*mass_t
              totcol(idep)=totcol(idep)+mass_t
              
           else

              gascol(idep)=gascol(idep)+0.d0
              totcol(idep)=totcol(idep)+0.d0

           end if

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

end subroutine deppeel
