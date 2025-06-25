module tauint_mod
  implicit none
contains

  subroutine tauint_cell(iphot,tsum,ir,it,ip,idust2)

    use tts_mod
    use grid_mod
    use taunum_mod
    use opacin_mod
    use messages
    use dust_mod
    use random
    use closest_wall, only : reset_t,find_next_wall
    use constants
    implicit none

    character(len=20),parameter :: routine = 'tauint_3d'

    real(8) :: tsum,dt,opac,dtau,ds,xran
        integer :: iface,iphot,it,ir,ip,idust2,id

    ! set optical depth sum to zero
!    tsum=0.d0

       call testgrid(xp,yp,zp,rsq,rtot,ir,it,ip,routine,iphot,r2p,uz)

       idust2=dustarr(ir,it,ip)

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

       !     optical depth to edge of cell
       if (tdust(ir,it,ip).lt.1599.d0/11605.d0) then
          opac=kapd(idust2)*densarr(ir,it,ip)
       else
          opac=0.d0
       end if

       dtau=opac*dt

       if ((tsum+dtau).gt.tau) then
          ds=(tau-tsum)/opac
          dtauabs(ir,it,ip)=dtauabs(ir,it,ip)+ds*kapd(idust2)*densarr(ir,it,ip)*(1.d0-albedo(idust2))
          xp=xp+ds*ux
          yp=yp+ds*uy
          zp=zp+ds*uz
          rsq=xp**2+yp**2+zp**2
          rtot=sqrt(rsq)
          on_wall = .false.
          ! we are done with optical depth integration
! 20090720, BAW added this in:
          tsum=tsum+dtau
          return
       end if

       dtauabs(ir,it,ip)=dtauabs(ir,it,ip)+dt*kapd(idust2)*densarr(ir,it,ip)*(1.d0-albedo(idust2))

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
       tsum=tsum+dtau

       ! if photon exists grid, exit integration
       if (ir == nrg) then
          exitflag=1
          return
       end if

      ! if photon intersects star, re-emit
       if (ir == 0) then
          aflux=aflux+1.d0
          rtot=sqrt(rsq)
! reemit thermal photon
          call initpabs
          call opacset(nub)
          wave=1.2398d0/nub
          i_orig=1
          iflag=0
          idust2=dustarr(1,1,1)
          do id=1,4
             kapd(id)=kappa(id)*rsol*rstar*kappav(id)
          end do
c          rlam(idust2)=albedo(idust2)
          ir=1
          tsum=0.d0
          xran=ran()
          tau=-log(xran)
          !set new random optical depth!   20090720,  JB and BAW
!save photon position, xpold, ypold,zpold?????
       end if

  end subroutine tauint_cell

end module tauint_mod
