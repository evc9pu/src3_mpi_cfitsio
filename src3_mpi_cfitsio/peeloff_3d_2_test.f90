module peeloff_mod
  implicit none
contains

!  subroutine peeloff(on_wall,xp1,yp1,zp1,sip,sqp,sup,svp,costp,sintp,cospp,&
!       &sinpp,phip,hit,htot,rsq1,rtot1,tsum,ir,it,ip,idust,idust2,iflag,iphot,peelid)
  subroutine peeloff(on_wall_in,xp1_in,yp1_in,zp1_in,sip_in,sqp_in,sup_in,&
       & svp_in,costp_in,sintp_in,cospp_in,sinpp_in,phip_in,hit_in,htot_in,&
       & rsq1_in,rtot1_in,tsum,ir_in,it_in,ip_in,idust_in,idust2_in,iflag_in,&
       & iphot_in,peelid_in)

    ! integrate optical depth until photon exits in a specified
    ! direction, or hits the star.

    use tts_mod
    use grid_mod
    use filt_mod
    use opacin_mod
    use messages
    use closest_wall, only : reset_t,find_next_wall
    use constants

    implicit none

    ! Input values - all copied as values instead of references
    !                (so as not to modify them in the parent program)
    logical(1) :: on_wall
    real(8)  :: xp1,yp1,zp1,sip,sqp,sup,svp
    integer :: idust,idust2
    real(8)  :: costp,sintp,cospp,sinpp,phip
    real(8)  :: rsq1,rtot1,hit,htot
    integer :: ir,it,ip
    integer :: peelid,iflag,iphot
    real(8),intent(out) :: tsum

    logical(1) :: on_wall_in
    real(8)  :: xp1_in,yp1_in,zp1_in,sip_in,sqp_in,sup_in,svp_in
    integer :: idust_in,idust2_in
    real(8)  :: costp_in,sintp_in,cospp_in,sinpp_in,phip_in
    real(8)  :: rsq1_in,rtot1_in,hit_in,htot_in
    integer :: ir_in,it_in,ip_in
    integer :: peelid_in,iflag_in,iphot_in

    ! Photon propagation
    real(8) :: ux1,uy1,uz1
    real(8) :: dtau
    integer :: iface,id
    logical :: inner_wall,outwards

    ! Constants
!    real(8),parameter :: taumax = 25.d0
    real(8),parameter :: taumax = 25.d0
    character(len=20),parameter :: routine = 'peeloff_3d'

    real(8) :: bmu,sinbm
    real(8) :: opac,phinew,ri1,factor,dt,wt,freq,rho2

    ! Image
    real(8) :: ximage,yimage
    integer :: ix,iy
    integer :: inub,iw,ift,ia,ifilt,nf,i_inter,indx(3)
    
    peeliter=peeliter+1

    on_wall=on_wall_in
    xp1=xp1_in
    yp1=yp1_in
    zp1=zp1_in
    sip=sip_in
    sqp=sqp_in
    sup=sup_in
    svp=svp_in
    idust=idust_in
    idust2=idust2_in
    costp=costp_in
    sintp=sintp_in
    cospp=cospp_in
    sinpp=sinpp_in
    phip=phip_in
    rsq1=rsq1_in
    rtot1=rtot1_in
    hit=hit_in
    htot=htot_in
    ir=ir_in
    it=it_in
    ip=ip_in
    peelid=peelid_in
    iflag=iflag_in
    iphot=iphot_in

    ! set optical depth sum to zero
    tsum=0.d0

    ! set multiple peeloff angles:
    cospe = cospe_arr(peelid)
    sinpe = sinpe_arr(peelid)
    coste = coste_arr(peelid)
    sinte = sinte_arr(peelid)

    ! set photon direction
    ux1=sinte*cospe
    uy1=sinte*sinpe
    uz1=coste

    ! need to check whether photon is going in opposite direction to cell it's in
    ! this matters because diffusion photons can start on a wall but be peeled off
    ! in direction opposite to photon travel

    if(on_wall) then
       inner_wall = abs(rsq1-r2arr(ir)) < abs(rsq1-r2arr(ir+1))
       outwards = ux1*xp1+uy1*yp1+uz1*zp1 > 0.d0
       if(inner_wall.and..not.outwards) then
          if(ir==1) then
             return
          else
             ir = ir - 1
          end if
       end if
       if(.not.inner_wall.and.outwards) then
          ir = ir + 1
       end if
    end if

    peeliter1=peeliter1+1

    do
       
       peeliter2=peeliter2+1

       call cpu_time(peeloff_start)

    !       call cpu_time(peeloff1_start)
       call testgrid(xp1,yp1,zp1,rsq1,rtot1,ir,it,ip,on_wall,routine)
!       call cpu_time(peeloff1_stop)
!       peeloff1_time=peeloff1_time+peeloff1_stop-peeloff1_start

       call reset_t()

       ! find shortest non-negative distance to two constant radius surfaces
!       call cpu_time(peeloff2_start)
       call radface(ux1,uy1,uz1,xp1,yp1,zp1,r2arr(ir),r2arr(ir+1))
!       call cpu_time(peeloff2_stop)
!       peeloff2_time=peeloff2_time+peeloff2_stop-peeloff2_start

       ! find shortest positive distance to two constant theta surfaces
!       call cpu_time(peeloff3_start)
       if (ntg.gt.2) call thetaface(ux1,uy1,uz1,xp1,yp1,zp1,tan2arr(it),tan2arr(it+1))
!       call cpu_time(peeloff3_stop)
1       peeloff3_time=peeloff3_time+peeloff3_stop-peeloff3_start

       ! find shortest non-negative distance to two constant phi surfaces
       if (npg.gt.2.and.ntg.gt.2) call phiface(aarr(ip),barr(ip),0.d0,0.d0,aarr(ip+1),barr(ip+1),0.d0,0.d0,ux1,uy1,uz1,xp1,yp1,zp1)

!       call cpu_time(peeloff4_start)
       call find_next_wall(on_wall,dt,iface)
!       call cpu_time(peeloff4_stop)
!       peeloff4_time=peeloff4_time+peeloff4_stop-peeloff4_start

    
       call cpu_time(peeloff_stop)
       peeloff_time=peeloff_time+peeloff_stop-peeloff_start
       
       opac=0.d0
       do id=1,ndg+2
          if (densarr(ir,it,ip,id).gt.0.d0) then
             opac=opac+kapd(findopac(ir,it,ip,id))*densarr(ir,it,ip,id)
          end if
       end do

       dtau=opac*dt

       tsum=tsum+dtau
       
       ! if photon weight is too small to matter, stop peeloff
       if (tsum.gt.taumax) then
          return
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

       xp1=xp1+dt*ux1
       yp1=yp1+dt*uy1
       zp1=zp1+dt*uz1
       on_wall = .true.
       rsq1=xp1*xp1+yp1*yp1+zp1*zp1
       rtot1=sqrt(rsq1)

       ! if photon exists grid, exit integration
       if (ir == nrg) exit

       ! if photon intersects star, stop peeloff
       if (ir == 0) then
          return
       end if

    end do
    
    peeliter3=peeliter3+1

    ! if you are here, you left the atmosphere and are heading to HST
    ! weight photon intensity by optical depth and
    ! account for probability of intercepting solid angle of scattered light.
 
    factor=exp(-tsum)
    sip=sip*factor
    sqp=sqp*factor
    sup=sup*factor
    svp=svp*factor

    if (iflag.eq.1) then

       call findangle3(costp,sintp,phip,cospp,sinpp,coste,sinte,phie,cospe,sinpe,bmu,sinbm,phinew,ri1)
       call stokespeel(idust,sip,sqp,sup,svp,coste,sinte,hit,htot,ri1,bmu,costp,sintp,findopac(ir,it,ip,idust2))

    end if

    yimage = zp1*sinte-yp1*coste*sinpe-xp1*coste*cospe
    ximage = yp1*cospe-xp1*sinpe

    ix = int((ximage+rmaxi)/fractxh)+1
    iy = int((yimage+rmaxi)/fractxh)+1

    ! change wave to frequency
    freq=2.9979d14/wave

    if(output_imfilt) then

       if (ix.le.nxhst.and.iy.le.nxhst.and.ix.gt.0.and.iy.gt.0) then

          ! see if wavelength falls into any of the filter functions
          do ifilt=1,nfilt
             ! nf=number of lines in this filter file
             nf=nwfilt(ifilt)-1
             if (freq.gt.filtwave(1,ifilt).and.freq.lt.filtwave(nf,ifilt))   then
                ! just be really stupid and loop through array!
                ! they aren't that big...
                iw=1
                do while (freq.gt.filtwave(iw,ifilt))
                   iw=iw+1
                end do
                iw=iw-1
                ift=ifilt
                if (iw.lt.1.or.iw.ge.nf) print*,'iw wrong',iw,ift
                wt=filtphi(iw,ifilt)

                image_b_i(peelid,ix,iy,ift) = image_b_i(peelid,ix,iy,ift) + real(sip*wt)
                image_b_q(peelid,ix,iy,ift) = image_b_q(peelid,ix,iy,ift) + real(sqp*wt)
                image_b_u(peelid,ix,iy,ift) = image_b_u(peelid,ix,iy,ift) + real(sup*wt)
                image_b_v(peelid,ix,iy,ift) = image_b_v(peelid,ix,iy,ift) + real(svp*wt)

                if (i_orig.eq.1) then
                   image_b_s(peelid,ix,iy,ift) = image_b_s(peelid,ix,iy,ift) + real(sip*wt)
                else if (i_orig.eq.2) then
                   image_b_d(peelid,ix,iy,ift) = image_b_d(peelid,ix,iy,ift) + real(sip*wt)
                else if (i_orig.eq.3) then
                   image_b_e(peelid,ix,iy,ift) = image_b_e(peelid,ix,iy,ift) + real(sip*wt)
                else if (i_orig.eq.4) then
                   image_b_o(peelid,ix,iy,ift) = image_b_o(peelid,ix,iy,ift) + real(sip*wt)
                end if

             end if
          end do
       end if

    end if

    rho2=ximage**2+yimage**2

    ! bin in frequency
    inub=min(max(int(dble(nfreq)*log(nub/numin)/lnurat)+1,1),nfreq)

    ! Monochromatic image

    if(output_imcube) then
       if (ix.le.nxhst.and.iy.le.nxhst.and.ix.gt.0.and.iy.gt.0) then
          image_m_i(peelid,ix,iy,inub) = image_m_i(peelid,ix,iy,inub) + real(sip)
          image_m_q(peelid,ix,iy,inub) = image_m_q(peelid,ix,iy,inub) + real(sqp)
          image_m_u(peelid,ix,iy,inub) = image_m_u(peelid,ix,iy,inub) + real(sup)
          image_m_v(peelid,ix,iy,inub) = image_m_v(peelid,ix,iy,inub) + real(svp)
       end if
    end if

    ! Sparse image

    if(output_imsparse) then
       call nested_image_bin(image_ms_i(peelid),real(ximage),real(yimage),real(log10(freq)),real(sip))
       call nested_image_bin(image_ms_q(peelid),real(ximage),real(yimage),real(log10(freq)),real(sqp))
       call nested_image_bin(image_ms_u(peelid),real(ximage),real(yimage),real(log10(freq)),real(sup))
       call nested_image_bin(image_ms_v(peelid),real(ximage),real(yimage),real(log10(freq)),real(svp))
    end if

    if(iflag==1) then
       i_inter=2
    else
       if(i_orig==1) then
          i_inter=1
       else
          i_inter=3
       end if
    end if

    indx = (/1,1+i_orig,5+i_inter/)

    do ia=1,nap
       if(rho2 < aperture2(ia)) then

          ! add photon to flux arrays
          ti(inub,peelid,ia,indx) = ti(inub,peelid,ia,indx) + real(sip)
          tq(inub,peelid,ia,indx) = tq(inub,peelid,ia,indx) + real(sqp)
          tu(inub,peelid,ia,indx) = tu(inub,peelid,ia,indx) + real(sup)
          tv(inub,peelid,ia,indx) = tv(inub,peelid,ia,indx) + real(svp)

          ! add photon to flux error arrays
          ti2(inub,peelid,ia,indx) = ti2(inub,peelid,ia,indx) + real(sip*sip)
          tq2(inub,peelid,ia,indx) = tq2(inub,peelid,ia,indx) + real(sqp*sqp)
          tu2(inub,peelid,ia,indx) = tu2(inub,peelid,ia,indx) + real(sup*sup)
          tv2(inub,peelid,ia,indx) = tv2(inub,peelid,ia,indx) + real(svp*svp)

          numt(inub,peelid,ia,indx) = numt(inub,peelid,ia,indx) + 1.

          peeliter4=peeliter4+1

       end if
    end do

    return
    
  end subroutine peeloff

end module peeloff_mod
