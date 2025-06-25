subroutine setdiffus(iter)

  use grid_mod
  use dust_mod
  use tts_mod
  use out_mod
  use opacin_mod
  use output_mod
  use log_mod

  implicit none

  integer :: ir,it,ip,id,nr0,itt,it1,iter,dcount
  integer :: iwfile

  real(8) :: rad,zmr,tau,tau1,lmfp,lmfp1,dt,tt
  real(8) :: chiR,kapP

  character*2 :: tfilename

  character(len=3) :: suffix

  dcount=0

  do ir=1,nrg
     do it=1,ntg
        do ip=1,npg
           diffus(ir,it,ip)=.false.
           diffdir(ir,it,ip)=0
        end do
     end do
  end do

  if (massd.gt.0.d0) then
     if (rddust.eq.rmin) then
        nr0=1
     else
        nr0=3
     end if
     ir=nr0
     rad=rarr(ir)
     do while (rad.lt.rmaxd)
        zmr=hdisk(ir)
!        print*,zmr
        do it=(ntg-1)/2,2,-1
           tau=0.d0
           itt=it
           do while (tau.lt.1.d0.and.itt.gt.1)
              dt=thetarr(itt+1)-thetarr(itt)
              tau1=tau
              tt=tdave(ir,itt,1)
              if (tt.gt.1400.d0/11605.d0.and.densarr(ir,itt,1,9)+densarr(ir,itt,1,10).eq.0.d0) then
                 print*,'error!, check temperature'
                 do id=1,ndg+2
                    if (densarr(ir,itt,1,id).gt.0.d0) print*,id,ir,itt
                 end do
                 do iwfile=1,10
                    write(suffix,'(I3.3)') iwfile
                    open(unit=24,file='WARNING_TEMPERATURE_'//suffix,status='unknown')
                    write(24,*)'error!, check temperature'
                    close(24)
                 end do
!                 stop
              end if
              if (iter.lt.2) tt=max(tt,50.d0/11605.d0)
              do id=1,ndg+2
                 if (is_thermal(id)) then
                    tau=tau+kapd(findopac(ir,itt,1,id))* &
                         & chiR(tt,findopac(ir,itt,1,id))* &
                         & densarr(ir,itt,1,id)*dt*rad
                 end if
              end do
              itt=itt-1
           end do
!           print*,itt
           if (tau.gt.1.d0) then
              lmfp1=rad*(thetarr(itt+2)-thetarr(itt+1))/(tau-tau1)*(1.d0-tau1)
           else
              lmfp1=rad*(thetarr(itt+2)-thetarr(itt+1))
           end if
           lmfp=rad*(thetarr(it+1)-thetarr(itt+2))+lmfp1
           
           if (tau.ge.1.d0.and.lmfp.lt.zmr/20.d0) then
              do it1=it,(ntg-1)/2
                 diffus(ir,it1,1)=.true.
                 diffdir(ir,it1,1)=-2
              end do
           end if
           
        end do
        ir=ir+1
        rad=rarr(ir)
     end do

     do ir=1,nrg-1
        do it=(ntg+1)/2,ntg-1
           diffus(ir,it,1)=diffus(ir,ntg-it,1)
           diffdir(ir,it,1)=-diffdir(ir,ntg-it,1)
        end do
     end do
     
     if (npg.gt.2) then
        do ir=1,nrg-1
           do it=1,ntg-1
              do ip=2,npg-1
                 diffus(ir,it,ip)=diffus(ir,it,1)
                 diffdir(ir,it,ip)=diffdir(ir,it,1)
              end do
           end do
        end do
     end if
     
  end if

  do ir=1,nrg-1
     do it=1,ntg-1
        do ip=1,npg-1
           if (diffus(ir,it,ip)) dcount=dcount+1
        end do
     end do
  end do

  if (ifprint) print*,'changed diffus grid'
  if (ifprint) print*,'number of diffusion cells in grid', dcount
  call diskdat_write('dcount',dcount,'number of diffusion cells in grid')

!  write(tfilename,'(i2)') iter
!  call output_grid('diffuse_'//trim(adjustl(tfilename)),diffus)
!  call output_grid('diffdir_'//trim(adjustl(tfilename)),diffdir)

  return
end subroutine setdiffus
