subroutine int_mean(nphot)

  ! BW & KW 2007.

  ! calculate mean intensity.

  use grid_mod
  use tts_mod
  use output_mod
  use draine_mod
  use isrf_mod
  use constants, only : pi

  implicit none

  integer :: i,j,k,nphot,im
  real(8) :: jave,count,ljrat_ave,jratio,vtot,jmeanave

  real(8),parameter :: lsol = 3.85d33

  if (ifprint) then
     call output_grid('jmeanbef',jmean)
     call output_grid('volume',vcell)
  end if
     
  if (ifprint) print*,'ltot',ltot/lsol

  jmeanave=0.d0
  do i=1,nrg-1
     do j=1,ntg-1
        do k=1,npg-1
           jmean(i,j,k)=jmean(i,j,k)/vcell(i,j,k)/dble(nphot)*ltot/pi/4.d0
           ! *lstars/pi/4.d0
           ! 20071126, use ltot instead of lstars.  I think... &

           jmeanave=jmeanave+jmean(i,j,k)
        end do
     end do
  end do
  if (ifprint) print*,'jmeanave',jmeanave

  if (ifprint) call output_grid('jmean',jmean)

  ! if (ialign.eq.1) then
  ! do i=1,nrg-1
  ! do j=1,ntg-1
  ! do k=1,npg-1
  ! do l=1,nl
  ! jml(i,j,k,l)=jml(i,j,k,l)/vcell(i,j,k)*lsun/nphot
  ! +        *ltot/pi/4.d0
  ! jmlx(i,j,k,l)=jmlx(i,j,k,l)/vcell(i,j,k)*lsun/nphot
  ! +        *ltot/pi/4.d0
  ! jmly(i,j,k,l)=jmly(i,j,k,l)/vcell(i,j,k)*lsun/nphot
  ! +        *ltot/pi/4.d0
  ! jmlz(i,j,k,l)=jmlz(i,j,k,l)/vcell(i,j,k)*lsun/nphot
  ! +        *ltot/pi/4.d0
  ! end do
  ! end do
  ! end do
  ! end do
  ! end if

  open(unit=15,file='jmean_radave.dat',status='unknown')
  write(15,*) &
       & 'i,r(rstar),r(au),  jmean (rad ave) log(jmean/0.02)'
  do i=1,nrg-1
     count=0.d0
     vtot=0.d0
     jave=0.d0
     do j=1,ntg-1
        do k=1,npg-1
           ! jratio(i,j,k)=0. ! zero out for main iteration loop
           jratio=jmean(i,j,k)/isrfkap_int ! Compare with average ISRF weighted by small grain opacity
 !          jratio=jmean(i,j,k)/2.d-2 ! Compare with average ISRF
           call locate(uratio_arr,nemit,jratio,im)
           ! if (i.eq.5)
           ! $            print*,'im,log10(jratio)',im,log10(jratio)
           if (im.gt.nemit) then
              print*,'jratio higher than 1.e7, ',log10(jratio)
              im=nemit
           end if
           if (im.lt.1) then
              im=1
           end if
           imean(i,j,k)=im
           vtot=vtot+vcell(i,j,k)
           jave=jratio*vcell(i,j,k)+jave
           count=count+1.d0
        end do
     end do
     jave=jave/count/vtot
     ljrat_ave=log10(jave)
     write(15,900) i,ravearr(i),ravearr(i)/autors &
          & ,jave,ljrat_ave
  end do
  close(15)
900 format(i10,f12.0,f10.7,f10.5,f10.5)

  return
end subroutine int_mean
