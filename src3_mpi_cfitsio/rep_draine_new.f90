subroutine rep_draine(im,freq)

  ! 2007, KW. modified for my code, BAW, 20070815

  use draine_mod
  use random

  implicit none

  integer :: ifreq,im
  real(8) :: xran,logfreq,freq

  ! need to find out if PAHs are destroyed above a certain energy density

  xran=ran()

  ! call locate(cdf_d,ndraine,xran,ifreq)
  call locate3(cdf_d,nemit,ndraine,im,xran,ifreq)

  if(ifreq.eq.ndraine)ifreq=ifreq-1

  if(ifreq.eq.0) then
     print*,'ifreq=0, rep_draine; problem?'
     ifreq=1
  end if

  logfreq=(log10(xran)-logcdf_d(im,ifreq))/ &
       & (logcdf_d(im,ifreq+1)-logcdf_d(im,ifreq))* &
       & (logfreq_d(ifreq+1)-logfreq_d(ifreq)) &
       & + logfreq_d(ifreq)
  freq=10.d0**logfreq

  ! lamran=10.d0**loglam_c(ifreq)

  return
end subroutine rep_draine

