subroutine wrimsg(csubrt,cmsgnm)
  implicit none
  character(len=*),intent(in) :: cmsgnm,csubrt
  write (12,'(" >",a," ",a)') csubrt, trim(cmsgnm)
end subroutine wrimsg
