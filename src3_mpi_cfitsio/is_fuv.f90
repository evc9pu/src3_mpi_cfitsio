logical function is_fuv(nu)

  implicit none

  real(8) :: nu

  if (nu.le.13.6d0.and.nu.ge.6.d0) then
     is_fuv=.true.
  else
     is_fuv=.false.
  end if

  return

end function is_fuv
