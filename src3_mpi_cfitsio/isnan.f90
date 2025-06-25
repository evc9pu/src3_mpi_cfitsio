logical function isnan(a)

  implicit none

  real(8) :: a

  if (a.ne.a) then
     isnan=.true.
  else
     isnan=.false.
  end if

  return
  
end function isnan

