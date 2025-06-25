! ***************************************************************

subroutine zerod(sxmu,sz,xmu0,iflag)

  ! finds zeros of cubic equation
  ! xmu0**3 +(z-1)*xmu0 - z*xmu
  ! this is a special form of the equation
  ! x**3 +a1*x**2 +a2*x +a3,
  ! which we find the real roots of following
  ! numerical recipes, p. 146.
  ! programmed and checked 5/29/89 lwh

  use constants

  implicit none

  real(8) :: a1,a2,a3,q,r,q3,xmu,z,qh,denom,theta
  real(8) :: xmu01,xmu02,xmu03,factor,factor3

  real(8) :: sxmu,sz,xmu0

  integer :: iflag

  xmu=(sxmu)
  z=(sz)

  a1 = 0.d0
  a2 = (z-1.d0)
  a3 = - z*xmu
  q = (a1**2 - 3.d0*a2)/9.d0
  r = (2.d0*a1**3 - 9.d0*a1*a2 + 27.d0*a3)/54.d0
  q3 = q**3
  if ((q3 - r**2).ge.(0.d0)) then
     qh = sqrt(q)
     denom = sqrt(q3)
     theta = acos(r/denom)
     xmu01 = -2.d0*qh*cos(theta/3.d0) - a1/3.d0
     xmu02 = -2.d0*qh*cos((theta + 2.d0*pi)/3.d0) - a1/3.d0
     xmu03 = -2.d0*qh*cos((theta + 4.d0*pi)/3.d0) - a1/3.d0
     if((xmu01.gt.(0.d0)).and.(xmu02.gt.(0.d0))) then
        iflag = 1
        write(6,200) z,xmu,xmu01,xmu02
200     format(1x,' trouble: 2 roots:z,mu,mu01,mu02', 4d11.4)
        return
     else if((xmu01.gt.(0.d0)).and.(xmu03.gt.(0.d0))) then
        iflag = 1
        write(6,201) z,xmu,xmu01,xmu03
201     format(1x,' trouble: 2 roots:z,mu,mu01,mu03', 4d11.4)
        return
     else if((xmu02.gt.(0.d0)).and.(xmu03.gt.(0.d0))) then
        iflag = 1
        write(6,202) z,xmu,xmu02,xmu03
202     format(1x,' trouble: 2 roots:z,mu,mu02,mu03', 4d11.4)
        return
     end if
     if(xmu01.gt.(0.d00)) xmu0 = (xmu01)
     if(xmu02.gt.(0.d00)) xmu0 = (xmu02)
     if(xmu03.gt.(0.d00)) xmu0 = (xmu03)
     return
  else
     factor = sqrt(r**2 - q3) + abs(r)
     factor3 = factor**(1.0d0/3.0d0)
     xmu01 = -1.d0*dsign(1.d0,r)*(factor3 + q/factor3) - a1/3.d0
     xmu0 = (xmu01)
  end if
  return
end subroutine zerod




