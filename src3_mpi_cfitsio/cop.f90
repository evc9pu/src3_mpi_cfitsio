subroutine cop(rho,T,aKext)

  use dust_mod
  
  implicit none

  real(8) :: rho,T,aKext,rho1,rhothr1
  real(8) :: chiR,temp
  real(8) :: gd(ndenmax)
  integer :: idopac,tt,ii,j,jn

  rhothr1=3.34d-14
  
  if (T.lt.1400.d0) then
     if (rho.gt.rhothr1) then 
        idopac=1 
     else 
        idopac=2
     end if
  else
     
     temp=log10(T)/0.025d0

     if (temp.lt.127.5d0) then
        tt=1
     else if (temp.lt.129.d0) then
        tt=2
     else if (temp.lt.131.d0) then
        tt=3
     else if (temp.lt.133.d0) then
        tt=4
     else if (temp.lt.135.d0) then
        tt=5
     else
        temp=(temp-138.d0)/4.d0
        tt=int(temp)+6
     end if

     if (tt.gt.ntemper) then
        tt=ntemper
        print*,'gas T is too high',tt,ntemper
     endif
     
     do ii=1,ndenmax
        gd(ii)=log10(gasden(ii,tt))
     end do
     
     if (rho.gt.0.d0) then
        rho1=log10(rho)
     else
        rho1=-100.
     end if
     call locate(gd,nden(tt),rho1,j)
     if (j.eq.0) then
        jn=1
     else if (j.gt.nden(tt)-1) then
        jn=nden(tt)
!        if (rho1.gt.gd(jn)+1.d0) print*,'gas density out of boundary',rho1,tt,gd(jn)
     else if ((rho1-gd(j)).lt.(gd(j+1)-rho1)) then
        jn=j
     else
        jn=j+1
     end if
     
     idopac=ndg+ntden(tt)-nden(tt)+jn
     
  end if

  aKext=kappav(idopac)*chiR(T/11605.d0,idopac)
  
!  print*,T,rho,idopac,aKext

  return

end subroutine cop
 
