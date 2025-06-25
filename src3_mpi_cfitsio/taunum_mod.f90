module taunum_mod

  implicit none
  save

  real(8) :: tau,rp,xp,yp,zp,ux,uy,cosmax,sp2
  real(8) :: rfraci,zfraci,rmax2,uz,rsq,rtot

  integer :: exitflag,aflag,iflag

  logical(1) :: on_wall

end module taunum_mod
