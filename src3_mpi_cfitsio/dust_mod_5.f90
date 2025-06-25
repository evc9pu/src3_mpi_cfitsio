module dust_mod

  use grid_mod, only: ndg
  implicit none
  save

  integer,parameter :: nopac1=179
  integer,parameter :: nopac=nopac1+ndg

  real(8),dimension(nopac) :: kappa,kappav,albedo,g,pl,g2,onemg2,g2p1,twog,p1maxinv,kappav1
  real(8) :: pc,sc,Tcond

  integer,parameter :: nlammax = 800
  integer :: nlambda(nopac)
  real(8),dimension(nlammax,nopac) :: lamdust,nudust,kappad,adust,gdust,pldust

  integer,parameter :: nT=500
  integer,parameter :: nnu=513
  real(8) :: lTint(nT),Tint(nT),nuint(nnu),kappap(nT,nopac),kappar(nT,nopac)
!  real(8) :: kapint(nT,nnu,nopac),dBint(nT,nnu,nopac),kappaf(nopac)
  real(8) :: kapint(nT,nnu,nopac),kappaf(nopac)
  logical :: is_used(nopac)

  integer,parameter :: ntemper=32
  integer,parameter :: ndenmax=9
  integer,parameter :: nden(ntemper)=(/8,8,8,8,8,8,8,7,9,9,8,7,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,2,2,2,1,1/)
  real(8) :: gasden(ndenmax,ntemper)
  integer :: ntden(ntemper)

end module dust_mod
