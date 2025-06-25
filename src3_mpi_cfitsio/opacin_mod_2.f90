module opacin_mod

  use dust_mod, only: nopac
  use grid_mod, only: ndg
  implicit none

  save

  real(8) :: kapd(nopac),rate
  real(8) :: opb,iopb,massc,rc,sonicp,rchole,rmine,ex1,z01,rhoconst1
  real(8) :: exf,thet1,ex2,z02,rhoconst2,thet2,zflowmin,zbub1,zbub2,nbub
  real(8) :: buboa,windmu0,rhoe0,c1e,c2e,cosbuboa,rhoamb,rminsq
  real(8) :: xymaxdens,rmine_in,rmind_in,mcore

  integer :: ihole,ipoly,istream,ibub,ifseed

end module opacin_mod
