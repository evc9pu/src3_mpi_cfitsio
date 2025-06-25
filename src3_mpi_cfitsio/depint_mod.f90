module depint_mod

  implicit none
  save

  integer :: ndep

  real(8),pointer :: S0_arr(:),Rdg_arr(:),crosssec_arr(:)
  real(8),pointer :: mu_arr(:),Eco_arr(:),abunddep_arr(:)
  real(8),pointer :: nu_arr(:),Agr_arr(:),Nb_arr(:),crrate_arr(:)
  character(len=10),pointer :: molname(:)

  real(8),pointer :: totmap(:,:,:),gasmap(:,:,:)
  real(8),pointer :: totcol(:),gascol(:)

end module depint_mod
