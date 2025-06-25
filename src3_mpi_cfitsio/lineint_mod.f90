module lineint_mod

  implicit none
  save

  integer :: nv,nspec

  real(8),pointer :: Ivthin(:,:),Ivthick(:,:),vnu(:)
  real(8),pointer :: Ivthind(:,:),Ivthine(:,:),Ivthino(:,:)
  real(8),pointer :: linemap_thin(:,:,:,:),linemap_thick(:,:,:,:)
  real(8),pointer :: linemap_thind(:,:,:,:),linemap_thine(:,:,:,:)
  real(8),pointer :: linemap_thino(:,:,:,:)

end module lineint_mod
