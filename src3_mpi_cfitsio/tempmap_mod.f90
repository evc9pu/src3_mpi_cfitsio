module tempmap_mod

  implicit none
  save

  real(8) :: Trhor,rhor,vrhor
  real(8),pointer :: Tmap(:,:),rhomap(:,:),vmap(:,:)

end module tempmap_mod
