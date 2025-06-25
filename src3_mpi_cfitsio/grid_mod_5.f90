module grid_mod

  implicit none
  save

  ! *****  set ntg to an odd number!!!!! *******
  ! will work on a fix later...

  ! Grid dimensions
  integer :: nrg, ntg, npg

  real(8) :: x_max0,zdmax

  ! number of opacity species (dust/gas/sgs)
  integer, parameter :: ndg = 8
  logical,parameter :: is_any(ndg+2) =      [1,1,1,1,1,1,1,1,1,1] == 1
  logical,parameter :: is_disk(ndg+2) =     [1,1,0,0,1,1,0,0,1,0] == 1
  logical,parameter :: is_envelope(ndg+2) = [0,0,1,0,0,0,1,0,0,0] == 1
  logical,parameter :: is_cavity(ndg+2) =   [0,0,0,1,0,0,0,1,0,1] == 1
  logical,parameter :: is_sg(ndg+2) =       [0,0,0,0,1,1,1,1,0,0] == 1
  logical,parameter :: is_thermal(ndg+2) = .not. is_sg
  
  integer :: fractmask(ndg)

  real(8),parameter :: texp = 1.5d0
  real(8),parameter :: rexp = 2.0d0
  ! real(8),parameter :: rexp = 2.0   ! constant density 1-D geometry

  integer,pointer :: killed(:,:,:)

  real(8),pointer :: densarr(:,:,:,:)
  real(8),pointer :: densvararr(:,:,:)
!  real(8),pointer :: massarr(:,:,:,:)
!  real(8),pointer :: tdust(:,:,:,:)
!  real(8),pointer :: dtauabs(:,:,:,:)
  real(8),pointer :: dtauabs(:,:,:)
!  real(8),pointer :: tdust2(:,:,:,:)
  real(8),pointer :: tdave(:,:,:)
  real(8),pointer :: tdust(:,:,:)
  real(8),pointer :: jmean(:,:,:)
  real(8),pointer :: vcell(:,:,:)
  real(8),pointer :: rarr(:),r2arr(:)
  real(8),pointer :: thetarr(:),tmptharr(:)
  real(8),pointer :: sintarr(:),costarr(:)
  real(8),pointer :: tan2arr(:),phiarr(:)
  real(8),pointer :: aarr(:),barr(:),tauarr(:)
  real(8),pointer :: ravearr(:),thetavearr(:),phiavearr(:)
  real(8),pointer :: Ldint(:),hdisk(:),zdisk(:),cyndisk(:),rhod(:),midrho(:)
  real(8),pointer :: x0arr(:,:,:),vrarr(:,:,:),vthetarr(:,:,:),vphiarr(:,:,:)
  real(8),pointer :: vzarr(:,:,:),vcynarr(:,:,:)
  real(8),pointer :: divvarr(:,:,:)
  real(8),pointer :: adiarr(:,:,:),advarr(:,:,:)
  real(8),pointer :: fuvmean(:,:,:)
  real(8),pointer :: sigmadotarrf(:)
  real(8),pointer :: avarr(:,:,:),avthetarr(:,:,:)

!  integer,pointer :: nabs(:,:,:,:),diffdir(:,:,:)
  integer,pointer :: nabs(:,:,:),diffdir(:,:,:)
  integer,pointer :: imean(:,:,:),dbig(:,:,:),dsg(:,:,:)
  integer :: iwarn
  integer,pointer :: findopac(:,:,:,:)
  integer,pointer :: compo(:,:,:)
  logical,pointer :: diffus(:,:,:)

  character(len=100) :: tabname(4)

contains

  subroutine grid_allocate()

    use ttsre_mpi_mod

    implicit none

    if (myid.eq.0) write(*,'("Setting grid dimensions to (",I0,",",I0,",",I0,")")') nrg, ntg, npg

    if(nrg*ntg*npg > 1e8) stop "Number of grid cells requested exceeds 100 million"

    if(nrg==0.or.ntg==0.or.npg==0) then
       stop "ERROR: grid dimensions not set correctly"
    end if

    allocate(killed(nrg, ntg, npg))

    ! Physical quantities
    allocate(densarr(nrg,ntg,npg,ndg+2))
    allocate(densvararr(nrg,ntg,npg))
!    allocate(massarr(nrg,ntg,npg,ndg+2))
!    allocate(tdust(nrg,ntg,npg,ndg+2))
!    allocate(dtauabs(nrg,ntg,npg,ndg+2))
    allocate(dtauabs(nrg,ntg,npg))
!    allocate(tdust2(nrg,ntg,npg,ndg+2))
    allocate(tdave(nrg,ntg,npg))
    allocate(tdust(nrg,ntg,npg))
    allocate(jmean(nrg,ntg,npg))
!    allocate(nabs(nrg,ntg,npg,ndg+2),diffdir(nrg,ntg,npg))
    allocate(nabs(nrg,ntg,npg),diffdir(nrg,ntg,npg))
    allocate(imean(nrg,ntg,npg),dbig(nrg,ntg,npg),dsg(nrg,ntg,npg))
    allocate(diffus(nrg,ntg,npg))

    ! Geometrical quantities
    allocate(vcell(nrg,ntg,npg))
    allocate(rarr(nrg),r2arr(nrg))
    allocate(thetarr(ntg),tmptharr(ntg))
    allocate(sintarr(ntg),costarr(ntg))
    allocate(tan2arr(ntg),phiarr(npg))
    allocate(aarr(npg),barr(npg),tauarr(nrg))
    allocate(ravearr(nrg-1),thetavearr(ntg-1),phiavearr(npg-1))
    allocate(Ldint(nrg-1),hdisk(nrg-1),zdisk(nrg-1),cyndisk(nrg-1),rhod(nrg-1),midrho(nrg-1))
    allocate(x0arr(nrg,ntg,npg),vrarr(nrg,ntg,npg),vthetarr(nrg,ntg,npg),vphiarr(nrg,ntg,npg))
    allocate(vzarr(nrg,ntg,npg),vcynarr(nrg,ntg,npg))
    allocate(divvarr(nrg,ntg,npg))
    allocate(adiarr(nrg,ntg,npg),advarr(nrg,ntg,npg))
    allocate(fuvmean(nrg,ntg,npg))
    allocate(compo(nrg,ntg,npg))
    allocate(avarr(nrg,ntg,npg),avthetarr(nrg,ntg,npg))

    allocate(findopac(nrg,ntg,npg,ndg+2))
    
    allocate(sigmadotarrf(nrg-1))

  end subroutine grid_allocate

end module grid_mod
