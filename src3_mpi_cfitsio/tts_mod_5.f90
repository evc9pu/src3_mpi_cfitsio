! make larger pixel array

module tts_mod

  use type_nested_image
  use grid_mod, only: ndg

  implicit none
  save

  character(len=10),parameter :: version = '20090821'

  ! STANDARD SEDS
  ! si/sq/su/sv     = used for fluxes
  ! si2/sq2/su2/sv2 = used for uncertainties
  ! nums            = number of photons used
  ! aveinc          = ???

  integer :: nfreq ! number of frequencies for SEDs
  ! 3-D outgoing angles
  ! integer,parameter :: nmu    = 20  ! number of mu viewing angles for SEDs
  ! integer,parameter :: nph    = 40  ! number of phi viewing angles for SEDs
  ! 2-D outgoing angles
  integer :: nmu ! number of mu viewing angles for SEDs
  integer :: nph ! number of phi viewing angles for SEDs
  ! 1=D outgoing angles
  ! integer,parameter :: nmu    = 2
  ! integer,parameter :: nph    = 1  ! number of phi viewing angles for SEDs
  !
  integer,parameter :: no     = 8   ! number of separate components
  ! (total/stellar/disk/envelope/outside/direct/scattered/thermal)

  logical :: diffusion = .true.
  logical :: partial_peeloff = .false.

  ! (nfreq,nmu,nph,nap,no)
  real,pointer,dimension(:,:,:,:,:) :: si,sq,su,sv,si2,sq2,su2,sv2,nums,aveinc

  ! PEELOFF ANGLES

  integer :: npeel ! number of peel-off angles

  real(8),pointer,dimension(:) :: thete_arr,phie_arr
  real(8),pointer,dimension(:) :: coste_arr,sinte_arr
  real(8),pointer,dimension(:) :: cospe_arr,sinpe_arr

  ! PEELOFF SEDS
  ! ti/tq/tu/tv     = used for fluxes
  ! ti2/tq2/tu2/tv2 = used for uncertainties
  ! numt            = number of photons used
  ! aveinc          = ???
  ! nfreq,npeel,nap,no

  real,pointer,dimension(:,:,:,:) :: ti,tq,tu,tv,ti2,tq2,tu2,tv2,numt

  ! PEELOFF IMAGES
  ! tihst/tqhst/tuhst/tvhst = used for fluxes

  ! nbnd needs to be equal to that in filt.txt
  integer,parameter :: nx    = 1
  integer,parameter :: nxhst = 601 ! image size
  integer,parameter :: nbnd  = 49  ! number of filters

  ! Broadband:
  ! (npeel,nxhst,nxhst,nbnd)

  real,pointer,dimension(:,:,:,:) :: image_b_i,image_b_q,image_b_u,image_b_v
  real,pointer,dimension(:,:,:,:) :: image_b_d,image_b_s,image_b_e,image_b_o

  ! Monochromatic:
  ! (npeel,nxhst,nxhst,nfreq)

  real,pointer,dimension(:,:,:,:) :: image_m_i,image_m_q,image_m_u,image_m_v

  ! Sparse

  type(nested_image),pointer,dimension(:) :: image_ms_i,image_ms_q,image_ms_u,image_ms_v

  real(8),pointer :: tauenv(:,:,:)

  real(8),pointer :: aperture2(:)
  real(8),pointer :: u(:)
  real(8) :: b(ndg),a(ndg),rho0(ndg),fmass(ndg),taur(ndg)
  real(8) :: z1(ndg)
  real(8) :: nscat,tot,zmax,rmax,lp,abflux,rcore
  real(8) :: rmin,rmind,rddust,flux,xmaxp,rmind2
  real(8) :: sflux,aflux,cosb,sinb,rmaxd,rmaxi
  real(8) :: rpolbeam,dmu,dmuhalf,thetmu0,autors,rstar,wave
  real(8) :: tstar,tstarave,normstar,thete,phie,coste,sinte,cospe,sinpe
  real(8) :: sfactor,fractxh,numax,numin,lnurat,nunorm
  real(8) :: nurat,nub,accfrac,accfrac2,alphad
  real(8) :: tsub,tshock,rtrunc,fspot,mdotdisk,apmax,apmin,ltot
  real(8) :: dthresh,l_isrf,isrf_Av,isrf_scl,tdiffavesave
  real(8) :: sigsave,isrf_frac,rmin_sg
  real(8) :: sigsave1,sigsave2,tdiffavesave1,tdiffavesave2
  real(8) :: rhodens1,envexp
  real(8) :: rgapd1,rgapd2,fractd,rhogapd,rgape1,rgape2,fracte,rhogape
  real(8) :: pitch,sw
  real(8) :: densratio,fractl,fw,fwesc
  real(8) :: ifuv,ifuv_ori,radmax
  real(8) :: alphamin,alphamax,alphafinal
  real(8) :: windext,clumpden,clumppow,sigmacl,betac
  real(8),target :: scount,dscount
  real :: peeloff_start,peeloff_stop,peeloff_time
  real :: peeloff1_start,peeloff1_stop,peeloff1_time
  real :: peeloff2_start,peeloff2_stop,peeloff2_time
  real :: peeloff3_start,peeloff3_stop,peeloff3_time
  real :: peeloff4_start,peeloff4_stop,peeloff4_time
 
  integer,target :: np
  integer :: npmin,n_iter_min,n_iter_max,npout,npbkg,iout
  integer :: npsav,limb,occult,nfinei,nfined
  integer :: nri,nzi,iwrite,isot,itherm,ipeel,ispot,i_orig,idiskacc
  integer :: ialpha,irminsub,izmin,iplanckst,nap,idiskwarp,ilucy
  integer :: ispiral,nspiral,ienvtype,ifractal
  integer :: igapd,igapddens,igape,igapedens,sexp
  integer :: np_fuv,nps_fuv,npd_fuv,nphs_fuv
  integer :: iterstart
  integer :: peeliter,peeliter1,peeliter2,peeliter3,peeliter4
  
  logical :: output_imfilt,output_imcube,output_imsparse
  
  integer :: sparse_dim, sparse_factor
  real :: sparse_rmin, sparse_rmax

  integer :: imave = 0
  integer :: imct = 0
  
  character(len=10) :: czmin,ifclump

  character(len=50) :: par_dir
  character(len=50) :: histdir

  logical :: ifprint,iffuv

end module tts_mod

