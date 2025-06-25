subroutine reapar(filename,atname,dustname,i1)

  use tts_mod
  use opacin_mod
  use out_mod
  use filt_mod, only : filter_dir
  use spot_mod
  use messages, only : set_debug, set_warnings, set_errors
  use random
  use grid_mod, only : nrg, ntg, npg, is_disk, is_envelope, is_cavity, fractmask, tabname
  use constants

  use configuration

  implicit none

  character(len=*),intent(in) :: filename
  character(len=*),intent(out) :: atname,dustname(ndg)

  character(len=50) :: cshape,cgapd,cgape,cspiral,cgapedens,cgapddens,cfractal,envtype

  integer :: iwarn

  real(8) :: fmassd1,fsg(4)
  real(8) :: rclump
  integer :: i1,id,fractmasktmp(4),itab

  call load_parameter_file(filename)
  
!  call get_parameter_value('PHIE',phie_arr,'Phi angle(s) (deg) of high S/N image(s)/SED(s)')
!  call get_parameter_value('THETE',thete_arr,'Theta angle(s) (deg) of high S/N image(s)/SED(s)')
 
  ! Preliminaries

  call get_parameter_value('NP',np,'Total number of photons')
  call get_parameter_value('NPMIN',npmin,'Minimum number of photons for Lucy temperature calculation')
  call get_parameter_value('NIMIN',n_iter_min,'Minimum number of iterations for Lucy temperature calculation')
  call get_parameter_value('NIMAX',n_iter_max,'Maxmimum number of iterations for Lucy temperature calculation')
  call get_parameter_value('NISTART',iterstart,'the iteration continued from last run')
  call get_parameter_value('CPEEL',ipeel,'Use peel-off algorithm for SEDs and images')
  call get_parameter_value('CLUCY',ilucy,'Use the Lucy method for temperature calculation')
  call get_parameter_value('I1',i1,'Random number seed')
  i1 = -abs(i1)
!  call set_random_seed(i1)
  call get_parameter_value('CWARN',iwarn,'Whether to print out all warnings and errors')
  call set_warnings(iwarn)
  call set_errors(iwarn)
  call set_debug(iwarn)
  call get_parameter_value('IWRITE',iwrite,'How often to print statistics')

  npout=0

  call get_parameter_value('PARDIR',par_dir,'Directory containing parameter files')
  filter_dir = par_dir
  call get_parameter_value('HISTDIR',histdir,'Directory containing history files')

  ! External illumination
  call get_parameter_value('COUTILLUM',iout,'Use outside illumination')
  call get_parameter_value('ISRF_SCL',isrf_scl,'Scale factor for the ISRF')
  call get_parameter_value('ISRF_AV',isrf_av,'Extinction of the ISRF')

  ! SGs
  call get_parameter_value('RMIN_SG',rmin_sg,'minimum radius of SGs/vsgs')

  ! Central source properties
  call get_parameter_value('CLIMB',limb,'Limb-darkening for central source')
  call get_parameter_value('CPLANCKST',iplanckst,'Planck function for stellar emission')
  call get_parameter_value('ATNAME',atname,'atmosphere file')
  atname = trim(par_dir)//'/'//trim(atname)
  call get_parameter_value('RSTAR',rstar,'Stellar radius in Solar radii')
  call get_parameter_value('TSTAR',tstar,'Blackbody temperature of central star (K)')
  call get_parameter_value('MASSC',massc,' Mass of central star (for TSC properties)')
  call get_parameter_value('CSPOT',ispot,'emission from star+spot')
  spotflag = ispot
  call get_parameter_value('SPOTLAT',spotlat,'latitude of spot (degrees)')
  spotlon=0.d0
  call get_parameter_value('NSPOT',nspot,'number of spots')
  !  if(nspot>2) stop "ERROR: nspot cannot be larger than 2"

  ! Disk properties
  call get_parameter_value('MASSD',massd,'Disk mass in solar masses')
  call get_parameter_value('RMAXD',rmaxd,'Disk outer radius in AU')
  call get_parameter_value('CRMINSUB',irminsub,'Scale inner radius to Rsub (instead of Rstar)')
  call get_parameter_value('RMIND',rmind_in,'Disk outer radius in Rstar or Rsub')

  rmind=rmind_in
  rddust=rmind_in

  call get_parameter_value('CZMIN',czmin,'Calculate disk scale height based on HSEQ at Rsub')

  select case(trim(czmin))
  case('RSTAR')
     if (ifprint) write(*,'("ZSCALE is the scaleheight in R_star at R_star")')
  case('RSUB')
     if (ifprint) write(*,'("ZSCALE is the fraction of HSEQ scaleheight at Rsub")')
  case('R100')
     if (ifprint) write(*,'("ZSCALE is the scaleheight in AU at 100AU")')
  case default
     stop "CZMIN should be one of RSTAR/RSUB/R100"
  end select

  call get_parameter_value('FMASSD1',fmassd1,' fraction of mass in disk 1 (Large grains settled disk)')

  fmass = -huge(1.d0)
  z1 = -huge(1.d0)
  a = -huge(1.d0)
  b = -huge(1.d0)

  ! gap in disk.  
  call get_parameter_value('CGAPD',cgapd,'gap in disk?')
  call get_parameter_value('RGAPD1',rgapd1,'Inner gap radius (AU)')
  call get_parameter_value('RGAPD2',rgapd2,'Outer gap radius (AU)')
  call get_parameter_value('CGAPDDENS',cgapddens,' scale density in gap (SCALE) or constant density (CONST')
  call get_parameter_value('FRACTD',fractd,'gap density scale factor (if CGAPDENS=SCALE)')
  call get_parameter_value('RHOGAPD',rhogapd,'density in disk gap if CGAPDENS=CONST')

 ! spiral arms in disk
  call get_parameter_value('CSPIRAL',cspiral,'spiral warps in disk (YES or NO)')
  call get_parameter_value('PITCH',pitch,'pitch angle of spiral arms (degrees)')
  call get_parameter_value('SEXP',sexp,'parameter that determines width of the arms')
  call get_parameter_value('SW',sw,'fraction of dust mass entrained in the arms')
 
  ! Disk warp

  call get_parameter_value('CDISKWARP',idiskwarp,'Include disk warp')
  call get_parameter_value('WSC',wsc,'Scale for disk warp (1 = double scaleheight)')
  call get_parameter_value('WEXP',wexp,' Exponent for disk warp (cos**wexp)')


  ! Large grains settled disk
  call get_parameter_value('DUSTNAME(1)',dustname(1),'dust file, thermal grains [disk 1]')
  call get_parameter_value('DUSTNAME(5)',dustname(5),'dust file, <200A grains [disk 1]')
  call get_parameter_value('FSG(1)',fsg(1),'fraction of the mass in <200A grains')
  call get_parameter_value('ZSCALE(1)',z1(1),'Scale height of disk at Rstar')
  call get_parameter_value('A(1)',a(1),'Disk density exponent (~ r^(-a))')
  call get_parameter_value('B(1)',b(1),'Disk scale height exponent (~ r^(b))')

  ! Normal disk
  call get_parameter_value('DUSTNAME(2)',dustname(2),'dust file, thermal grains [disk 2]')
  call get_parameter_value('DUSTNAME(6)',dustname(6),'dust file, <200A grains [disk 2]')
  call get_parameter_value('FSG(2)',fsg(2),'fraction of the mass in <200A grains')
  call get_parameter_value('ZSCALE(2)',z1(2),'Scale height of disk at Rstar')
  call get_parameter_value('A(2)',a(2),'Disk density exponent (~ r^(-a))')
  call get_parameter_value('B(2)',b(2),'Disk scale height exponent (~ r^(b))')

  ! SG disks
  z1(5:6) = z1(1:2)
  a(5:6) = a(1:2)
  b(5:6) = b(1:2)

  fmass(1) = fmassd1 * (1. - fsg(1))
  fmass(5) = fmassd1 * fsg(1)
  fmass(2) = (1. - fmassd1) * (1. - fsg(2))
  fmass(6) = (1. - fmassd1) * fsg(2)

  ! Disk accretion

  call get_parameter_value('CDISKACC',idiskacc,'Include disk accretion luminosity')
  if(idiskacc==1) then
     call get_parameter_value('CALPHA',ialpha,'using alpha parameter to calculate disk Mdot')
     if (ialpha.eq.1) then
        call get_parameter_value('ALPHA',alphad,'Disk alpha parameter')
     else
        call get_parameter_value('ALPHA',mdotdisk,'Disk accretion rate')
     end if
  end if

  call get_parameter_value('ALPHAMIN',alphamin,'lower boundary of alpha')
  call get_parameter_value('ALPHAMAX',alphamax,'upper boundary of alpha')

  call get_parameter_value('RTRUNC',rtrunc,'Magnetosphere co-rotation radius')
  call get_parameter_value('FSPOT',fspot,'Fractional area of spot')
!  print*,'ispot',ispot
  if (ispot.eq.0) then
    fspot=1.d0
    if (ifprint) print*,'no stellar hotspots, so fspot is reset to 1'
  endif

  if (fspot.eq.0.d0.and.idiskacc.eq.1) then
     print*,'ERROR: fspot=0 and CDISKACC=yes - 0 area hotspot is a singularity'
     print*,'       maybe you meant fspot=1? '
     stop
  end if

  ! reset some variables if disk mass is 0
  if (massd.eq.0.d0) then
     czmin='RSUB'
     z1=1.d0
     idiskacc=0
     idiskwarp=0
     igapd=0
     ispiral=0
     if (ifprint) print*,'since disk mass is zero, setting:'
     if (ifprint) print*,'CZMIN=RSUB, zscale=1.0, CDISKACC=NO '
  end if

  if (cgapd.eq.'YES') then
    igapd=1
    if (cgapddens.eq.'SCALE') then
      igapddens=1
    else if (cgapddens.eq.'CONST') then
      igapddens=0
    else
	  print*,'ERROR:  CGAPdDENS should be set to SCALE or CONST; stopping program'
	  stop
    endif
  else
    igapd=0
    igapddens=0
    fractd=1.d0  
  endif

  if (cspiral.eq.'YES')  then
    ispiral=1
  else
    ispiral=0
  endif

  if (ispiral.eq.1.and.idiskwarp.eq.1) then
    print*,'WARNING:  you have both spiral warps and an inner warp.  do you want that?'
  endif

  ! Envelope properties

  call get_parameter_value('DUSTNAME(3)',dustname(3),'dust file, thermal grains [envelope]')
  call get_parameter_value('FSG(3)',fsg(3),'fraction of the mass in <200A grains [envelope]')
  call get_parameter_value('DUSTNAME(7)',dustname(7),'dust file, <200A grains [envelope]')

  fmass(3) = 1. - fsg(3)
  fmass(7) = fsg(3)

  call get_parameter_value('RADTAB',tabname(1),'radii of the input density table for the envelope')
  call get_parameter_value('THETATAB',tabname(2),'theta of the input density table for the envelope')
  call get_parameter_value('DENTAB',tabname(3),'density of the input density table for the envelope')
  call get_parameter_value('VTAB',tabname(4),'velocity of the input density table for the envelope')
  

  do itab=1, 4
     tabname(itab) = trim(par_dir)//'/'//trim(tabname(itab))
  end do

  call get_parameter_value('RCORE',rcore,'Core outer radius')
  call get_parameter_value('RMINE',rmine_in,'Envelope inner radius in Rstar or Rsub')
  call get_parameter_value('MCORE',mcore,'initial mass of the core')
  call get_parameter_value('SIGMACL',sigmacl,'mean surface density of the clump')
  call get_parameter_value('BETAC',betac,'rotation-to-gravition energy ratio')
  ! gap in envelope.  


  rmine=rmine_in

  call get_parameter_value('ENVTYPE',envtype,'ULRICH or POWLAW')
  call get_parameter_value('RATE',rate,'Mass infall rate for ULRICH env (M_O/yr)')
  call get_parameter_value('RC',rc,'centrifugal radius for ULRICH env')
  call get_parameter_value('SONICP',sonicp,'sonic point for ULRICH env')
  call get_parameter_value('RHODENS1',rhodens1,'fiducial density at 1 AU for POWLAW env')
  call get_parameter_value('ENVEXP',envexp,'exponent for POWLAW env')

  if (ENVTYPE.eq.'ULRICH')  then
    ienvtype=1
  else if (ENVTYPE.eq.'POWLAW') then
    ienvtype=0
  else
    print*,'ERROR:  ENVTYPE should be set to ULRICH or POWLAW; stopping program'
    stop
  endif

  ! gap in envelope.  
  call get_parameter_value('CGAPE',cgape,'gap in envelope?')
  call get_parameter_value('RGAPE1',rgape1,'Inner gap radius (AU)')
  call get_parameter_value('RGAPE2',rgape2,'Outer gap radius (AU)')
  call get_parameter_value('CGAPEDENS',cgapedens,' scale density in gap (SCALE) or constant density (CONST)')
  call get_parameter_value('FRACTE',fracte,'gap density scale factor (if CGAPDENS=SCALE)')
  call get_parameter_value('RHOGAPE',rhogape,'density in envelope gap if CGAPDENS=CONST')

  if (cgape.eq.'YES') then
    igape=1
    if (cgapedens.eq.'SCALE') then
      igapedens=1
    else if (cgapedens.eq.'CONST') then
      igapedens=0
    else
	  print*,'ERROR:  CGAPEDENS should be set to SCALE or CONST; stopping program'
	  stop
    endif
  else
    igape=0
    igapedens=0
    fracte=1.d0  
  endif

  ! Fractal density variations
  call get_parameter_value('CFRACTAL',cfractal,'YES for fractal density variations, or NO')
  call get_parameter_value('DENSRATIO',densratio,'ratio of background to fractal density')
  call get_parameter_value('FRACTMASK',fractmasktmp,'apply to [disk1,disk2,enve,outflow], 0=no 1=yes')
  call get_parameter_value('IFSEED',ifseed,'random number seed, negative integer')
!  call get_parameter_value('FRACTL',fractl,'1/length scale.') 
  fractl=3.792
  if (cfractal.eq.'YES') then
    ifractal=1
  else
	ifractal=0
  endif

  fractmask(1)=fractmasktmp(1)
  fractmask(2)=fractmasktmp(2)
  fractmask(3)=fractmasktmp(3)
  fractmask(4)=fractmasktmp(4)
  fractmask(5)=fractmasktmp(1)
  fractmask(6)=fractmasktmp(2)
  fractmask(7)=fractmasktmp(3)
  fractmask(8)=fractmasktmp(4)
  
  if (ifprint) print*,'fractmask',fractmask

  ! Cavity properties

  call get_parameter_value('CHOLE',ihole,'Bipolar cavity')
  call get_parameter_value('DUSTNAME(4)',dustname(4),'dust file, thermal grains [cavity]')
  call get_parameter_value('FSG(4)',fsg(4),'fraction of the mass in <200A grains [cavity]')
  call get_parameter_value('DUSTNAME(8)',dustname(8),'dust file, <200A grains [cavity]')

  fmass(4) = 1. - fsg(4)
  fmass(8) = fsg(4)

  if (ihole==1) then

     call get_parameter_value('FW',fw,'wind massloss rate/stellar accretion rate')
     call get_parameter_value('FWESC',fwesc,'escape fraction in the wind')

     call get_parameter_value('CSHAPE',cshape,'')
     if (cshape.eq.'STREAM') then
        ipoly=0
        istream=1
        call get_parameter_value('RCHOLE',rchole,'Streamline hole size in AU')
!        call get_parameter_value('THETMU0',thetmu0,'Opening angle of streamline wall (deg)')
     else
        ipoly=1
        istream=0
        call get_parameter_value('EX1',ex1,'Cavity wall exponent')
     end if

     if (cshape.eq.'STREAM'.and.ENVTYPE.eq.'POWLAW') then
        print*,'ERROR: CSHAPE=STREAM and ENVTYPE=POWLAW not allowed'
        print*,'For ENVTYPE=POWLAW, set CSHAPE=POLYN'
        print*,'stopping program'
        stop
     endif

     call get_parameter_value('THET1',thet1,'Opening angle of cavity wall')
     call get_parameter_value('Z01',z01,'Height of wall at w=0')

     call get_parameter_value('EXF',exf,'exponent for cavity density power-law')
     call get_parameter_value('RHOCONST1',rhoconst1,'Coefficient for cavity density')

     thetmu0=thet1
     thet2=thet1
     ex2=ex1
     z02=z01
     rhoconst2=rhoconst1

  end if

   call get_parameter_value('RHOAMB',rhoamb,'Ambient density')
   call get_parameter_value('IFCLUMP',ifclump,'IF include clump')
   call get_parameter_value('WINDEXT',windext,'Wind extension')
   call get_parameter_value('CLUMPDEN',clumpden,'density ratio of clump to core at core boundary')
   call get_parameter_value('CLUMPPOW',clumppow,'power-law index of clump density profile')

   if (ifprint) then
      print*,''
      print*,''
      print*,''
      print*,'======================'
   end if
   if (IFCLUMP.eq.'NO'.or.IFCLUMP.eq.'no') then
      RMAX=RCORE*WINDEXT
   else
      clumpden=clumpden*3.d0*mcore*msol/8.d0/pi/(rcore*au2cm)**3
      if (ifprint) print*,'clumpden',clumpden
      if (clumppow.eq.1.d0) then 
         rclump=exp(sigmacl/2.d0/rcore/au2cm/clumpden)
      else
         rclump=(1.d0+(sigmacl/2.d0/rcore/au2cm/clumpden)*(1.d0-clumppow))**(1.d0/(1.d0-clumppow))
      end if
      if (ifprint) print*,'rclump',rclump
      rmax=rclump*rcore
   end if
   if (ifprint) print*,'rmax/rcore',rmax/rcore
   if (ifprint) then
      print*,'======================='
      print*,''
      print*,''
      print*,''
   end if

  nbub=4.d0
  zbub1=250.d0
  zbub2=300.d0
  buboa=0.0d0

  ! Output parameters

  call get_parameter_value('IMFILT',output_imfilt,'Outputing images convolved with filters')
  call get_parameter_value('IMCUBE',output_imcube,'Outputing multi-wavelength data cubes')
  call get_parameter_value('IMSPARSE',output_imsparse,'Outputing multi-wavelength sparse matrice')

  call get_parameter_value('NPEEL',npeel,'Number of peel-off angles for output')
  call get_parameter_value('NAP',nap,'Number of apertures for output')

  allocate(thete_arr(npeel),coste_arr(npeel),sinte_arr(npeel))
  allocate(phie_arr(npeel),cospe_arr(npeel),sinpe_arr(npeel))

  allocate(aperture2(nap))

  call get_parameter_value('RMAXI',rmaxi,'Image half-size in AU')
  call get_parameter_value('APMIN',apmin,'radius of smallest aperture in AU')
  call get_parameter_value('APMAX',apmax,'radius of largest aperture in AU')

  rmaxi=rmax*1.1
  
  call get_parameter_value('THETE',thete_arr,'Theta angle(s) (deg) of high S/N image(s)/SED(s)')
  call get_parameter_value('PHIE',phie_arr,'Phi angle(s) (deg) of high S/N image(s)/SED(s)')

  phie=0.d0

  if(any(thete_arr == 0.)) stop "ERROR: THETE has to be <> 0"
    
  ! Advanced parameters

  call get_parameter_value('PARTIALPEEL',partial_peeloff,"Do partial peeling-off (only accretion and scattered photons)")
  call get_parameter_value('DIFFUSION',diffusion,"Whether to use the diffusion approximation in very optically thick regions")

  call get_parameter_value('SPARSEF',sparse_factor,"Sparse image factor")
  call get_parameter_value('SPARSED',sparse_dim,"Sparse image size")
  call get_parameter_value('SPARSERMIN',sparse_rmin,"Sparse image inner image size")
  call get_parameter_value('SPARSERMAX',sparse_rmax,"Sparse image outer image size")

  call get_parameter_value('FUV',iffuv,"Whether to calculate the fuv field")

  if (iffuv) then
     diffusion=.false.
     ipeel=0
  end if

  ! GRID

  call get_parameter_value('NRG',nrg,'Number of radial cells')
  call get_parameter_value('NTG',ntg,'Number of theta cells')
  call get_parameter_value('NPG',npg,'Number of phi cells')
  
  ! IMAGES/SEDS
  
  call get_parameter_value('NMU',nmu,'Number of theta bins')
  call get_parameter_value('NPH',nph,'Number of phi bins')
  call get_parameter_value('NFREQ',nfreq,'Number of frequencies')

  do id=1,ndg
     dustname(id) = trim(par_dir)//'/'//trim(dustname(id))
  end do

  if (ifprint) then
     print *,'-----------------'
     print *,'Checking fmass'
     print *,'Total disk: ',sum(fmass,mask=is_disk(1:8))
     print *,'Total envelope: ',sum(fmass,mask=is_envelope(1:8))
     print *,'Total cavity: ',sum(fmass,mask=is_cavity(1:8))
     print *,'-----------------'
  end if

end subroutine reapar




