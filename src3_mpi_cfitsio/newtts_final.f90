program tts

  ! radiative transfer in a disk+envelope around a star.
  ! height of disk, z, goes as z1*r**b.
  ! density goes as c*r**-a.
  ! opacity is dust with either rayleigh or H-G phase function.
  ! output is stokes parameters as a function of angle
  !
  ! history:
  ! 00/03/19 (mjw):  update version stamp for changes on 03/19.  Also,
  ! change version name to VGER (from mctherm)

  use tts_mod
  use grid_mod
  use output_mod
  use constants
  use user_set_generator
  use ttsre_mpi_mod

  implicit none

  integer :: iwav, ipeelorig,nporig,iter,ncells,i1
  real(8) :: isrf_int,cpusec,cpusec_old,cpusec_start,delta_cpusec
  character(len=80) :: atname,dustname(ndg),cfllog,cmsgnm
  character(len=3) :: converge
  integer :: ir,it,ip

! ... parallel version stamp and processor info
  character :: CPAR_STAMP*40,CPROC*40
  integer :: i,len_proc
  data CPAR_STAMP/'PARALLEL VERSION 1.1a [2012.03.22]'/

! ... misc scalars
  integer :: np_reset

! define a few useful things
  call ttsre_mpi_init()

! ... Initialize MPI for this process, getting MPI_COMM_WORLD process 
!     ID and total number of processes.
  call MPI_INIT( ierr )
  call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )


! create file tag from MYID (for child processes) and open log file for output
  if (myid .eq. 0) then ! Parent PROCESS
     write(cfllog,'(a)') 'p_disk-parent.log'
! create new communicator for child processes, so pass undefined color for parent
     call MPI_COMM_SPLIT(MPI_COMM_WORLD,MPI_UNDEFINED,&
          myid-1,COMM_REDUCE,ierr)
  else ! Child PROCESSES
     if (myid .lt. 10000) then
        write(ext,'(i4.4)') myid
     else
        stop 'HA HA HA...increase logic block of CFLLOG'
     endif
     cfllog = 'p_disk-child'//ext//'.log'
! Create new communicator for child processes
! myid=1 will become myid_reduce=0
     call MPI_COMM_SPLIT(MPI_COMM_WORLD,0,myid-1,COMM_REDUCE,ierr)
  endif

  open(unit=12,file=cfllog,status='unknown')
  write(12,'(a,a)') version,' - ttsre version'
  write(12,'(a,i4,a,i4,a)') &
       '(MPI_COMM_WORLD): Process number ', myid, ' of ', &
       numprocs,' total (proc. no. start with 0)'

  if (numprocs.gt.31) then
     write(0,*) 'CANNOT USE MORE THAN 31 PROCESSES AT THIS TIME.'
     goto 0990
  end if

  if (myid.gt.0) then
     call MPI_COMM_RANK( COMM_REDUCE, myid_reduce, ierr )
     call MPI_COMM_SIZE( COMM_REDUCE, numprocs_reduce, ierr )
     write(12,'(a,i4,a,i4,a)') &
          '(COMM_REDUCE): Process number ', myid_reduce, ' of ',&
          numprocs_reduce,' total (proc. no. start with 0)'
  endif

! ... get processor info and writeout version stamp
  call MPI_GET_PROCESSOR_NAME(CPROC, len_proc, ierr)
  CPROC =  'Processor = '//CPROC(1:28)
  write(12,'(a)') CPROC
  write(12,'(a)') CPAR_STAMP


!**********************************************************************
!**********************************************************************
!**********************************************************************
!**********************************************************************
!**********************************************************************

  if (myid.eq.0) then
     ifprint=.true.
  else
     ifprint=.false.
  end if

  ! input and setup
  call setup(atname,dustname,i1)

  write(cmsgnm,'(a,i4)') 'Setting Random Number Generator for stream ',myid+1
  call user_set_all(12345678-i1,87654321,myid+1)

  call wrimsg('TTS',cmsgnm)

  ipeelorig=ipeel
  nporig=np

  write(cmsgnm,'(a,i9)') 'Reducing NP so that NP_new * num_processors = NP_old =  ', np
  call wrimsg('TTS',cmsgnm)
  np = np / numprocs

  npsav=np
  if (ifprint) print*,'np,nporig',np,nporig
  write(12,*) 'np,nporig',np,nporig
  iterstart=1
  iter=0

  if (ifprint) print*,'calling dustinit'
  call dustinit(dustname)
  if (ifprint) print*,'done with dustinit, calling draine_cdf'

  call draine_cdf() ! small grains
  call isrf_cdf(isrf_int,isrf_Av)   !need this for energy density normalization
  if (ifprint) print*, 'tts, done with draine_cdf'
  if (ifprint) print*,'tts, done with isrf_cdf'
  if (iout.eq.1) then
     l_isrf=isrf_int*4.d0*pi*(rmax*rsol*rstar)**2/2.d0*isrf_scl
     if (ifprint) print*,'isrtf_int,rmax*rsol*rstar,isrf_scl'
     if (ifprint) print*,isrf_int,rmax*rsol*rstar,isrf_scl
     if (ifprint) print*,'l_isrf (lsun)',l_isrf/lsun
  else
     l_isrf=0.d0
  end if

  call atmosinit(atname)
  if (ifprint) print*,'done with atmosinit'

  iwav=1

  write(12,*) 'Random number seed barrier...'
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if (myid.eq.0) then !Parent process

     write(12,'(a,i8)') 'PARENT>> base random number seed = ', i1
     do i=1,numprocs-1
        call MPI_RECV(iready,1,MPI_INTEGER,MPI_ANY_SOURCE,&
             tag_ready,MPI_COMM_WORLD, status, ierr)
        sender = status(MPI_SOURCE)
        ibuf(1) = i1*((sender+1)**3)
        write(cmsgnm,'(a,i11,a,i4)') &
             'PARENT>> Sending random number seed (i1) = ',&
             ibuf(1),' to process ',sender
        call wrimsg('TTS',cmsgnm)
        call MPI_SEND(ibuf,4,MPI_INTEGER,sender,&
             tag_start,MPI_COMM_WORLD,ierr)
     end do

  else !Child processes

     call MPI_SEND(iready,1,MPI_INTEGER,0,tag_ready,&
          MPI_COMM_WORLD, ierr)
     write(12,'(a)') 'CHILD>> Sending IREADY signal to Parent'
     call MPI_RECV(ibuf,4,MPI_INTEGER,0,&
          tag_start,MPI_COMM_WORLD, status, ierr)
     i1 = -abs(ibuf(1))
     write(cmsgnm,'(a,i11)') &
          'CHILD>> Receiving random number seed (i1) from parent: ',i1
     call wrimsg('TTS',cmsgnm)

  end if

  call cpu_time(cpusec_start)

  call initarr_final(.true.)
  call setup_wave()

  converge='yes'

  ! final iteration (or first if ilucy=0)
  ipeel=ipeelorig
  np=nporig
  write(cmsgnm,'(a,i9)') 'Reducing NP so that NP_new * num_processors = NP_old =  ', np
  call wrimsg('TTS',cmsgnm)
  np = np / numprocs

  iter=iter+1

  call cpu_time(cpusec)
  write(cmsgnm,'(a,f11.2)') 'Initializing images and calling DISK at CPU time ',cpusec
  call wrimsg('TTS',cmsgnm)
  if (ifprint) print*,'Initializing images and calling DISK at CPU time ',cpusec

  call initimages()

  call disk(converge,iter,i1)

  cpusec_old = cpusec
  call cpu_time(cpusec)
  write(cmsgnm,'(a,f11.2,a)') 'Elapsed time for final RT call to DISK',cpusec-cpusec_old,' (sec CPU time)'
  call wrimsg('TTS',cmsgnm)
  if (ifprint) print*,'Elapsed time for final RT call to DISK',cpusec-cpusec_old,' (sec CPU time)'

  ! only parent process needs to continue, all others go to barrier and wait
  if (myid > 0) goto 0990
 

  ! convert stokes arrays to flux arrays
  call diskflux()
  write(cmsgnm,'(a)') 'tts, done with diskflux, calling output.'
  call wrimsg('TTS',cmsgnm)
  if (ifprint) print*,'tts, done with diskflux, calling output'

  call output()

  call cpu_time(cpusec)
  write(cmsgnm,'(a,f11.2,a)') 'Elapsed CPU time for execution of main program: ',&
       cpusec-cpusec_start,' (sec)'
  call wrimsg('TTS',cmsgnm)
  if (ifprint) print*,'Elapsed CPU time for execution of main program: ',&
       cpusec-cpusec_start,' (sec)'

! ... Clean up MPI baggage for this process
0990 call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call MPI_FINALIZE(irc)
  if (myid.eq.0) then
     WRITE(12,'(a)') 'Parent process terminating....'
  else
     WRITE(12,'(a)') 'Child process terminating....'
  endif
  WRITE(12,'(a)') 'Have a nice day.'
  close(12)

end program tts
