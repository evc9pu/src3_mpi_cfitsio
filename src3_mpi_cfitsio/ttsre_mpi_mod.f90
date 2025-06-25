module ttsre_mpi_mod

  implicit none
  include 'mpif.h'

! ... variables added for communication and log control in MPI
  character :: cfltmp*30,ext*4
  integer,target :: ibuf(4),ierr,iready,myid,numprocs,irc
  integer :: sender,status(MPI_STATUS_SIZE)
  integer :: COMM_REDUCE,myid_reduce,numprocs_reduce

! ... constants which control message tags
  integer :: tag_ready,tag_start

! ... some work variables
  integer :: icountbuf
!  real,pointer,dimension(:,:,:,:) :: Rscratch4

CONTAINS

   subroutine ttsre_mpi_init
     implicit none
    
     iready = 1
     tag_ready = 0
     tag_start = 1
   end subroutine ttsre_mpi_init
 
end module ttsre_mpi_mod
