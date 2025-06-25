module smatrix

  type sparse_matrix
     integer :: m = 0
     integer :: N = 0
     integer :: rowlen = 0
     integer,pointer :: ia(:)
     integer,pointer :: ja(:)
     real,pointer    :: a(:)
  end type sparse_matrix

  logical :: debug = .false.

contains

  subroutine sparse_matrix_destroy(matrix)
    implicit none
    type(sparse_matrix),intent(inout) :: matrix
    deallocate(matrix%ia)
    deallocate(matrix%ja)
    deallocate(matrix%a)
  end subroutine sparse_matrix_destroy

  subroutine sparse_matrix_initialize(matrix,size,rowlen) ! initialize the matrix

    implicit none

    type(sparse_matrix),intent(inout) :: matrix
    integer,intent(in)                :: size
    integer,intent(in),optional       :: rowlen
    if(debug) print *,"[sparse_matrix] initializing"

    ! Set dimension to 'size'
    matrix%N      = size

    ! Set initial row length to 'size'
    if(present(rowlen)) then
       matrix%rowlen = rowlen
    else
       matrix%rowlen = 1
    end if

    ! Allocate the arrays
    allocate(matrix%ia(size+1))        ; matrix%ia = 1
    allocate(matrix%ja(matrix%rowlen)) ; matrix%ja = 0
    allocate(matrix%a(matrix%rowlen))  ; matrix%a  = 0

  end subroutine sparse_matrix_initialize

  subroutine grow(matrix) ! double the size of the JA and A arrays

    implicit none

    type(sparse_matrix),intent(inout) :: matrix
    integer :: ja_temp(matrix%rowlen)
    real    :: a_temp(matrix%rowlen)

    if(debug) print *,"[sparse_matrix] doubling array sizes"

    ! temporarily store contents of JA and A arrays
    ja_temp = matrix%ja
    a_temp  = matrix%a

    ! double the size of the JA array
    deallocate(matrix%ja)
    allocate(matrix%ja(matrix%rowlen*2))
    matrix%ja = 0
    matrix%ja(1:matrix%rowlen) = ja_temp

    ! double the size of the A array
    deallocate(matrix%a)
    allocate(matrix%a(matrix%rowlen*2))
    matrix%a  = 0.
    matrix%a(1:matrix%rowlen)  =  a_temp

    ! set new row length
    matrix%rowlen = matrix%rowlen * 2

  end subroutine grow

  subroutine sparse_matrix_insert(matrix,ix,iy,value) ! insert a value into the sparse matrix

    implicit none

    type(sparse_matrix),intent(inout) :: matrix
    integer,intent(in)                :: ix,iy
    real,intent(in)                   :: value

    integer :: istart,iend,i,ipos

    logical :: new

    if(ix < 1 .or. ix > matrix%N .or. iy < 1 .or. iy > matrix%N) return

    ! Find first and last element in the row
    istart = matrix%ia(iy)
    iend   = matrix%ia(iy+1)-1
    ipos   = istart
    new    = .true.

    ! Cycle through, and stop when greater than current element (insert at i+1)
    do i=iend,istart,-1
       if(ix==matrix%ja(i)) then
          ipos = i
          new = .false.
          exit
       end if
       if(ix > matrix%ja(i)) then
          ipos = i+1
          new = .true.
          exit
       end if
    end do

    if(new) then

       if(matrix%m==matrix%rowlen) call grow(matrix)

       if(debug) print *,"[sparse_matrix] inserting element"

       ! shift JA array to make place for new value
       do i=matrix%m,ipos,-1
          matrix%ja(i+1) = matrix%ja(i)
          matrix%a(i+1)  = matrix%a(i)
       end do

       ! insert new value
       matrix%ja(ipos) = ix
       matrix%a(ipos)  = value

       ! increment element number
       matrix%m = matrix%m + 1

       ! increment IA array
       do i=iy+1,matrix%N+1
          matrix%ia(i) = matrix%ia(i) + 1
       end do

    else

       if(debug) print *,"[sparse_matrix] adding element"

       ! add value to existing element
       matrix%a(ipos) = matrix%a(ipos) + value

    end if

  end subroutine sparse_matrix_insert

  subroutine show(matrix)
    implicit none
    type(sparse_matrix),intent(inout) :: matrix
    print '('//repeat('I2,1X,',size(matrix%ia))//'1X)',matrix%ia(:)
    print '('//repeat('I2,1X,',matrix%m)//'1X)',matrix%ja(1:matrix%m)
    print *,matrix%a(1:matrix%m)
  end subroutine show

  subroutine show_full(matrix)
    implicit none
    type(sparse_matrix),intent(inout) :: matrix
    integer :: ix,iy
    integer :: istart,iend
    real :: line(matrix%N)
    do iy=1,matrix%N
       line = 0.
       istart = matrix%ia(iy)
       iend   = matrix%ia(iy+1)-1
       do ix=istart,iend
          line(matrix%ja(ix)) = matrix%a(ix)
       end do
       print '('//repeat('F4.2,1X,',matrix%N)//'1X)',line
    end do
  end subroutine show_full

  subroutine element(matrix,pos,ix,iy,value)
    implicit none
    type(sparse_matrix),intent(in) :: matrix
    integer,intent(in) :: pos
    integer,intent(out) :: ix,iy
    real,intent(out) :: value
    integer :: i
    do i=1,matrix%N+1
       if(matrix%ia(i) > pos) then
          iy = i-1
          exit
       end if
    end do
    ix    = matrix%ja(pos)
    value = matrix%a(pos)
  end subroutine element

end module smatrix



