module type_smart_image

  use smatrix
  implicit none

  type smart_image
     integer :: n,m
     logical :: sparse
     type(sparse_matrix) :: s_array
     real,pointer    :: array(:,:)
  end type smart_image

contains

  subroutine smart_image_initialize(image,imsize)

    implicit none

    type(smart_image),intent(inout) :: image
    integer,intent(in) :: imsize

    image%n = imsize
    image%m = 0
    image%sparse = .true.
    call sparse_matrix_initialize(image%s_array,imsize)

  end subroutine smart_image_initialize

  subroutine smart_image_insert(image,ix,iy,value)

    implicit none

    type(smart_image),intent(inout) :: image
    integer,intent(in) :: ix,iy
    real,intent(in) :: value
    integer :: i,ixt,iyt
    real :: valuet

    if(image%sparse) then

       call sparse_matrix_insert(image%s_array,ix,iy,value)
       image%m = image%s_array%m

       if(image%m > 0.50 * image%n**2) then
          print *,"[smart_image] switching to image array"
          allocate(image%array(image%n,image%n))
          image%array(image%n,image%n) = 0.
          do i=1,image%m
             call element(image%s_array,i,ixt,iyt,valuet)
             image%array(ixt,iyt) = valuet
          end do
          image%sparse = .false.
          call sparse_matrix_destroy(image%s_array)
       end if

    else

       image%array(ix,iy) = image%array(ix,iy) + value

    end if

  end subroutine smart_image_insert

end module type_smart_image

