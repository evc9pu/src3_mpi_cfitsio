module type_nested_image

  implicit none
  save

  type nested_image
     logical(1) :: disabled = .true.
  end type nested_image

contains

  subroutine nested_image_setup(image,imsize,l_max,l_min,n_nu,nu_min,nu_max,factor)

    implicit none

    type(nested_image),intent(in) :: image
    integer,intent(in) :: imsize,n_nu,factor
    real,intent(in) :: l_max,l_min,nu_min,nu_max
    logical,save :: warned = .false.
    if(.not.warned) then
       print *,'Warning - nested image is only available with FITS output'
       warned = .true.
    end if
  end subroutine nested_image_setup

  subroutine nested_image_bin(image,x,y,nu,flux)
    implicit none
    type(nested_image),intent(in) :: image
    real,intent(in) :: x,y,nu,flux
  end subroutine nested_image_bin

  subroutine nested_image_write(image,filename)
    implicit none
    character(len=*),intent(in) :: filename
    type(nested_image),intent(in),target :: image
    logical,save :: warned = .false.
    if(.not.warned) then
       print *,'Warning - nested image is only available with FITS output'
       warned = .true.
    end if
  end subroutine nested_image_write

  subroutine nested_image_read(image,filename)
    implicit none
    character(len=*),intent(in) :: filename
    type(nested_image),intent(in),target :: image
    logical,save :: warned = .false.
    if(.not.warned) then
       print *,'Warning - nested image is only available with FITS output'
       warned = .true.
    end if
  end subroutine nested_image_read

  subroutine get_single_image(image,level,inu,array)
    implicit none
    type(nested_image),intent(in),target :: image
    integer,intent(in) :: level,inu
    real,intent(in) :: array(:,:)
    logical,save :: warned = .false.
    if(.not.warned) then
       print *,'Warning - nested image is only available with FITS output'
       warned = .true.
    end if
  end subroutine get_single_image

end module type_nested_image

