module type_nested_image

  use lib_cfitsio
  use type_smart_image

  implicit none
  save

  type nested_image

     integer :: imsize
     real :: l_min,l_max

     integer :: n_nu
     real :: nu_min,nu_max

     integer :: n_levels
     real,pointer :: levels(:)

     integer :: factor

     type(smart_image),pointer :: image(:,:)

  end type nested_image

contains

  subroutine nested_image_setup(image,imsize,l_max,l_min,n_nu,nu_min,nu_max,factor)

    implicit none

    type(nested_image),intent(inout) :: image
    integer,intent(in) :: imsize,n_nu,factor
    real,intent(in) :: l_max,l_min,nu_min,nu_max
    integer :: j,l

    if(mod(imsize,factor).ne.0) stop "imsize must be multiple of factor"
    if(mod(imsize/2,2).ne.0)    stop "imsize/2 must be multiple of 2"

    image%factor = factor

    image%imsize = imsize
    image%l_min  = l_min
    image%l_max  = l_max

    image%n_nu   = n_nu
    image%nu_min = nu_min
    image%nu_max = nu_max

    call find_levels(l_min,l_max,factor,image%n_levels,image%levels)

    allocate(image%image(image%n_levels,n_nu))
    do l=1,image%n_levels
       do j=1,n_nu
          call smart_image_initialize(image%image(l,j),imsize)
       end do
    end do

  end subroutine nested_image_setup

  subroutine find_levels(l_min,l_max,factor,n_levels,levels)
    implicit none
    real,intent(in) :: l_min,l_max
    integer,intent(in) :: factor
    integer,intent(out) :: n_levels
    real,pointer,intent(out) :: levels(:)
    integer :: l
    n_levels = ceiling(log10(l_max/l_min) / log10(real(factor)))+1
    allocate(levels(n_levels))
    do l=1,n_levels
       levels(n_levels+1-l) = l_max / real(factor)**(l-1)
    end do
  end subroutine find_levels

  subroutine nested_image_bin(image,x,y,nu,flux)

    use math_binning

    implicit none

    type(nested_image),intent(inout) :: image
    real,intent(in) :: x,y,nu,flux
    integer :: l,ix,iy,inu

    if(max(abs(x),abs(y)) > image%l_max) return

    l = locate(image%levels,max(abs(x),abs(y)))
    if(l < 0) then
       l=1
    else
       l=l+1
    end if

    ix  = ipos(-image%levels(l),+image%levels(l),x,image%imsize)
    iy  = ipos(-image%levels(l),+image%levels(l),y,image%imsize)
    inu = ipos(image%nu_min,image%nu_max,nu,image%n_nu)
    if(inu>=1.and.inu<=image%n_nu) then
       call smart_image_insert(image%image(l,inu),ix,iy,flux)
    end if

  end subroutine nested_image_bin

  subroutine nested_image_write(image,filename)

    implicit none

    character(len=*),intent(in) :: filename

    type(nested_image),intent(inout),target :: image
    type(smart_image),pointer :: im
    integer,pointer :: n_levels,n_nu,imsize

    integer :: n_sparse,n_array
    real,pointer :: cube(:,:,:)
    integer,pointer :: ia_array(:,:)
    integer :: ca,cs,row,l,j
    integer :: unit

    n_sparse = count(image%image(:,:)%sparse)
    n_array  = size(image%image)-n_sparse

    imsize   => image%imsize
    n_nu     => image%n_nu
    n_levels => image%n_levels

    allocate(cube(imsize,imsize,n_array))
    allocate(ia_array(imsize+1,n_sparse))

    ca = 0
    cs = 0
    do l=1,n_levels
       do j=1,n_nu
          im => image%image(l,j)
          if(im%sparse) then
             cs = cs + 1
             ia_array(:,cs) = im%s_array%ia(:)
          else
             ca = ca + 1
             cube(:,:,ca) = im%array(:,:)
          end if
       end do
    end do

    call fits_open_new(unit,filename,confirm=.false.)

    call fits_write_primary_header(unit,-32,(/imsize,imsize,n_array/),.true.)
    call fits_write_array(unit,cube)

    call fits_write_keyword(unit,'IMSIZE',imsize)
    call fits_write_keyword(unit,'N_NU',n_nu)
    call fits_write_keyword(unit,'N_ARRAY',n_array)
    call fits_write_keyword(unit,'N_SPARSE',n_sparse)
    call fits_write_keyword(unit,'N_LEVELS',n_levels)
    call fits_write_keyword(unit,'LOGNUMIN',image%nu_min)
    call fits_write_keyword(unit,'LOGNUMAX',image%nu_max)
    call fits_write_keyword(unit,'LMIN',image%l_min)
    call fits_write_keyword(unit,'LMAX',image%l_max)
    call fits_write_keyword(unit,'FACTOR',image%factor)

    call fits_create_hdu(unit,2)
    call fits_write_primary_header(unit,16,(/imsize+1,n_sparse/),.true.)
    call fits_write_array(unit,ia_array)

    call fits_create_hdu(unit,3)
    call fits_table_write_header(unit,0,4,(/'L     ','NU    ','SPARSE','M     '/),(/'J','J','L','J'/),(/'','','',''/),"Summary")

    row = 0
    do l=1,n_levels
       do j=1,n_nu
          im => image%image(l,j)
          row = row + 1
          call fits_table_write_column(unit,'L',l,row=row)
          call fits_table_write_column(unit,'NU',j,row=row)
          call fits_table_write_column(unit,'SPARSE',im%sparse,row=row)
          if(im%sparse) then
             call fits_table_write_column(unit,'M',im%m,row=row)
          else
             call fits_table_write_column(unit,'M',0)
          end if
       end do
    end do

    call fits_create_hdu(unit,4)
    call fits_table_write_header(unit,0,2,(/'JA','A '/),(/'I','E'/),(/'',''/),"Sparse elements")
    row = 1
    do l=1,n_levels
       do j=1,n_nu
          im => image%image(l,j)
          if(im%sparse) then
             call fits_table_write_column(unit,'JA',im%s_array%ja(1:im%m),row=row)
             call fits_table_write_column(unit,'A',im%s_array%a(1:im%m),row=row)
             row = row + im%m
          end if
       end do
    end do

    call fits_close(unit)

  end subroutine nested_image_write

  subroutine nested_image_read(image,filename)

    implicit none

    character(len=*),intent(in) :: filename

    type(nested_image),intent(out),target :: image
    type(smart_image),pointer :: im
    integer,pointer :: n_levels,n_nu,imsize
    integer :: n_sparse,n_array

    real,pointer :: cube(:,:,:)
    integer,pointer :: ia_array(:,:)
    integer :: ca,cs,row,l,j,unit
    call fits_open_read(unit,filename)

    imsize   => image%imsize
    n_nu     => image%n_nu
    n_levels => image%n_levels

    call fits_read_keyword(unit,'IMSIZE',imsize)
    call fits_read_keyword(unit,'N_NU',n_nu)
    call fits_read_keyword(unit,'N_ARRAY',n_array)
    call fits_read_keyword(unit,'N_SPARSE',n_sparse)
    call fits_read_keyword(unit,'N_LEVELS',n_levels)
    call fits_read_keyword(unit,'LOGNUMIN',image%nu_min)
    call fits_read_keyword(unit,'LOGNUMAX',image%nu_max)
    call fits_read_keyword(unit,'LMIN',image%l_min)
    call fits_read_keyword(unit,'LMAX',image%l_max)
    call fits_read_keyword(unit,'FACTOR',image%factor)

    allocate(image%image(n_levels,n_nu))
    call find_levels(image%l_min,image%l_max,image%factor,image%n_levels,image%levels)

    allocate(cube(imsize,imsize,n_array))
    allocate(ia_array(imsize+1,n_sparse))

    call fits_read_array(unit,cube)

    call fits_move_hdu(unit,2)
    call fits_read_array(unit,ia_array)

    call fits_move_hdu(unit,3)

    row = 0
    cs = 0
    ca = 0

    do l=1,n_levels
       do j=1,n_nu

          im => image%image(l,j)

          row = row + 1
          call fits_table_read_column(unit,'SPARSE',im%sparse,row=row)
          call fits_table_read_column(unit,'M',im%m,row=row)

          if(im%sparse) then
             cs = cs + 1
             call sparse_matrix_initialize(im%s_array,imsize,im%m)
             im%s_array%ia(:) = ia_array(:,cs)
          else
             ca = ca + 1
             allocate(im%array(imsize,imsize))
             im%array(:,:) = cube(:,:,ca)
          end if

       end do
    end do

    call fits_move_hdu(unit,4)
    row = 1
    do l=1,n_levels
       do j=1,n_nu
          im => image%image(l,j)
          if(im%sparse) then
             call fits_table_read_column(unit,'JA',im%s_array%ja,row=row)
             call fits_table_read_column(unit,'A',im%s_array%a,row=row)
             row = row + im%m
          end if
       end do
    end do

    call fits_close(unit)

  end subroutine nested_image_read

  subroutine get_single_image(image,level,inu,array)

    implicit none

    type(nested_image),intent(in),target :: image
    type(smart_image),pointer :: im
    integer,intent(in) :: level,inu
    real,intent(out) :: array(:,:)
    real,pointer :: array_sub(:,:)
    integer :: l,i,ixt,iyt,imin,imax,jmin,jmax
    real :: valuet

    allocate(array_sub(image%imsize/image%factor,image%imsize/image%factor))

    do l=1,level

       if(l > 1) call rebin(array,array_sub,image%factor)

       im => image%image(l,inu)

       if(im%sparse) then
          array = 0.
          do i=1,im%m
             call element(im%s_array,i,ixt,iyt,valuet)
             array(ixt,iyt) = valuet
          end do
       else
          array = im%array
       end if

       if(l > 1) then

          imin = image%imsize/2-image%imsize/image%factor/2+1
          imax = image%imsize/2+image%imsize/image%factor/2
          jmin = imin
          jmax = imax

          array(imin:imax,jmin:jmax) = array_sub

       end if

    end do

  end subroutine get_single_image

  subroutine rebin(array_in,array_out,factor)
    implicit none
    real,intent(in)  :: array_in(:,:)
    real,intent(out) :: array_out(:,:)
    integer,intent(in) :: factor
    integer :: i,j,ii,jj
    array_out = 0.
    do i=1,size(array_out,1)
       do j=1,size(array_out,2)
          do ii=1+(i-1)*factor,i*factor
             do jj=1+(j-1)*factor,j*factor
                array_out(i,j) = array_out(i,j) + array_in(ii,jj)
             end do
          end do
       end do
    end do
  end subroutine rebin

end module type_nested_image

