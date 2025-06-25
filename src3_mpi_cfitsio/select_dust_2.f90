subroutine select_dust(ir,it,ip,mask,idust2)

  use random
  use grid_mod
  use opacin_mod

  implicit none

  integer,intent(in) :: ir,it,ip
  integer,intent(out) :: idust2
  integer :: id,nn
  logical,intent(in) :: mask(ndg)

  real(8) :: xran,opac,frac

  xran=ran()
  opac=0.d0
  frac=0.d0
  idust2=0

  nn=0
  do id=1,ndg+2
     if (.not.is_sg(id)) then
        if (densarr(ir,it,ip,id).gt.0.d0) then
           nn=nn+1
           idust2=id
        end if
     end if
  end do
  if (nn.eq.1) then
     return
  else
     print*,'more than one non-zero density',ir,it,nn
  end if

  print*,'select_dust'
  idust2=0
  do id=1,ndg+2
     if(mask(id)) &
          & opac=opac+kapd(findopac(ir,it,ip,id))*densarr(ir,it,ip,id)
  end do

  if(opac.eq.0) then
    print *,'ERROR: cannot select dust type from cell, total opacity is zero'
    print *,ir,it,ip
    print *,'kapd: ',kapd
    print *,'dens: ',densarr(ir,it,ip,:)
    stop
  end if

  do id=1,ndg+2
     if(mask(id)) then
     frac=frac+kapd(findopac(ir,it,ip,id))*densarr(ir,it,ip,id)/opac
     if (xran.le.frac + 1.d-10) then
        idust2=id
        exit
     end if
     end if
  end do

  if (idust2.eq.0) then
     print*,'idust2 not set'
     print*,'frac,opac',frac,opac,ir,it,ip
     print*,'kapd',kapd
     print*,'density',densarr(ir,it,ip,1:5)
     print*,'if it is roundoff, okay comment out the stop command'
     print*,'stopping program for now'
     frac = 0.d0
     do id=1,ndg+2
        if(mask(id)) then
        frac=frac+kapd(findopac(ir,it,ip,id))*densarr(ir,it,ip,id)/opac
        print '(F40.35)',xran
        print '(F40.35)',frac
        if (xran.le.frac + 1.e-10) then
           idust2=id
           exit
        end if
        end if
     end do
     print*,idust2
     stop
  end if

end subroutine select_dust
