subroutine findopacid()

  use grid_mod
  use dust_mod
  use tts_mod, only:ifprint
  implicit none

  integer :: ir,it,ip,id,tt,ii,j,jn,ntt,itemp,dcount
  real(8) :: temp,dens,temp1,temp2,temp3
  real(8) :: gd(ndenmax)
  real(8) :: Tmin
!  integer :: ntden(ntemper)

  Tmin=1400.d0

  do ir=1,nrg
     do it=1,ntg
        do ip=1,npg
           do id=1,ndg+2

              if (id.gt.ndg) then

!                 temp=tdust2(ir,it,ip,id)*11605.d0
                 temp=tdave(ir,it,ip)*11605.d0
                 if (temp.lt.Tmin.and.densarr(ir,it,ip,9).gt.0.d0) then
                    print*,'check temperature in findopacid'
                    stop
                 end if
                 temp=log10(temp)/0.025d0
                 
                 if (id.eq.9) then

! 1:127, 2:128, 3:130, 4:132, 5:134, 6:136, 7:140

                    if (temp.lt.127.5d0) then
                       tt=1
                    else if (temp.lt.129.d0) then
                       tt=2
                    else if (temp.lt.131.d0) then
                       tt=3
                    else if (temp.lt.133.d0) then
                       tt=4
                    else if (temp.lt.135.d0) then
                       tt=5
                    else
                       temp=(temp-138.d0)/4.d0
                       tt=int(temp)+6
                       if (tt.gt.ntemper) then
                          print*,'gas T is too high',tt,ntemper
                          tt=ntemper
                       end if
                    end if

                 else
                    
                    if (temp.lt.133.d0) then
                       tt=4 !make 132 the lowest temp for gas opacity
                    else if (temp.lt.135.d0) then
                       tt=5
                    else 
                       temp=(temp-138.d0)/4.d0
                       tt=int(temp)+6
                       if (tt.gt.ntemper) then
                          print*,'gas T is too high',tt,ntemper
                          tt=ntemper
                       end if
                    end if

                 end if

                 do ii=1,ndenmax
                    gd(ii)=log10(gasden(ii,tt))
                 end do

                 dens=densarr(ir,it,ip,id)
                 dens=log10(dens)
                 call locate(gd,nden(tt),dens,j)
                 if (j.eq.0) then
                    jn=1
                 else if (j.gt.nden(tt)-1) then
                    jn=nden(tt)
!                    if (dens.gt.gd(jn)+1.d0) print*,'gas density out of boundary',dens,tt
                 else if ((dens-gd(j)).lt.(gd(j+1)-dens)) then
                    jn=j
                 else
                    jn=j+1
                 end if

!                 ntt=0
!                 do itemp=1,ntemper
!                    ntt=ntt+nden(itemp)
!                    ntden(itemp)=ntt
!                 end do
     
                 findopac(ir,it,ip,id)=ndg+ntden(tt)-nden(tt)+jn
     
              else
                 
                 findopac(ir,it,ip,id)=id
                 
              end if

              if (findopac(ir,it,ip,id).gt.nopac) then
                 print*,'opac assignment error',ir,it,ip,id,findopac(ir,it,ip,id)
                 stop
              end if

           end do
        end do
     end do
  end do

  is_used=.false.
  do id=1,ndg+2
     if (id.le.ndg) then
        is_used(id)=.true.
     else
        do ir=1,nrg
           do it=1,ntg
              do ip=1,npg
                 is_used(findopac(ir,it,ip,id))=.true.
              end do
           end do
        end do
     end if
  end do
  
  dcount=0
  do id=1,nopac
     if (is_used(id)) dcount=dcount+1 
  end do
  if (ifprint) print*,dcount,'of dust files are used among',nopac
  
  return
     
end subroutine findopacid
