subroutine linemap() 

  use grid_mod
  use dust_mod
  use tauint_mod
  use taunum_mod
  use tts_mod
  use constants
  use tts_mod
  use opacin_mod
  use lineint_mod
  use output_mod
  use ttsre_mpi_mod

  implicit none

  real(8) :: peelthet,sint,cost,thet,rad
  real(8) :: vrange,radmin,dr
  real(8) :: xmap,ymap,rmaxm,pixsize

  real(8),pointer :: Aul(:),abund(:),b0(:)
  real(8),pointer :: jlevel(:),Eul(:),gammaul(:)
  real(8),pointer :: Inum(:,:)

  integer :: ipeelthet,iphot
  integer :: ir,it,ip,ir1,it1,ip1,ir0
  integer :: inv,ix,iy,ispec,peelid
  integer :: nxmap,npixel,ipixel
  integer :: npixmod,npixpro,ipixstart,ipixend
  integer :: icenter,ii

  character(len=30),pointer :: molname(:)
  character(len=10) :: ctheta,depsuffix

  logical :: ifout,ifdep

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if (ifprint) print*,'doing line radiation transfer now'

  nspec=2

  nv=201

  vrange=20.d0 ! in km/s

  nxmap=501

  allocate(molname(nspec),Aul(nspec),abund(nspec),b0(nspec))
  allocate(jlevel(nspec),Eul(nspec),gammaul(nspec))
  allocate(Ivthin(nv,nspec),Ivthick(nv,nspec),vnu(nv))
  allocate(linemap_thin(nxmap,nxmap,nv,nspec))
!  allocate(linemap_thick(nxmap,nxmap,nv,nspec))
  allocate(Inum(nxmap,nxmap))

  do ispec=1,nspec
     
     if (ispec.eq.1) then

        !C18O1-0
!!$        molname(ispec)='C18O1-0'
!!$        Aul(ispec)=6.266e-8 ! Einstein A
!!$        abund(ispec)=1.d-4/62./5.5 ! abundance
!!$        b0(ispec)=5.4891e10*h ! rotation constant in erg
!!$        jlevel(ispec)=1.d0
!!$        Eul(ispec)=(5.27-0.)*k ! energy
!!$        gammaul(ispec)=6.266e-8 ! see page 57 of Draine's book     

        !C18O2-1
        molname(ispec)='C18O2-1'
        Aul(ispec)=6.011e-7 ! Einstein A
        abund(ispec)=1.d-4/62./5.5 ! abundance
        b0(ispec)=5.4891e10*h ! rotation constant in erg
        jlevel(ispec)=1.d0
        Eul(ispec)=(15.81-5.27)*k ! energy
        gammaul(ispec)=6.266e-8+6.011e-7 ! see page 57 of Draine's book, not used  
        
     end if

     if (ispec.eq.2) then

        !C18O3-2
        molname(ispec)='C18O3-2'
        Aul(ispec)=2.172e-6 ! Einstein A
        abund(ispec)=1.d-4/62./5.5 ! abundance
        b0(ispec)=5.4891e10*h ! rotation constant in erg
        jlevel(ispec)=1.d0
        Eul(ispec)=(31.61-15.81)*k ! energy
        gammaul(ispec)=6.266e-8+6.011e-7+6.011e-7 ! see page 57 of Draine's book, not used 
        
     end if
  
  end do

  do inv=1,nv
     vnu(inv)=vrange/dble(nv-1)*dble(inv-(nv-1)/2-1)
  end do

  rmaxm=5000.*autors
  pixsize=2.d0*rmaxm/dble(nxmap)

  icenter=(nxmap+1)/2

  npixel=nxmap*nxmap

  do ii=1,2

     if (ii.eq.1) then
        ifdep=.true.
        depsuffix='_dep'
     else
        ifdep=.false.
        depsuffix='_nodep'
     end if
        
     do peelid=1,npeel

        peelthet=thete_arr(peelid)
        
        write(ctheta,'(F5.1)') peelthet*rad2deg

        linemap_thin=0.d0
        linemap_thick=0.d0
        Inum=0

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        npixmod=mod(npixel,numprocs)
        npixpro=(npixel-npixmod)/numprocs

        if (myid.lt.npixmod) then
           ipixstart=myid*(npixpro+1)+1
           ipixend=(myid+1)*(npixpro+1)
        else
           ipixstart=npixmod*(npixpro+1)+(myid-npixmod)*npixpro+1
           ipixend=npixmod*(npixpro+1)+(myid-npixmod+1)*npixpro
        end if

        do ipixel=ipixstart,ipixend

           ix=int(ipixel/nxmap)+1
           iy=ipixel-(ix-1)*nxmap

           xmap=pixsize*dble(ix-icenter)
           ymap=pixsize*dble(iy-icenter)

           if (xmap**2+ymap**2.gt.1.d0) then

              xp=-ymap*cos(peelthet)
              yp=xmap
              zp=ymap*sin(peelthet)

              rsq=xp**2+yp**2+zp**2
              rtot=sqrt(rsq)

              call locate(rarr,nrg,rtot,ir)
              if (rarr(ir).gt.rtot) ir=ir-1

              cost=zp/rtot
              thet=acos(cost)
              call locate(thetarr,ntg,thet,it)
              if (thetarr(it).gt.thet) it=it-1
              
              ip=1
                             
              ifout=.false.
    
              ir1=ir
              it1=it
              ip1=ip
              
              ux=-sin(peelthet)
              uy=0.d0
              uz=-cos(peelthet)

              call linepeel(ir1,it1,ip1,Aul,abund,b0,jlevel,Eul,gammaul,ifout,ifdep)

              ux=sin(peelthet)
              uy=0.d0
              uz=cos(peelthet)

              ifout=.true.

              Ivthin=0.d0
              Ivthick=0.d0

              call linepeel(ir1,it1,ip1,Aul,abund,b0,jlevel,Eul,gammaul,ifout,ifdep)

              Inum(ix,iy)=Inum(ix,iy)+1
              do inv=1,nv
                 linemap_thin(ix,iy,inv,:)=linemap_thin(ix,iy,inv,:)+Ivthin(inv,:)
!                 linemap_thick(ix,iy,inv,:)=linemap_thick(ix,iy,inv,:)+Ivthick(inv,:)
              end do
        
           end if

        end do
        
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        icountbuf=nxmap*nxmap
        call MPI_ALLREDUCE(MPI_IN_PLACE,Inum,icountbuf,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        icountbuf=nxmap*nxmap*nv*nspec
        call MPI_ALLREDUCE(MPI_IN_PLACE,linemap_thin,icountbuf,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!        call MPI_ALLREDUCE(MPI_IN_PLACE,linemap_thick,icountbuf,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

        do ix=1,nxmap
           do iy=1,nxmap
              if (Inum(ix,iy).gt.0) then
                 do inv=1,nv
                    do ispec=1,nspec
                       linemap_thin(ix,iy,inv,ispec)=linemap_thin(ix,iy,inv,ispec)/dble(Inum(ix,iy))
!                       linemap_thick(ix,iy,inv,ispec)=linemap_thick(ix,iy,inv,ispec)/dble(Inum(ix,iy))
                    end do
                 end do
              end if
           end do
        end do

        if (myid.eq.0) then

           do ispec=1,nspec

              open(unit=15,file='linemap_thin_'//trim(molname(ispec))//'_' &
                   //trim(adjustl(ctheta))//trim(adjustl(depsuffix))//'.unf',form='unformatted')
              write(15) nxmap,nxmap,nv
              write(15) rmaxm/autors,pixsize/autors
              write(15) vnu
              write(15) linemap_thin(:,:,:,ispec)
              close(15)

!!$              open(unit=15,file='linemap_thick_'//trim(molname(ispec))//'_' &
!!$                   //trim(adjustl(ctheta))//trim(adjustl(depsuffix))//'.unf',form='unformatted')
!!$              write(15) nxmap,nxmap,nv
!!$              write(15) rmaxm/autors,pixsize/autors
!!$              write(15) vnu
!!$              write(15) linemap_thick(:,:,:,ispec)
!!$              close(15)

           end do
        
        end if

     end do

  end do

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  return

end subroutine linemap
