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

  integer :: ndirec,idirec,ipeelthet,iphot
  integer :: ir,it,ip,ir1,it1,ip1,ir0
  integer :: inv,nr,nthet,iir,iit,ix,iy
  integer :: ispec
  integer :: nrmod,nrpro,irstart,irend
  integer :: peelid
  integer :: nxmap

  character(len=30),pointer :: molname(:)
  character(len=10) :: ctheta

  logical :: ifout,ifdep

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ifdep=.false.
  
  if (ifprint) print*,'doing line radiation transfer now'

  nspec=1

  nv=201

  vrange=500.d0 ! in km/s

  nxmap=1001 !

  allocate(molname(nspec),Aul(nspec),abund(nspec),b0(nspec))
  allocate(jlevel(nspec),Eul(nspec),gammaul(nspec))
  allocate(Ivthin(nv,nspec),Ivthick(nv,nspec),vnu(nv))
  allocate(linemap_thin(nxmap,nxmap,nv,nspec))
  allocate(linemap_thick(nxmap,nxmap,nv,nspec))
  allocate(Inum(nxmap,nxmap))

  do ispec=1,nspec
     
     if (ispec.eq.1) then

!!$        !C18O1-0
!!$        molname(ispec)='C18O1-0'
!!$        Aul(ispec)=6.266e-8 ! Einstein A
!!$        abund(ispec)=1.d-4/62./5.5 ! abundance
!!$        b0(ispec)=5.4891e10*h ! rotation constant in erg
!!$        jlevel(ispec)=1.d0
!!$        Eul(ispec)=(5.27-0.)*k ! energy
!!$        gammaul(ispec)=6.266e-8 ! see page 57 of Draine's book     

        !CO2-1
        molname(ispec)='CO2-1'
        Aul(ispec)=6.910d-7 
        abund(ispec)=1.d-4
        b0(ispec)=2.77*k
        jlevel(ispec)=3.d0
        Eul(ispec)=(16.60-5.53)*k
        gammaul(ispec)=6.910d-7+7.203d-08

!!$     else if (ispec.eq.2) then
!!$        
!!$        !13CO2-1
!!$        molname(ispec)='13CO3-2'
!!$        Aul(ispec)=2.181d-6
!!$        abund(ispec)=1.d-4/50.d0
!!$        b0(ispec)=2.64*k
!!$        jlevel(ispec)=3.d0
!!$        Eul(ispec)=(31.73-15.87)*k
!!$        gammaul(ispec)=2.181d-6+6.038d-7
!!$
!!$     else if (ispec.eq.3) then
!!$
!!$        !C18O2-1
!!$        molname(ispec)='C18O3-2'
!!$        Aul(ispec)=2.172d-6
!!$        abund(ispec)=1.d-4/327.d0
!!$        b0(ispec)=2.6355*k
!!$        jlevel(ispec)=3.d0
!!$        Eul(ispec)=(31.61-15.81)*k
!!$        gammaul(ispec)=2.172d-6+6.011d-7
!!$        
!!$     else if (ispec.eq.4) then
!!$
!!$        !SiO5-4
!!$        molname(ispec)='SiO7-6'
!!$        Aul(ispec)=1.4638d-3
!!$        abund(ispec)=1.d-10
!!$        b0(ispec)=1.042*k
!!$        jlevel(ispec)=7.d0
!!$        Eul(ispec)=(58.35-43.76)*k
!!$        gammaul(ispec)=1.4638d-3+9.1161d-4
!!$        
!!$     else if (ispec.eq.5) then
!!$        
!!$        !SiO8-7
!!$        molname(ispec)='SiO8-7'
!!$        Aul(ispec)=2.2030d-3
!!$        abund(ispec)=1.d-10
!!$        b0(ispec)=1.042*k
!!$        jlevel(ispec)=8.d0
!!$        Eul(ispec)=(75.02-58.35)*k
!!$        gammaul(ispec)=2.203d-3+1.4638d-3
        
     end if
  
  end do

  do inv=1,nv
     vnu(inv)=vrange/dble(nv-1)*dble(inv-(nv-1)/2-1)
  end do

  rmaxm=5000.*autors
  pixsize=2.d0*rmaxm/dble(nxmap)

  do peelid=1,npeel

     peelthet=thete_arr(peelid)

     write(ctheta,'(F5.1)') peelthet*rad2deg

     radmin=1.d0/sin(peelthet)
     call locate(ravearr,nrg-1,radmin,ir0)
     ir0=ir0+1
     if (ifprint) print*,'thet,ir0',peelthet*180.d0/pi,radmin,ir0

     linemap_thin=0.d0
     linemap_thick=0.d0
     Inum=0

     nr=int(rmaxm/sin(peelthet)/autors/100.d0)

     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

     nrmod=mod(nr,numprocs)
     nrpro=(nr-nrmod)/numprocs

     if (myid.lt.nrmod) then
        irstart=myid*(nrpro+1)+1
        irend=(myid+1)*(nrpro+1)
     else
        irstart=nrmod*(nrpro+1)+(myid-nrmod)*nrpro+1
        irend=nrmod*(nrpro+1)+(myid-nrmod+1)*nrpro
     end if

     print*,myid,irstart,irend

     do iir=irstart,irend

        print*,'myid,iir',myid,iir

        rad=100.d0*dble(iir)*autors

        call locate(rarr,nrg,rad,ir)

        if (rarr(ir).gt.rad) ir=ir-1

        nthet=int(180.d0/1.d0)

        do iit=1,nthet-1
           
           thet=pi/180.d0*dble(iit)*1.d0

           call locate(thetarr,ntg,thet,it)
           
           if (thetarr(it).gt.thet) it=it-1

           ip=1
           
!========================================
           xp=0.d0
           sint=sin(thet)
           cost=cos(thet)
           yp=rad*sint
           zp=rad*cost
           
           rsq=xp**2+yp**2+zp**2
           rtot=sqrt(rsq)
                  
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

           ymap=zp*sin(peelthet)-xp*cos(peelthet)
           xmap=yp

           ix=int((xmap+rmaxm)/pixsize)+1
           iy=int((ymap+rmaxm)/pixsize)+1

           if (ix.ge.1.and.ix.le.nxmap.and.iy.ge.1.and.iy.le.nxmap) then
              Inum(ix,iy)=Inum(ix,iy)+1
              do inv=1,nv
                 linemap_thin(ix,iy,inv,:)=linemap_thin(ix,iy,inv,:)+Ivthin(inv,:)
                 linemap_thick(ix,iy,inv,:)=linemap_thick(ix,iy,inv,:)+Ivthick(inv,:)
              end do
           end if
!========================================
           xp=0.d0
           sint=sin(thet)
           cost=cos(thet)
           yp=-rad*sint
           zp=rad*cost
           
           rsq=xp**2+yp**2+zp**2
           rtot=sqrt(rsq)
                  
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

           ymap=zp*sin(peelthet)-xp*cos(peelthet)
           xmap=yp

           ix=int((xmap+rmaxm)/pixsize)+1
           iy=int((ymap+rmaxm)/pixsize)+1

           if (ix.ge.1.and.ix.le.nxmap.and.iy.ge.1.and.iy.le.nxmap) then
              Inum(ix,iy)=Inum(ix,iy)+1
              do inv=1,nv
                 linemap_thin(ix,iy,inv,:)=linemap_thin(ix,iy,inv,:)+Ivthin(inv,:)
                 linemap_thick(ix,iy,inv,:)=linemap_thick(ix,iy,inv,:)+Ivthick(inv,:)
              end do
           end if
!========================================

        end do
        
     end do

     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

     if (myid.eq.0) then 
        print*,'done calculation, output line map'
     end if

     icountbuf=nxmap*nxmap
     call MPI_ALLREDUCE(MPI_IN_PLACE,Inum,icountbuf,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
     icountbuf=nxmap*nxmap*nv*nspec
     call MPI_ALLREDUCE(MPI_IN_PLACE,linemap_thin,icountbuf,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE,linemap_thick,icountbuf,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

     do ix=1,nxmap
        do iy=1,nxmap
           if (Inum(ix,iy).gt.0) then
              do inv=1,nv
                 do ispec=1,nspec
                    linemap_thin(ix,iy,inv,ispec)=linemap_thin(ix,iy,inv,ispec)/dble(Inum(ix,iy))
                    linemap_thick(ix,iy,inv,ispec)=linemap_thick(ix,iy,inv,ispec)/dble(Inum(ix,iy))
                 end do
              end do
           end if
        end do
     end do


     if (myid.eq.0) then

        do ispec=1,nspec

           open(unit=15,file='linemap_thin_'//trim(molname(ispec))//'_'//trim(adjustl(ctheta))//'.unf',form='unformatted')
           write(15) nxmap,nxmap,nv
           write(15) rmaxm/autors,pixsize/autors
           write(15) vnu
           write(15) linemap_thin(:,:,:,ispec)
           close(15)

           open(unit=15,file='linemap_thick_'//trim(molname(ispec))//'_'//trim(adjustl(ctheta))//'.unf',form='unformatted')
           write(15) nxmap,nxmap,nv
           write(15) rmaxm/autors,pixsize/autors
           write(15) vnu
           write(15) linemap_thick(:,:,:,ispec)
           close(15)

        end do
        
     end if

  end do

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  return

end subroutine linemap
