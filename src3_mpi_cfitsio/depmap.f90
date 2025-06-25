subroutine depmap()

  use grid_mod  
  use dust_mod
  use tauint_mod
  use taunum_mod
  use tts_mod
  use constants
  use tts_mod
  use opacin_mod
  use depint_mod
  use output_mod
  use ttsre_mpi_mod

  implicit none

  integer :: idep,peelid,ir0,nr
  integer :: nrmod,nrpro,irstart,irend,nthet
  integer :: iir,iit,ix,iy,ir,it,ip,ir1,it1,ip1
  
  real(8) :: rmaxm,pixsize,xmap,ymap
  real(8) :: peelthet,sint,cost,thet,rad
  real(8) :: radmin,maxfd,minfd
  real(8),pointer :: fdmap(:,:,:)

  character(len=10) :: ctheta

  logical :: ifout

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if (ifprint) print*,'calulate depletion maps now'

  ndep=1

  allocate(S0_arr(ndep),Rdg_arr(ndep),crosssec_arr(ndep))
  allocate(mu_arr(ndep),Eco_arr(ndep),abunddep_arr(ndep))
  allocate(nu_arr(ndep),Agr_arr(ndep),Nb_arr(ndep),crrate_arr(ndep))
  allocate(molname(ndep))
  allocate(fdmap(nxhst,nxhst,ndep))
  allocate(totmap(nxhst,nxhst,ndep),gasmap(nxhst,nxhst,ndep))
  allocate(totcol(ndep),gascol(ndep))

  totmap=0.d0
  gasmap=0.d0
  fdmap=0.d0

  idep=1 ! CO, values form Visser et al. 2009
  molname(idep)='CO'
  S0_arr(idep)=1.d0
  Rdg_arr(idep)=1.d-12
  crosssec_arr(idep)=3.14d-10
  mu_arr(idep)=28.d0
  Eco_arr(idep)=1100.d0 !850K for CO to CO ice surface, 1100K for CO to water ice surface
  nu_arr(idep)=7.d26
  Agr_arr(idep)=1.26d-9
  Nb_arr(idep)=1.d6
  abunddep_arr(idep)=1.d-4
  crrate_arr(idep)=3.d-17

  rmaxm=1.2*rcore*autors
  pixsize=2.d0*rmaxm/dble(nxhst)

  do peelid=1,npeel

     peelthet=thete_arr(peelid)

     write(ctheta,'(F5.1)') peelthet*rad2deg

     radmin=1.d0/sin(peelthet)
     call locate(ravearr,nrg-1,radmin,ir0)
     ir0=ir0+1
     if (ifprint) print*,'thet,ir0',peelthet*180.d0/pi,radmin,ir0

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

!     print*,myid,irstart,irend

     do iir=irstart,irend

!        print*,'myid,iir',myid,iir

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

           call deppeel(ir1,it1,ip1,ifout)

           ux=sin(peelthet)
           uy=0.d0
           uz=cos(peelthet)

           ifout=.true.

           totcol=0.d0
           gascol=0.d0

           call deppeel(ir1,it1,ip1,ifout)

           ymap=zp*sin(peelthet)-xp*cos(peelthet)
           xmap=yp

           ix=int((xmap+rmaxm)/pixsize)+1
           iy=int((ymap+rmaxm)/pixsize)+1

           if (ix.ge.1.and.ix.le.nxhst.and.iy.ge.1.and.iy.le.nxhst) then
              totmap(ix,iy,:)=totmap(ix,iy,:)+totcol(:)
              gasmap(ix,iy,:)=gasmap(ix,iy,:)+gascol(:)
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

           call deppeel(ir1,it1,ip1,ifout)

           ux=sin(peelthet)
           uy=0.d0
           uz=cos(peelthet)

           ifout=.true.

           totcol=0.d0
           gascol=0.d0

           call deppeel(ir1,it1,ip1,ifout)

           ymap=zp*sin(peelthet)-xp*cos(peelthet)
           xmap=yp

           ix=int((xmap+rmaxm)/pixsize)+1
           iy=int((ymap+rmaxm)/pixsize)+1

           if (ix.ge.1.and.ix.le.nxhst.and.iy.ge.1.and.iy.le.nxhst) then
              totmap(ix,iy,:)=totmap(ix,iy,:)+totcol(:)
              gasmap(ix,iy,:)=gasmap(ix,iy,:)+gascol(:)
           end if

!========================================
        end do

     end do

     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

     icountbuf=nxhst*nxhst*ndep
     call MPI_ALLREDUCE(MPI_IN_PLACE,totmap,icountbuf,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE,gasmap,icountbuf,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

     if (myid.eq.0) then

        maxfd=0.d0
        minfd=10000.d0
        do idep=1,ndep
           do ix=1,nxhst
              do iy=1,nxhst
                 if (totmap(ix,iy,idep).gt.0) then
                    fdmap(ix,iy,idep)=totmap(ix,iy,idep)/gasmap(ix,iy,idep)
                    if (fdmap(ix,iy,idep).gt.maxfd) then
                       maxfd=fdmap(ix,iy,idep)
                    end if
                    if (fdmap(ix,iy,idep).lt.minfd) then 
                       minfd=fdmap(ix,iy,idep)
                    end if
                 end if
              end do
           end do
        end do
        print*,'maxfd,minfd',maxfd,minfd

        do idep=1,ndep
           
           open(unit=15,file='depmap_'//trim(molname(idep))//'_'//trim(adjustl(ctheta))//'.unf',form='unformatted')
           write(15) nxhst,nxhst
           write(15) rmaxm/autors,pixsize/autors
           write(15) fdmap(:,:,idep)
           write(15) totmap(:,:,idep)
           write(15) gasmap(:,:,idep)
           close(15)

        end do
        
     end if

  end do

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  return
  
end subroutine depmap
  
