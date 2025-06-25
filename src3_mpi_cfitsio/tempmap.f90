subroutine tempmap() 

  use ttsre_mpi_mod
  use grid_mod
  use dust_mod
  use tauint_mod
  use taunum_mod
  use constants
  use tts_mod
  use opacin_mod
  use tempmap_mod
  use output_mod

  implicit none

  real(8) :: peelthet,radmin,rmaxm,rad,thet,sint,cost
  real(8) :: pixsize,xmap,ymap,ncrit

  integer :: peelid,ir0,nr,iir,ir,nthet,iit,it,ip
  integer :: ir1,it1,ip1,ix,iy,irl
  integer :: nrmod,nrpro,irstart,irend

  character(len=4) :: suffix

  logical :: ifout

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if (ifprint) print*,'doing temp now'

  allocate(Tmap(nxhst,nxhst),rhomap(nxhst,nxhst),vmap(nxhst,nxhst))

!  ncrit=1.d3
  ncrit=0.d0
  
  rmaxm=20000.d0*autors
  pixsize=2.d0*rmaxm/dble(nxhst)

  do peelid=1,npeel

     peelthet=thete_arr(peelid)
     
     write(suffix,'("_",I3.3)') int(peelthet*180.d0/pi)
  
     radmin=1.d0/sin(peelthet)
     call locate(ravearr,nrg-1,radmin,ir0)
     ir0=ir0+1
     if (ifprint) print*,'thet,ir0',peelthet*180.d0/pi,radmin,ir0
  
     Tmap=0.d0
     rhomap=0.d0
     vmap=0.d0

     nr=int(rmaxm*2.d0/sin(peelthet)/autors/100.d0)

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

     do iir=irstart,irend
        
        rad=100.d0*dble(iir)*autors
     
        call locate(rarr,nrg,rad,ir)
     
        if (rarr(ir).gt.rad) ir=ir-1
        
        print*,'iir,ir',iir,ir,rad/autors
        
        nthet=int(180.d0/1.d0)

        do iit=1,nthet-1
           
           thet=pi/180.d0*dble(iit)*1.d0
        
           call locate(thetarr,ntg,thet,it)
        
           if (thetarr(it).gt.thet) it=it-1

           ip=1
        
!=====================================

           do irl=1,2
        
              xp=0.d0
              sint=sin(thet)
              cost=cos(thet)
              if (irl.eq.1) yp=rad*sint
              if (irl.eq.2) yp=-rad*sint
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
              
              call tempeel(ir1,it1,ip1,ncrit,ifout)
              
              ux=sin(peelthet)
              uy=0.d0
              uz=cos(peelthet)
              
              ifout=.true.
              
              Trhor=0.d0
              rhor=0.d0
              vrhor=0.d0
              
              call tempeel(ir1,it1,ip1,ncrit,ifout)
              
              ymap=zp*sin(peelthet)-xp*cos(peelthet)
              xmap=yp
              
              ix=int((xmap+rmaxm)/pixsize)+1
              iy=int((ymap+rmaxm)/pixsize)+1

              if (ix.ge.1.and.ix.le.nxhst.and.iy.ge.1.and.iy.le.nxhst) then
                 Tmap(ix,iy)=Tmap(ix,iy)+Trhor
                 vmap(ix,iy)=vmap(ix,iy)+vrhor
                 rhomap(ix,iy)=rhomap(ix,iy)+rhor
              end if
              
           end do
           !=====================================
        end do
        
     end do

     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

     icountbuf=nxhst*nxhst
     call MPI_ALLREDUCE(MPI_IN_PLACE,Tmap,icountbuf,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE,vmap,icountbuf,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE,rhomap,icountbuf,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    
     do ix=1,nxhst
        do iy=1,nxhst
           if (rhomap(ix,iy).gt.0.d0) then
              Tmap(ix,iy)=Tmap(ix,iy)/rhomap(ix,iy)
              vmap(ix,iy)=vmap(ix,iy)/rhomap(ix,iy)
           else
              Tmap(ix,iy)=0.d0
           end if
        end do
     end do
     
     if (myid.eq.0) then 

        open(unit=15,file='tempmap'//suffix//'.unf',form='unformatted')
        write(15) nxhst,nxhst
        write(15) rmaxm/autors,pixsize/autors
        write(15) Tmap
        close(15)
        
        open(unit=15,file='vmap'//suffix//'.unf',form='unformatted')
        write(15) nxhst,nxhst
        write(15) rmaxm/autors,pixsize/autors
        write(15) vmap
        close(15)

     end if

  end do

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  
  return
  
end subroutine tempmap
