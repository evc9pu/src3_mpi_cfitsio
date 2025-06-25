subroutine depletion_co(ir,it,ip,cogas)

  use tts_mod
  use grid_mod
  use constants

  implicit none

  real(8) :: S0,Rdg,crosssec,mu,Eco,abunddep,nu,Agr,Nb,crrate
  real(8) :: co_on,co_off,gas_crit,tau_on,tau_off1,tau_off2,tau_offcr
  real(8) :: dens,temp,nh_t,t_t,cogas,yr,gas_t

  integer :: ir,it,ip

  logical :: isnan

  ! input parameters
  S0=1.d0 ! sticking coefficient
  Rdg=1.d-12 ! dust to H number density ratio
  crosssec=1.d-12 ! cross section of the dust grain
  mu=28.d0 ! CO molecule weight
  Eco=1100.d0 ! binding energy in K
  abunddep=1.d-4 ! CO abundance (to calculate saturated CO desorption rate)
  nu=7.d26 ! thermal desorption related
  Agr=1.26d-9 ! thermal desorption related
  Nb=1.26d-9 ! number of CO binding positions on dust grain
  crrate=3.d-17 ! CR rate
  
  yr=365.25d0*24.d0*3600.d0

  ! some parameters
  co_on=S0*Rdg*crosssec*sqrt(8.*k/pi/mH)
  co_off=Agr*Rdg
  gas_crit=1.-Nb*Rdg/abunddep

  ! density and temperature
  dens=(densarr(ir,it,ip,1)+densarr(ir,it,ip,2)+densarr(ir,it,ip,3)+densarr(ir,it,ip,4))/2.34d-24
  temp=tdave(ir,it,ip)*11605.d0

  T_t=temp
  nH_t=dens

  gas_t=1.d0

  if (dens.gt.0.d0) then

     tau_on=1./(co_on*nH_t*sqrt(T_t/mu))/yr
     tau_off1=1./(co_off*nu*exp(-Eco/T_t)/abunddep)/yr
     tau_off2=1./(co_off*nu*exp(-Eco/T_t)/Nb/Rdg)/yr
     tau_offcr=3.3e6/(crrate/1.e-17)*exp(Eco/70.)/yr
     gas_t=(tau_on/tau_off1)*(tau_offcr+tau_off1)/(tau_offcr+tau_on)
     if (gas_t.gt.gas_crit) then
        gas_t=(tau_on*tau_off2+tau_on*tau_offcr)/(tau_on*tau_off2+tau_on*tau_offcr+tau_offcr*tau_off2)
        if (gas_t.lt.gas_crit) then
           print*,'could not find proper gas_t in eq.'
           stop
        end if
     end if
     if (T_t.lt.10.d0) then
        gas_t=tau_on/(tau_on+tau_offcr)
     end if

  endif

  cogas=gas_t

  if (isnan(cogas)) then
     print*,'cogas is nan'
     print*,ir,it,ip
     print*,dens,T_t
     print*,tau_on,tau_off1,tau_off2,tau_offcr
     stop
  end if

end subroutine depletion_co
