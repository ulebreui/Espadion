program espadion
 
  use ionisation_commons
  implicit none

  !================================================================
  ! This routine is the main routine of espadion. It is just 
  ! meant as an example to learn how to use the modules of the code
  ! It can also be use for exploration.
  !================================================================


  ! General parameters linked to ionisation and variables

  real(dp):: mu_ions        = 25.0d0    ! Mean ion mass (in units of mh)
  real(dp):: stickeff_el    = 0.5d0     ! Sticking coefficient of electrions
  real(dp):: zeta_ionis     = 5d-17     ! CR Ionisation rate
  real(dp):: epsilon_ionis  = 1d-6      ! Tolerance of the ionisation scheme
  integer :: nitermax_ionis = 1000      ! Maximum number of iterations
  real(dp) :: psi0          = -2.5d0    ! Initial guess for psi
  real(dp) :: epsone        = 0.99999d0 ! Maximal electron fraction /!\ 1 is a singularity
  real(dp) :: psi_new       = -2.5d0    ! new guess for psi
 
  real(dp), dimension(:), allocatable      :: ni
  real(dp), dimension(:), allocatable      :: ne
  real(dp), dimension(:,:), allocatable    :: Zd
  real(dp), dimension(:), allocatable      :: Zd_loc

  real(dp), dimension(:), allocatable      :: psi_old


  ! Gas related parameters and variables

  real(dp) :: nHmin     = 1e4    ! Min density
  real(dp) :: nHmax     = 1e15   ! Max density
  integer  :: nnH       = 100    ! Number of density bins
  real(dp) :: Temp      = 10.0d0 ! 10 Kelvin
  real(dp) :: mu_gas    = 2.31d0
  integer  :: inH
  real(dp) :: zeta_nh

  real(dp), dimension(:), allocatable      :: nH  ! Gas number density
  real(dp), dimension(:), allocatable      :: rho ! Gas density

  ! Dust related parameters and variables
  real(dp) :: rhograin   = 2.3d0         ! Grain intrinsic density in g/cm3
  real(dp) :: smin       = 5d-7          ! Grain min size in cm
  real(dp) :: smax       = 250d-7        ! Grain max size in cm
  real(dp) :: ice_mantle = 8.7d-7        ! Ice mantle size -> it changes the grain size but not its mass 
  real(dp) :: dust2gas   = 0.01d0        ! Dust-to-gas ratio
  real(dp) :: mrn        = -3.5d0        ! MRN dust distribution power law

  integer  :: ndust      = 10            ! Number of dust sizes
  integer  :: idust
  real(dp) :: zeta_d,eta_d,m_min

  real(dp), dimension(:), allocatable      :: aminus, aplus, mminus, mplus, adust, mdust
  real(dp), dimension(:,:), allocatable    :: rhodust, n_dust
  real(dp), dimension(:), allocatable      :: ndust_loc


  ! Resistivity parameters and variables

  real(dp), dimension(:), allocatable      :: eta_a
  real(dp), dimension(:), allocatable      :: eta_o
  real(dp), dimension(:), allocatable      :: eta_h
  real(dp), dimension(:), allocatable      :: sigma_o
  real(dp), dimension(:), allocatable      :: sigma_p
  real(dp), dimension(:), allocatable      :: sigma_h

 


  print *, "Entering Espadion, let's get ionised !"


  ! Initialisation of the gas properties
  allocate(nH(1:nnH))
  allocate(rho(1:nnH))
  zeta_nh = (nHmax/nHmin)**(1.0d0/nnH)
  do inH =1,nnH
    nH(inh)  = nHmin  *zeta_nh ** (inH)
    rho(inh) = nH(inh)*mu_gas*mH
  end do

  ! Initialisation of the dust size distribution
  allocate(aplus(1:ndust))
  allocate(aminus(1:ndust))
  allocate(mplus(1:ndust))
  allocate(mminus(1:ndust))
  allocate(adust(1:ndust))
  allocate(mdust(1:ndust))
  allocate(rhodust(1:nnH,1:ndust))
  allocate(n_dust(1:nnH,1:ndust))

  allocate(ni(1:nnH))
  allocate(ne(1:nnH))
  allocate(Zd(1:nnH,1:ndust))
  allocate(Zd_loc(1:ndust))

  zeta_d = (smax/smin)**(1.0d0/ndust)
  eta_d  = zeta_d**3.
  m_min = 4./3.*pi*rhograin*smin**3.
  
  do idust =1,ndust
    ! Bin borders
    aminus(idust) = smin  * zeta_d ** (idust-1)
    aplus (idust) = smin  * zeta_d ** (idust)
    mminus(idust) = m_min * eta_d  ** (idust-1)
    mplus(idust)  = m_min * eta_d  ** (idust)

    ! Bin centers
    adust (idust) = sqrt(aminus(idust)   * aplus(idust))
    mdust (idust) = 4./3. *pi *rhograin  * adust(idust)**3
  enddo

  ! Initialise the MRN distribution
  do inH=1,nnH
    do idust =1,ndust
      rhodust(inH,idust)= dust2gas*rho(inh)*(aplus(idust)**(4.0+mrn)-aminus(idust)**(4.0+mrn))/(aplus(ndust)**(4.0+mrn)-aminus(1)**(4.0+mrn))
      n_dust(inH,idust) = rhodust(inH,idust)/mdust(idust)
    end do
  end do


  call init_PSI(psi0,mu_ions,stickeff_el,epsone,epsilon_ionis)

  ! Actual charge computation
  do inH=1,nnH
    call compute_charges(ni(inH),ne(inH),zd(inH,:),nH(inh), n_dust(inh,:), Temp, adust, zeta_ionis, mu_ions,stickeff_el, psi0, psi_new, ndust,epsone,nitermax_ionis,epsilon_ionis)
    psi0 = psi_new !change the guess to something as close as possible from the solution ! In ramses take the old value for e.g. or maybe we tabulate it ?
  end do

  call output(nh,ni,ne,zd,nnh,ndust)


end program espadion


subroutine init_PSI(psi0,mu_ions,stickeff_el,epsone,epsilon_ionis)
  use ionisation_commons
  implicit none

  real(dp) :: psi0,mu_ions,epsone,stickeff_el,epsilon_ionis
  real(dp) :: psi_loc, fpsi, dfdpsi,convergence_ionis, thetai

  integer :: niter_ionis

  convergence_ionis = 1d4

  niter_ionis       = 0

  !Initial guess

  thetai=stickeff_el*dsqrt(mu_ions*mH/m_el)
  do while(convergence_ionis>epsilon_ionis.and.niter_ionis<1000)
     psi_loc=psi0
     !f(psi)
     fpsi=(1.d0-psi_loc)/(epsone*thetai)*exp(-psi_loc)-1.0d0
     ! df/dpsi
     dfdpsi=-exp(-psi_loc)/(epsone*thetai)-(1.d0-psi_loc)/(epsone*thetai)*exp(-psi_loc)
     !Newton-Raphson step
     psi_loc=psi_loc-fpsi/dfdpsi
     !Convergence check
     convergence_ionis=abs(fpsi)
     niter_ionis=niter_ionis+1
     psi0=psi_loc
  end do

end subroutine init_PSI


subroutine compute_charges(ni,ne,zd,nH_loc, n_k, T, l_grain_loc, zeta_ionis, mu_ions,stickeff_el, psi0, psi_new, ndust,epsone,nitermax_ionis,epsilon_ionis)
  use ionisation_commons
  implicit none

  integer :: ndust,idust

  real(dp) :: ni,ne
  real(dp),dimension(1:ndust) :: zd
  real(dp) :: nH_loc, T,psi0,mu_ions,thetai,stickeff_el,zeta_ionis,psi_new
  real(dp),dimension(1:ndust) :: l_grain_loc, n_k,tau_k,alpha_k

  real (dp) :: epsone,psi_loc,sigmav_ie,vi
  real (dp) :: convergence_ionis,epsilon_ionis

  real(dp) ::fpsi, eps_psi, eps_theta,n_i_loc,dfdpsi,dni_dpsi
  integer :: niter_ionis,nitermax_ionis,lowT

  !Limit for epsilon cannot be 1.
  epsone     = 0.99999d0

  convergence_ionis = 1d4
  niter_ionis       = 0

  thetai=stickeff_el*dsqrt(mu_ions*mH/m_el)
  ! Grain size in cm
  ! Ion electron collisional cross-section
  sigmav_ie=(2d-7)/dsqrt(T/300.0d0)
  ! Reduced dust temperature
  tau_k = l_grain_loc*kB*T/e_el_stat**2.
  ! Ion thermal velocity
  vi=dsqrt(8.0d0*kB*T/(pi*mu_ions*mH))
  ! Variable alpha_k for convenience (used in Marchand et al., 21)
  alpha_k=dsqrt(8.0d0/(pi*tau_k))

  !Psi0 determined, serves as a guess for psi unless it's not the first dt
  convergence_ionis=1d4
  niter_ionis=0
  psi_loc=psi0

  lowT=0

  do while(convergence_ionis>epsilon_ionis .and.niter_ionis<nitermax_ionis)

      ! Electron fraction
      eps_psi  = (1.0d0-psi_loc)/thetai*exp(-psi_loc)

      ! Variable added for convenience. Should in principle be replaced everywhere but the expression are already horrible
      eps_theta= eps_psi*thetai
      ! Ion number density
      n_i_loc  = sum(-(psi_loc*tau_k(:)+(1.0d0-eps_theta**2.0)/(1.0d0+eps_theta*alpha_k(:)+eps_theta**2.0))*n_k(:)/(1.0d0-eps_psi))

      !f(psi). This is the function we try to zero.
      fpsi=sigmav_ie*eps_psi*n_i_loc**2.0/(zeta_ionis*nH_loc)+n_i_loc*vi/(zeta_ionis*nH_loc)*(sum(n_k(:)*pi*l_grain_loc(:)**2.0*((1.0d0-psi_loc)+(2.0d0/tau_k(:))*(eps_theta**2.0+eps_theta)/(1.0d0+eps_theta*alpha_k(:)+eps_theta**2.0))))-1.0d0


      ! dni_dpsi /!\ trust me it's correct :). Yes. I suffered
      dni_dpsi   = -1.0d0/(1.0d0-eps_psi)**2.*(-eps_psi*(2.0d0-psi_loc)/(1.0d0-psi_loc))&
           &*(sum(n_k(:)*(-eps_theta**4.0 + (4.0d0*thetai**2.-thetai**3.*alpha_k(:))*eps_psi**2.&
           &+(2.0d0*thetai*alpha_K(:)-4.0d0*thetai**2.0)*eps_psi +1.0d0 - thetai*alpha_k(:))/(1.0d0 + eps_theta*alpha_k(:) + eps_theta**2.0d0)**2.0))&
           &-(psi_loc/(1.0d0-eps_psi)*(-eps_psi*(2.0d0-psi_loc)/(1.0d0-psi_loc))+1.0d0)*sum(n_k(:)*tau_k(:))/(1.0d0-eps_psi)

      ! Actual derivative
      dfdpsi=sigmav_ie*n_i_loc/(zeta_ionis*nH_loc)*(n_i_loc*(-eps_psi*(2.0d0-psi_loc)/(1.0d0-psi_loc))+2.0d0*eps_psi*dni_dpsi)+2.0d0*n_i_loc*vi/&
           &(zeta_ionis*nH_loc)*(eps_theta**2.+eps_theta)*((dni_dpsi/n_i_loc+(-eps_psi*(2.0d0-psi_loc)/(1.0d0-psi_loc))*(2.0d0*eps_theta+1.0d0)/(eps_psi*(1.0+eps_theta)))&
           &*(sum(n_k(:)*pi*l_grain_loc(:)**2.0/(tau_k(:)*(1.0d0+eps_theta*alpha_k(:)+eps_theta**2.0d0))))&
           &+((-eps_psi*(2.0d0-psi_loc)/(1.0d0-psi_loc))*sum(n_k(:)*pi*l_grain_loc(:)**2.0*(2.0d0*eps_psi*thetai**2.+thetai*alpha_k(:))&
           &/(tau_k(:)*(1.0d0+eps_theta*alpha_k(:)+eps_theta**2.0d0)**2.0d0))))+n_i_loc*vi/(zeta_ionis*nH_loc)*(dni_dpsi/n_i_loc*(1.0d0-psi_loc)-1.0d0)*(sum(n_k(:)*pi*l_grain_loc(:)**2.0))

      !Psi^n+1=Psi^n-f(psi)/dfdpsi(psi)
      psi_loc=psi_loc-fpsi/dfdpsi
      if(psi_loc<psi0) then
         psi_loc = psi0*0.9999999
         lowT=1              
         if(convergence_ionis==1d5) then
            convergence_ionis=epsilon_ionis*0.1 !To not loop indefinitely
         else
            convergence_ionis=1d5  !We cheated so we shouldn't converge at this step
         endif 
      else if (psi_loc>0.0d0) then
         psi_loc=-1d-5
         convergence_ionis=1d4 !We cheated so we shouldn't converge at this step
      end if
      convergence_ionis=abs(fpsi)
      niter_ionis=niter_ionis+1
  end do

 eps_psi=(1.0d0-psi_loc)/thetai*exp(-psi_loc)
 if(lowT>0)eps_psi=epsone
 ni=0.0d0
 do idust=1,ndust
    zd(idust)=psi_loc*tau_k(idust)+(1.0d0-eps_psi**2.0*thetai**2.0)/(1.0d0+eps_psi*thetai*alpha_k(idust)+eps_psi**2*thetai**2.0)
    ni= ni-zd(idust)*n_k(idust)/(1.0d0-eps_psi)
    if(lowT>0) ni=sqrt(eps_psi*zeta_ionis*nH_loc/sigmav_ie)
 end do
 ne=eps_psi*ni
 psi_new=psi_loc
 if(lowT>0)then
    psi_new=psi0
 end if

end subroutine compute_charges

!   as_He_dust = 1.28
!   as_He_ions = 1.14
!   as_He_el   = 1.16


!   real(dp) :: B_field

!      ! Resistivity computation

!      if (res_Marchand) then

!         !B_gauss = min(B_0_lee*sqrt(u_prim(i,irho)*unit_nh/1d4),B_threshold)! Magnetic field
!         B_gauss = B_0_lee*sqrt(u_prim(i,irho)*unit_nh/1d4)! Magnetic field


!         sigmav_dust(:)=pi*(l_grain_loc(:))**2.*dsqrt(8.0*kB*T/(pi*2.0d0*mH))*(1.0d0+dsqrt(pi/(2.*tau_k(:))))

!         t_sdust(:)=dsqrt(pi*gamma/8.0d0)*(rhograin)*(sdust(i,:)*unit_l)/(u_prim(i,irho)*unit_d*cs_eos(T)*unit_v)

!         mu_i=2.0d0*mH*mu_ions*mH/(2.0d0*mH+mu_ions*mH)
!         mu_e=2.0d0*mH*m_el/(m_el+2.0d0*mH)

!         vrms_i =dsqrt(8.0d0*kB*T/(pi*mu_i))*1d-5  ! /!\ /!\ These velocities need to be in km/s for the Pinto & Galli 2008 fit
!         vrms_el=dsqrt(8.0d0*kB*T/(pi*mu_e))*1d-5

!         sigmav_el  = 3.16d-11*vrms_el**1.3
!         sigmav_ions= 2.4d-9  *vrms_i**0.6

!         t_sel  = 1.0d0/as_He_el*((m_el+2.0d0*mH)/(2.0d0*mH))/sigmav_el/nH_loc
!         t_sions= 1.0d0/as_He_ions*((mu_ions*mH+2.0d0*mH)/(2.0d0*mH))/sigmav_ions/nH_loc

!         sigmas_el     = (ne(i))*e_el_stat**2.*t_sel/m_el
!         sigmas_ions   = (ni(i))*e_el_stat**2*t_sions/(mu_ions*mH)
!         sigmas_dust(:)= n_k(:)*(zd(i,:)*e_el_stat)**2.*t_sdust(:)/(mdust(i,:)*unit_m)

!         omegas_el     = -e_el_stat*B_gauss/clight/m_el
!         omegas_ions   = e_el_stat *B_gauss/clight/(mu_ions*mH)
!         omegas_dust   = zd(i,:)*e_el_stat*B_gauss/clight/(mdust(i,:)*unit_m)

!         ! Conductivities
!         sigma_o(i)=sum(sigmas_dust(:))+sigmas_ions+sigmas_el
!         sigma_p(i)=sigmas_ions/(1.0d0+(omegas_ions*t_sions)**2.0)+sigmas_el/(1.0d0+(omegas_el*t_sel)**2.0)+sum(sigmas_dust(:)/(1.0d0+(omegas_dust(:)*t_sdust(:))**2.0))
!         sigma_H(i)=-sigmas_ions*(omegas_ions*t_sions)/(1.0d0+(omegas_ions*t_sions)**2.0)-sigmas_el*(omegas_el*t_sel)/(1.0d0+(omegas_el*t_sel)**2.0)-sum(sigmas_dust(:)*(omegas_dust(:)*t_sdust(:))/(1.0d0+(omegas_dust(:)*t_sdust(:))**2.0))

!         ! Resistivities
!         eta_o(i)=1.0d0/sigma_o(i)
!         eta_H(i)=sigma_H(i)/(sigma_p(i)*2.+sigma_h(i)**2.)
!         eta_a(i)=sigma_p(i)/(sigma_p(i)**2.+sigma_H(i)**2.)-1.0d0/sigma_o(i)
!         !Hall factors
!         do idust=1,ndust
!            gamma_d(i,idust)=t_sdust(idust)*omegas_dust(idust)
!         end do
!       end if


!      !endif
! end do


! end subroutine compute_charges