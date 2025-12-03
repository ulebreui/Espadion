program espadion
 
  use ionisation_commons
  use ionisation_functions
  implicit none

  !================================================================
  ! This routine is the main routine of espadion. It is just 
  ! meant as an example to learn how to use the modules of the code
  ! It can also be use for exploration.
  !================================================================
  character(len=70):: nmlfile,infile
  integer :: io_n
  logical::nml_ok

  real(dp) :: barotrop

  real(dp) :: B_vs_rho

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
  real(dp) :: nHmax     = 1e12   ! Max density
  integer  :: nnH       = 100    ! Number of density bins
  real(dp) :: Temp      = 10.0d0 ! 10 Kelvin
  real(dp) :: mu_gas    = 2.31d0
  real(dp) :: B_0       = 3d-5   ! Value of the B field at 10^4

  integer  :: inH
  real(dp) :: zeta_nh

  real(dp), dimension(:), allocatable      :: nH  ! Gas number density
  real(dp), dimension(:), allocatable      :: rho ! Gas density

  real(dp), dimension(:), allocatable      :: B_field  ! Magnetic field

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

  real(dp), dimension(:,:), allocatable      :: tau_js
  real(dp), dimension(:,:), allocatable      :: Omega_js_B
  real(dp), dimension(:,:), allocatable      :: Sigma_js


   namelist/params/mu_ions,stickeff_el,zeta_ionis,epsilon_ionis,nitermax_ionis&
   &,psi0,epsone,nHmin,nHmax,nnh,temp,mu_gas,B_0,rhograin,smin,smax,ice_mantle,dust2gas,mrn,ndust

   print *, "Entering Espadion, let's get ionised !"
   print *, "########################################################################################################################################"
   print *, "########################################################################################################################################"
   CALL getarg(1,infile)
   nmlfile=TRIM(infile)
   INQUIRE(file=infile,exist=nml_ok)

   if(.not. nml_ok)then
      write(*,*)'File '//TRIM(infile)//' does not exist'
      stop
   end if
   print *, "Namelist reading : ", nmlfile
   
   open(13,file=nmlfile)
   print *, "Namelist reading  !"
   read(13,params,IOSTAT=io_n)
   rewind(13)
   if (io_n/=0) then
      write(*,*) 'Invalid line in the namelist'
      stop
   end if

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

  allocate(tau_js(1:nnH,1:ndust+2))
  allocate(Omega_js_B(1:nnH,1:ndust+2))
  allocate(Sigma_js(1:nnH,1:ndust+2))

  allocate(sigma_o(1:nnH))
  allocate(sigma_p(1:nnH))
  allocate(sigma_h(1:nnH))

  allocate(eta_o(1:nnH))
  allocate(eta_h(1:nnH))
  allocate(eta_a(1:nnH))


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
    adust (idust)  = sqrt(aminus(idust)   * aplus(idust))
    mdust (idust)  = 4./3. *pi *rhograin  * adust(idust)**3
    adust (idust) =  adust (idust) + ice_mantle
  enddo


  ! Initialise the MRN distribution
  do inH=1,nnH
    do idust =1,ndust
      rhodust(inH,idust)= dust2gas*rho(inh)*(aplus(idust)**(4.0+mrn)-aminus(idust)**(4.0+mrn))/(aplus(ndust)**(4.0+mrn)-aminus(1)**(4.0+mrn))
      n_dust(inH,idust) = rhodust(inH,idust)/mdust(idust)
    end do
  end do

  
  ! First guess for psi. To put in initial conditions. PSI0 should be a global variable
  call init_PSI(psi0,mu_ions,stickeff_el,epsone,epsilon_ionis)

  ! Actual charge/resistivity computation
  do inH=1,nnH
    call compute_charges(ni(inH),ne(inH),zd(inH,:),nH(inh), n_dust(inh,:), barotrop(Temp,nH(inh)), adust, zeta_ionis, mu_ions,stickeff_el, psi0, psi_new, ndust,epsone,nitermax_ionis,epsilon_ionis)
    psi0 = psi_new !change the guess to something as close as possible from the solution ! In ramses take the old value for e.g. or maybe we tabulate it ?
    
    ! Now we compute the quantities needed for the resistivities

    ! For the dust
    do idust =1,ndust
      tau_js(inH,idust)     =  tau_d(n_dust(inh,idust),adust(idust),rhograin,nH(inH),barotrop(Temp,nH(inh)),mu_gas*mH)
      Sigma_js(inH,idust)   =  sigma_j(n_dust(inh,idust),zd(inh,idust),tau_js(inH,idust) ,mdust(idust))
      Omega_js_B(inh,idust) =  Omega_j_over_B(mdust(idust),zd(inh,idust))
    end do

    ! Now the ions
    tau_js(inH,ndust+1)     = tau_i(nH(inh),barotrop(Temp,nH(inh)),mu_ions)
    Sigma_js(inH,ndust+1)   = sigma_j(ni(inh),1.0d0,tau_js(inH,ndust+1),mu_ions*mH)
    Omega_js_B(inh,ndust+1) =  Omega_j_over_B(mu_ions*mH,1.0d0)


    ! Now the electron
    tau_js(inH,ndust+2)     = tau_el(nH(inh),barotrop(Temp,nH(inh)))
    Sigma_js(inH,ndust+2)   = sigma_j(ne(inh),-1.0d0,tau_js(inH,ndust+2),m_el)
    Omega_js_B(inh,ndust+2) = Omega_j_over_B(m_el,-1.0d0)


    sigma_o(inH) = 0.0d0
    sigma_P(inH) = 0.0d0
    sigma_H(inH) = 0.0d0

    do idust =1,ndust+2
      sigma_o(inH) = sigma_o(inH) + Sigma_js(inH,idust)
      sigma_P(inH) = sigma_P(inH) + SigmaP(Sigma_js(inH,idust),Omega_js_B(inh,idust),tau_js(inH,idust),B_vs_rho(B_0,nh(inh)))
      sigma_H(inH) = sigma_H(inH) + SigmaH(Sigma_js(inH,idust),Omega_js_B(inh,idust),tau_js(inH,idust),B_vs_rho(B_0,nh(inh)))
    end do

    eta_o(inH)=1.0d0/sigma_o(inH)
    eta_H(inH)=sigma_H(inH)/(sigma_p(inH)*2.+sigma_h(inH)**2.)
    eta_a(inH)=sigma_p(inH)/(sigma_p(inH)**2.+sigma_H(inH)**2.)-1.0d0/sigma_o(inH)

  end do


  call output(nh,ni,ne,zd,tau_js,Sigma_js,Omega_js_B,sigma_o,sigma_p,sigma_h,eta_o,eta_H,eta_a,nnh,ndust)


end program espadion

double precision function barotrop(Temp,nH)
   use ionisation_commons
   implicit none

   real(dp) :: nH,n_star,n1,n2,n3 , n_star2,gammaloc,gamma2,temp
   gammaloc= 5./3.
   gamma2=7./5.
   n_star = 1e10
   n_star2=30.*n_star

   barotrop = Temp * ( 1.0d0 + (nH/n_star)**(gammaloc-1.0d0) / (1.+(nH/n_star2)**(gammaloc-gamma2)) )

end function barotrop

double precision function B_vs_rho(B0,nH)
   use ionisation_commons
   implicit none

   real(dp) :: nH,B0

   B_vs_rho = B0*sqrt(nH/1d4)! Magnetic field

end function B_vs_rho





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