module phys_const
  use precision
  implicit none

  real(dp),parameter:: e_el_stat     = 4.8d-10      ! Electron charge /!\ in statcoulomb here.
  real(dp),parameter:: m_el          = 9.109d-28    ! Electron mass in grams
  real(dp),parameter:: clight        = 3d10         ! Speed of light
  real(dp),parameter:: mH            = 1.67d-24     ! Mass of an hydrogen atom
  real(dp),parameter:: kB            = 1.38d-16     ! Boltzmann constant
  real(dp),parameter:: Grav          = 6.67d-8      ! Gravity constant
  real(dp),parameter:: au            = 1.5d13       ! Astronomical unit
  real(dp),parameter:: pc            = 3d18         ! Parsec
  real(dp),parameter:: Msun          = 2d33         ! Solar mass
  real(dp),parameter:: pi            = 3.141592653589793238462643383279d0  ! Pi

end module phys_const



module ionisation_commons
  use precision
  use phys_const

end module ionisation_commons


module ionisation_functions

  use ionisation_commons

  contains


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


end module ionisation_functions