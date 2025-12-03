!===========================================================================
subroutine title(n,nchar) ! Shamelessely stolen from RAMSES (Teyssier, 2002)
!===========================================================================
  implicit none
  integer::n
  character(LEN=5)::nchar

  character(LEN=1)::nchar1
  character(LEN=2)::nchar2
  character(LEN=3)::nchar3
  character(LEN=4)::nchar4
  character(LEN=5)::nchar5

  if(n.ge.10000)then
     write(nchar5,'(i5)') n
     nchar = nchar5
  elseif(n.ge.1000)then
     write(nchar4,'(i4)') n
     nchar = '0'//nchar4
  elseif(n.ge.100)then
     write(nchar3,'(i3)') n
     nchar = '00'//nchar3
  elseif(n.ge.10)then
     write(nchar2,'(i2)') n
     nchar = '000'//nchar2
  else
     write(nchar1,'(i1)') n
     nchar = '0000'//nchar1
  endif


end subroutine title


! Writes the simulation output. The routine is a bit ugly and will be rewritten soon
subroutine output(nh,ni,ne,zd,tau_js,Sigma_js,Omega_js_B,sigma_o,sigma_p,sigma_h,eta_o,eta_H,eta_a,nnh,ndust)
  use ionisation_commons
  implicit none
  integer  :: inh,idust,nnh,ndust,ilun
  real(dp), dimension(1:nnh) :: nh,ni,ne,sigma_o,sigma_p,sigma_h,eta_o,eta_H,eta_a
  real(dp), dimension(1:nnh,1:ndust)   :: zd
  real(dp), dimension(1:nnh,1:ndust+2) :: tau_js,Sigma_js,Omega_js_B
  character(LEN = 5) :: nchar
  character(len=80)  :: path, makedirectory,format_out

  format_out=trim("unformatted")
  ilun=20
  !call title(1,nchar)
  path='ion_files'
  makedirectory = 'mkdir ' //  trim(path)
  call system(makedirectory)

  open(ilun,file=trim(path) //trim('/info.dat'))
  write(ilun,'("nnh        =",I11)')nnh
  write(ilun,'("ndust        =",I11)')NDUST


   open(ilun,file=trim(path) //trim('/nh'), form=format_out,access='stream')
   do inh = 1,nnh
      write(ilun) nh(inh)
   end do
   close(ilun)   
   open(ilun,file=trim(path) //trim('/ni'), form=format_out,access='stream')
   do inh = 1,nnh
      write(ilun) ni(inh)
   end do
   close(ilun)
   open(ilun,file=trim(path) //trim('/ne'), form=format_out,access='stream')
   do inh = 1,nnh
      write(ilun) ne(inh)
   end do
   close(ilun)
   open(ilun,file=trim(path) //trim('/zd'), form=format_out,access='stream')
   do idust=1,ndust
      do inh = 1,nnh
          write(ilun) zd(inh,idust)
      end do
   end do


   open(ilun,file=trim(path) //trim('/tau'), form=format_out,access='stream')
   do idust=1,ndust+2
      do inh = 1,nnh
          write(ilun) tau_js(inh,idust)
      end do
   end do
   open(ilun,file=trim(path) //trim('/sigma'), form=format_out,access='stream')
   do idust=1,ndust+2
      do inh = 1,nnh
          write(ilun) Sigma_js(inh,idust)
      end do
   end do

   open(ilun,file=trim(path) //trim('/omega_B'), form=format_out,access='stream')
   do idust=1,ndust+2
      ! All cells have the same valu
      do inh = 1,nnh
      write(ilun) Omega_js_B(inh,idust)
      end do
   end do


   open(ilun,file=trim(path) //trim('/eta_o'), form=format_out,access='stream')
   do inh = 1,nnh
      write(ilun) eta_o(inh)
   end do
   close(ilun)

   open(ilun,file=trim(path) //trim('/eta_h'), form=format_out,access='stream')
   do inh = 1,nnh
      write(ilun) eta_h(inh)
   end do
   close(ilun)

   open(ilun,file=trim(path) //trim('/eta_a'), form=format_out,access='stream')
   do inh = 1,nnh
      write(ilun) eta_a(inh)
   end do
   close(ilun)

   open(ilun,file=trim(path) //trim('/sigma_o'), form=format_out,access='stream')
   do inh = 1,nnh
      write(ilun) sigma_o(inh)
   end do
   close(ilun)

   open(ilun,file=trim(path) //trim('/sigma_h'), form=format_out,access='stream')
   do inh = 1,nnh
      write(ilun) sigma_h(inh)
   end do
   close(ilun)

   open(ilun,file=trim(path) //trim('/sigma_p'), form=format_out,access='stream')
   do inh = 1,nnh
      write(ilun) sigma_p(inh)
   end do
   close(ilun)
 end subroutine output