Module mod_afssh
!! Hammes-Schiffer, Tully, JCP 101, 4657 (1994)
implicit none
real*8, parameter :: clight=2.99792458D10,av=6.0221367D23,hbar=1.05457266D-34
real*8, parameter :: kb=1.3806503d-23
real*8, parameter :: pascal_to_atm=9.86923267d-6,kcal_to_J=6.95d-21
real*8, parameter :: amu2kg=1.66053892d-27
real*8, parameter :: au2kg=9.10938291d-31,au2J=4.35974417d-18,au2m=5.2917721092d-11,au2s=2.418884326505d-17
real*8, parameter :: q_el=2.307075d-28
complex*16,parameter :: iota = (0.d0,1.d0)
real*8 pi,wave_to_J

!! rate estimates
real*8 lambda,fac,omg1_prime,q0

!! Potential
integer nquant
real*8 V_coup,V_exothermicity
real*8 omg1,V_reorg1,g_coup1
real*8 omg2,V_reorg2,g_coup2
real*8 gamma_B,VER_rate
real*8 temperature,beta
real*8 omg_c,omg_scaled
real*8 s01,s02,x_cr,V_barrier
real*8,allocatable :: mass(:),omg(:),ck(:)

!! Input/Output
real*8 cnt_frust,cnt_collapse,cnt_init,cnt_term
real*8,allocatable :: pop(:,:),pop_surf(:,:),pop_amp(:,:)
complex*16,allocatable :: rho(:,:,:)

!! Classical
integer nclass,idistribution
real*8,allocatable :: x(:),v(:),acc(:)
real*8,allocatable :: x_old(:),v_old(:),acc_old(:),x_hop(:)
real*8 tim_hop
integer iforward,flag_terminate,flag_frust,flag_reactant,flag_hop,flag_ortho
complex*16,allocatable :: delr(:,:,:),delp(:,:,:),delacc(:,:,:)
complex*16,allocatable :: delr_old(:,:,:),delp_old(:,:,:)

!! Quantized vibration
integer ncl_site
integer nb_vib,n_dvr
real*8,allocatable ::si_sho(:,:,:),sho_overlap(:,:,:,:),q_exp(:,:,:)
real*8,allocatable ::qsq_exp(:,:,:)
real*8,allocatable ::en_sho(:,:),fc_init(:)
real*8,allocatable ::Hamil_diab_0(:,:)

!! Quantum
integer state,nbasis,state_tentative
integer state_old
real*8,allocatable :: si_adiab(:,:),V_k(:),d_ij(:,:,:),vdotd(:,:),V_k_old(:)
real*8,allocatable :: Hamil_site(:,:),Hamil_diab(:,:),delH_dels(:,:,:),delH_dels_ad(:,:,:)
real*8,allocatable :: pot(:,:),force(:,:,:),force_old(:,:,:),delF(:,:,:)
complex*16,allocatable :: ci(:),ci_old(:),sigma(:,:)
real*8,allocatable :: si_adiab_prev(:,:)
complex*16,allocatable :: mat(:,:),mat_adiab(:,:)
real*8,allocatable :: hop_prob(:),W_overlap(:,:),hop_prob_net(:)

!! Evolution
integer n_traj,nsteps,nsteps_eq,nstep_write,iwrite,nstep_avg
real*8 dtc,total_time,curr_time,traj_num,tim_eq
real*8 energy,pot_en,KE_en,energy_cutoff,energy_old
real*8 ensq_avg,en_avg
integer nst_av
integer ihop,icollapse,iterminate,ifriction,iaverage
real*8 prob_tun

!! Parallelization
integer iparallel,iflow,iproc,iwait,ifolder

!! Misc
integer nold,cnt_rate
real*8 tim_tot,tim_ev_cl,tim_diag,tim_cl,tim_rattle,tim_pbc,tim_LJ_tot,tim_solv_solv,tim_check,tim_check2
real*8 tim_T_jk
integer,allocatable:: seed(:)
real*8,allocatable:: work(:)
complex*16,allocatable:: cwork(:)
integer,allocatable:: iwork(:),isuppz(:)

contains
!---------------------------------------------------------- 
!---------------------------------------------------------- 

subroutine setup
  implicit none
  character st_ch
  integer i,size_seed,seed2(2)
  real*8 rnd,c_0,c_e,kt,w

  pi=dacos(-1.d0)
  wave_to_J=2*pi*clight*hbar
  nold=0

  open(10,file="AFSSH.inp")
  read(10,*) iflow
  read(10,*) iproc
  read(10,*) iparallel
  read(10,*) iwait
  read(10,*) N_traj
  read(10,*) dtc
  read(10,*) total_time
  read(10,*) iwrite
  read(10,*) nstep_write
  read(10,*) nstep_avg
  read(10,*) idistribution
  read(10,*) flag_frust
  read(10,*) flag_ortho
  read(10,*) energy_cutoff
  read(10,*) nclass
  read(10,*) nquant
  read(10,*) n_dvr
  read(10,*) nb_vib
  read(10,*) V_coup
  read(10,*) V_exothermicity
  read(10,*) omg1
  read(10,*) V_reorg1
  read(10,*) omg2
  read(10,*) gamma_B
  read(10,*) VER_rate
  read(10,*) temperature
  read(10,*) iforward
  read(10,*) icollapse
  read(10,*) ifriction
  read(10,*) seed2
  read(10,*) st_ch
  close(10)
  !----------------------------------------------------------

  if(st_ch.ne.'x') then
    write(6,*) "problem in reading input file"
    stop
  endif
  !---------------------------------------------------------- 

  nbasis=2*nb_vib
  !nquant=nquant*nb_vib
  !nbasis=nquant
  ncl_site=nclass

  energy_cutoff=energy_cutoff*wave_to_J
  !temperature=temperature*wave_to_J/kb
  omg1=omg1*(2*pi*clight)
  omg2=omg2*(2*pi*clight)
  gamma_B=gamma_B*(2*pi*clight)
  V_exothermicity=V_exothermicity*wave_to_J
  V_coup=V_coup*wave_to_J
  V_reorg1=V_reorg1*wave_to_J
  beta=1.d0/(kb*temperature)

!write(6,*) omg_B/(2*pi*clight),gamma_B/(2*pi*clight)
!stop
  nsteps=nint(total_time/dtc)+1
  beta=1.d0/(kb*temperature)
  !lambda_D=V_reorg/4.d0

  !-----------------------------------------------------------------  
  i=nsteps/nstep_avg+1
  allocate(pop(nquant,i),pop_surf(nquant,i),pop_amp(nquant,i))
  allocate(rho(nbasis,nbasis,i))
  allocate(x(nclass),v(nclass),acc(nclass))
  allocate(x_old(nclass),v_old(nclass),acc_old(nclass),x_hop(nclass))
  allocate(mass(nclass),omg(nclass),ck(nclass))
  allocate(delr(nquant,nquant,nclass),delp(nquant,nquant,nclass),delacc(nquant,nquant,nclass))
  allocate(delr_old(nquant,nquant,nclass),delp_old(nquant,nquant,nclass))
  allocate(si_adiab(nbasis,nquant),ci(nquant),V_k(nquant),V_k_old(nquant),sigma(nquant,nquant))
  allocate(Hamil_site(2,2),Hamil_diab(nbasis,nbasis),delH_dels(nbasis,nbasis,nclass),delH_dels_ad(nquant,nquant,nclass))
  allocate(pot(nquant,nquant),force(nquant,nquant,nclass),force_old(nquant,nquant,nclass),delf(nquant,nquant,nclass))
  allocate(mat(nbasis,nbasis),mat_adiab(nquant,nquant))
  allocate(d_ij(nquant,nquant,nclass),vdotd(nquant,nquant),hop_prob(nquant),W_overlap(nquant,nquant))
  allocate(hop_prob_net(nquant))
  allocate(ci_old(nquant),si_adiab_prev(nbasis,nquant))
  allocate(si_sho(n_dvr,nb_vib,2))
  allocate(sho_overlap(nb_vib,nb_vib,2,2))
  allocate(q_exp(nb_vib,nb_vib,2),qsq_exp(nb_vib,nb_vib,2))
  allocate(en_sho(nb_vib,2),fc_init(nb_vib))
  allocate(Hamil_diab_0(nbasis,nbasis))

  call random_seed(size=size_seed)
  allocate(seed(size_seed))
  do i=1,size_seed/2
    seed(i)=seed2(1)*(2*i+i*i-7)
  enddo
  do i=size_seed/2+1,size_seed
    seed(i)=seed2(2)*(i/2+34-i**3)
  enddo
  call random_seed(put=seed)
  call system_clock(count_rate=cnt_rate)
  !-----------------------------------------------------------------  

  if(iflow==2) then
    open(10,file="ifolder.inp")
    read(10,*) ifolder
    close(10)
    N_traj=N_traj/iparallel
    if(ifolder>1) then
      do i=1,(ifolder-1)*N_traj
        seed=seed+1
!        call init_cond
      enddo
      call random_seed(put=seed)
    endif
  else
    ifolder=1
  endif

  tim_tot=0.d0

end subroutine setup
!---------------------------------------------------------- 

subroutine main
  implicit none
  integer i,j,k,n
  real*8 t1,t2
  real*8 Vc_save

  call cpu_time(t1)

  call setup_parameters
  call rate_estimates

  call cpu_time(t2);tim_tot=tim_tot+t2-t1

end subroutine main
!---------------------------------------------------------- 

subroutine setup_parameters
  implicit none
  integer i
  real*8 c_0,c_e
  real*8 omg_max,delw,w

  mass=1836.d0*au2kg
!  omg_max=3*omg_B
!  omg_c=2*omg_B
!  delw=omg_max/real(nclass)
  w=omg1

  g_coup1=dsqrt(V_reorg1*mass(1)*omg1**2/2.d0)
  !V_reorg2=VER_rate*(1-dexp(-beta*hbar*omg1)) * mass(1)*omg1
  V_reorg2=VER_rate*mass(1)*omg1
  V_reorg2=V_reorg2 * 2 * ((w*w-omg2*omg2)**2+(gamma_B*w)**2)/(omg2**2*gamma_B*w)

  g_coup2=dsqrt(V_reorg2*mass(1)*omg2**2/2.d0)


  Hamil_site=0.d0
  Hamil_site(2,2)=-V_exothermicity

  Hamil_site(1,2)=V_coup; Hamil_site(2,1)=V_coup

!  Hamil_site=Hamil_site*wave_to_J

  omg=omg2
  ck=g_coup2

 omg_scaled=omg1!dsqrt(omg_b**2+sum(ck(1:ncl_site)**2/(mass(1)**2*omg(1:ncl_site)**2)))

end subroutine setup_parameters
!-----------------------------------------------------------------  

subroutine rate_estimates
  real*8 kk(10,10),rate,rate_tot
  real*8 eps,Vc,kt,mm
  real*8 ci_xi,overlap
  real*8 prob
  integer i,n,j

  open(20,file="rates_kinetic_model.out")
  write(6,*) "VER_rate (ps-1)","      rate from kinetic model (ps-1)"
  write(20,*) "VER_rate (ps-1)","     rate from kinetic model (ps-1)"
  do j=1,20
    VER_rate=1.d12+39.d12*(j-1)/19.d0
    call setup_parameters

    mm=mass(1)
    kt=kb*temperature
    q0=g_coup1/(mm*omg1**2)
    omg1_prime=dsqrt(omg1**2+V_reorg2/4.d0*2.d0/mm)
    lambda=V_reorg2*q0**2*(1-V_reorg2/(2*mm*omg1_prime**2))

    rate_tot=0.d0
    do n=1,5
      prob=(1-exp(-hbar*omg1_prime/kt))*exp(-hbar*omg1_prime*(n-1)/kt)
      rate=0.d0
      do i=1,5
        eps=Hamil_site(2,2)-Hamil_site(1,1) +(i-n)*hbar*omg1_prime
        ci_xi=eps/(2*q0*(1-V_reorg2/(2*mm*omg1_prime**2)))
        call calculate_overlap(n,i,ci_xi,overlap)
        Vc=Hamil_site(1,2)*overlap

        kk(n,i)=2*pi/hbar * Vc**2/dsqrt(4*pi*lambda*kt)*dexp(-(lambda+eps)**2/(4*lambda*kt))
        rate=rate+kk(n,i)
      enddo
      rate_tot=rate_tot+rate*prob
    enddo
    write(6,*) VER_rate/1.d12,rate_tot/1.d12
    write(20,*) VER_rate/1.d12,rate_tot/1.d12

  enddo
  close(20)

!do n=1,5
!  do i=1,5
!     write(6,'(es12.3$)') 1.d12/kk(n,i)
!  enddo
!  write(6,*) 
!enddo

end subroutine rate_estimates
!-----------------------------------------------------------------  

subroutine calculate_overlap(state1,state2,ci_xi,overlap)
  implicit none
  integer,intent(in)::state1,state2
  real*8,intent(in)::ci_xi
  real*8,intent(out)::overlap
  integer i,j,k1,k2
  real*8 H_dvr(n_dvr,n_dvr),ke_dvr(n_dvr,n_dvr),x_dvr(n_dvr),q,delq
  real*8 ens(n_dvr),vect(n_dvr,n_dvr)
  real*8 potential
  real*8 tmp(4,4),en_tmp(4),vec_tmp(4,4)

  do i=1,n_dvr
    x_dvr(i)=-5.d-10+10.d-10*(i-1)/real(n_dvr-1)
  enddo
  delq=x_dvr(2)-x_dvr(1)

  call compute_KE_matrix_dvr(KE_dvr,n_dvr,delq,mass(1))
  H_dvr=KE_dvr
  do i=1,n_dvr
    q=x_dvr(i)
    potential = 0.5*mass(1)*omg1_prime**2*(q+q0-1.d0/(mass(1)*omg1_prime**2)*(ci_xi+V_reorg2*q0/2.d0))**2
    H_dvr(i,i)=H_dvr(i,i)+potential
  enddo
  call diag(H_dvr,n_dvr,ens,vect,n_dvr)
  si_sho(:,:,1)=vect(:,1:nb_vib)
  en_sho(:,1)=ens(1:nb_vib)


  H_dvr=KE_dvr
  do i=1,n_dvr
    q=x_dvr(i)
    potential = 0.5*mass(1)*omg1_prime**2*(q-q0-1.d0/(mass(1)*omg1_prime**2)*(ci_xi-V_reorg2*q0/2.d0))**2
    H_dvr(i,i)=H_dvr(i,i)+potential
  enddo
  call diag(H_dvr,n_dvr,ens,vect,n_dvr)
  si_sho(:,:,2)=vect(:,1:nb_vib)
  en_sho(:,2)=ens(1:nb_vib)

  overlap=sum(si_sho(:,state1,1)*si_sho(:,state2,2))

end subroutine calculate_overlap
!-----------------------------------------------------------------  

subroutine compute_KE_matrix_dvr(KE_matrix,ngrid,delq,mass)
  !! computes KE matrix in DVR basis
  !! Appendix A of JCP 96, 1982 (1991)
  implicit none
  integer,intent(in) :: ngrid       !! size of DVR grid
  real*8,intent(inout) :: KE_matrix(ngrid,ngrid)
  real*8,intent(in) :: delq         !! step size in position of DVR basis
  real*8,intent(in) :: mass         !! mass
  integer i,j
  real*8 pi,hbar

  pi=dacos(-1.d0);hbar=1.05457266D-34

  KE_matrix=hbar**2/(2*mass*delq**2)
  do i=1,ngrid
    do j=1,ngrid
      KE_matrix(i,j)=KE_matrix(i,j)*(-1.d0)**(i-j)
      if(i==j) then
        KE_matrix(i,j)=KE_matrix(i,j)*pi**2/3.d0
      else
        KE_matrix(i,j)=KE_matrix(i,j)*2.d0/real(i-j)**2
      endif
    enddo
  enddo
end subroutine compute_KE_matrix_dvr
!---------------------------------------------------------- 
subroutine diag(mat,n,eigen_value,eigen_vect,m_values)
  !! Diaganalizing matrix using dsyevr. First m_values eigen values and eigenvectors computed.
  !! The module's common variables should contain:

  !! Initialize nold=0 

  !! nold makes sure that everytime value of n changes, work and iwork are re-allocated for optimal performance.
  !! mat is destroyed after use.

  implicit none
  integer,intent(in) :: n,m_values
  real*8,intent(out) :: eigen_value(n),eigen_vect(n,m_values)
  real*8,intent(inout) :: mat(n,n)
  real*8 vl,vu,abstol
  integer il,iu,info,m,AllocateStatus
  integer lwork,liwork

  vl=0.d0;vu=0.d0   !! not referenced
  il=1;iu=m_values
  abstol=0.d0
  info=0

  if(nold.ne.n .or. .not.allocated(work) .or. .not.allocated(iwork) .or. .not.allocated(isuppz)) then
  !if(nold.ne.n) then
    lwork=-1;liwork=-1
    if(allocated(isuppz))deallocate(isuppz)
    if(allocated(work))deallocate(work)
    if(allocated(iwork))deallocate(iwork)
    allocate(isuppz(2*m_values),work(n),iwork(n))
    call dsyevr('V','I','U',n,mat,n,vl,vu,il,iu,abstol,m,eigen_value,eigen_vect,n,isuppz,work,lwork,iwork,liwork,info)
    lwork=nint(work(1)); liwork=iwork(1)
    deallocate(work,iwork)
    allocate(work(lwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    allocate(iwork(liwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    nold=n
  endif

  lwork=size(work)
  liwork=size(iwork)

  call dsyevr('V','I','U',n,mat,n,vl,vu,il,iu,abstol,m,eigen_value,eigen_vect,n,isuppz,work,lwork,iwork,liwork,info)
  if(info.ne.0) then
    write(6,*) "problem in diagonalization",info
    stop
  endif

end subroutine diag
!---------------------------------------------------------- 



End Module mod_afssh
