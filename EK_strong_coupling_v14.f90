!##############################################################################
!######              Eguchi-Kawai model with continuous time          #########
!######                                                               #########
!######                 written by Masanori Hanada                    #########
!######                                                               #########
!######            ver.1 modified from "tetrahedron_v3.f90"           #########
!######                  strong coupling limit                        #########
!######            ver.2 measurement of Polyakov loop correlator added. #######
!######            ver.3 normalization of <Energy> corrected.         #########
!######           ver.11 gauge fix to confined and deconfined sectors #########
!######           ver.12 a little bit of cleaning up                  #########  
!######           ver.13 Mersenne twister                             #########
!######           ver.14 OpenMP                                       #########
!##############################################################################
!Mersenne twister.
include 'mt19937.f90'
program Eguchi_Kawai

  use mtmod !Mersenne twistor
  implicit none

  include 'size_EK_v14.inc'
  !---------------------------------
  character(150) input_config,data_output,output_config,output_phase,&
       output_lambda_dec,output_cor_con,output_cor_dec,output_cor_mix
  !----------------------------------------
  double precision LatticeSpacing !lattice spacing a
  double precision temperature,mass
  double complex condensate_sum
  !----------------------------------------
  !     Iteration
  integer iteration
  integer NoSkip
  !----------------------------------------
  !for correlator of Pol's
  integer nprod
  parameter(nprod=15)
  !    initial configuration
  !    IniConfig=0 -> new start
  !    IniConfig=1 -> use old configuration
  integer IniConfig
  !----------------------------------------
  !   ndirac = 1 -> calculate the Dirac spectrum
  integer ndirac
  !----------------------------------------
  !   nseed=0 -> use mersenne_seed in the input file (always so if iniconfig=0)
  integer nseed
  !----------------------------------------
  !     Parameters for Molecular Evolution 
  !     Ntau: # of steps, Dtau: width of each step    
  integer Ntau
  doubleprecision Dtau_U,Dtau_alpha
  !----------------------------------------
  doublecomplex U(1:NSite,1:NDIM,1:NMAT,1:NMAT)
  double precision alpha(1:NMAT)
  doublecomplex Pol,MAT_test(1:NMAT,1:NMAT),pol_dec,pol_con
  double precision phase(1:NMAT)
  integer seedsize,mersenne_seed
  integer,allocatable:: seed(:)
  double precision P_fix,P_fix_width,g_coeff,wilson,&
       av_plaq,ransu,g_coeff_alpha
  !min_diag,deviation_from_one,deviation_max

  !---------------------------------------
  integer i,j,k,l,ia,ja,a,p,q,t,time,i_measure_dirac,iperm,perm_count,ncorrelator
  doubleprecision hantei,HamStart,HamEnd,tot_ene,dH

  integer acceptance,trial,NumCalObservable,rejection,CheckHam,ite,info,nphase
  double complex eigenvalues(1:4*nmat*nsite),cond_dec,cond_con

  double precision cor_dec(1:nprod),cor_con(1:nprod),cor_mix(1:nprod)
  
  !call random_seed(size=seedsize)
  !allocate(seed(seedsize))
  

  open(unit=10,status='OLD',file='input_v14.dat',action='READ')
  read(10,*)input_config
  read(10,*)data_output
  read(10,*)output_phase
  read(10,*)output_cor_dec
  read(10,*)output_cor_con
  read(10,*)output_cor_mix
  read(10,*)output_config
  read(10,*)IniConfig
  read(10,*)nseed
  read(10,*)temperature
  read(10,*)ncorrelator
  read(10,*)nphase
  read(10,*)iteration
  read(10,*)NoSkip
  read(10,*)Ntau
  read(10,*)Dtau_U
  read(10,*)Dtau_alpha
  read(10,*)g_coeff
  read(10,*)P_fix
  read(10,*)P_fix_width
  read(10,*)g_coeff_alpha
  read(10,*)mersenne_seed
  close(10)
  
  LatticeSpacing=1d0/temperature/dble(NSite)
  perm_count=0
  !*************************************
  !*** Set the initial configuration ***
  !*************************************
  !new start
  if(IniConfig.EQ.0)then
     !call random_seed()
     call sgrnd(mersenne_seed)!use the seed in the input file
     U=(0d0,0d0)
     do t=1,NSite
        do p=1,NDIM
           do i=1,NMAT
              U(t,p,i,i)=(1d0,0d0)
           end do
        end do
     end do
     do i=1,NMAT
        !call random_number(ransu)
        ransu=grnd()
        alpha(i)=ransu
     end do
     ! Initial configuration specified   
  else if(IniConfig.EQ.1)then
     open(unit=9,status='OLD',file=input_config,action='READ')
     !read(9,*) seed
     read(9,*) U
     read(9,*) alpha
     IF(nseed.eq.1)then
        call mtgetu(9)!use the seed in the config file
     else
        call sgrnd(mersenne_seed)!use the seed in the input file
     end IF
     close(9)
     !call random_seed(put=seed)       
  end if
  !write(*,*)"seed=", seed
  !**************************************************
  !**************************************************
  acceptance=0
  trial=0

  NumCalObservable=0

  !************************************
  !************************************
  !     Make the output file
  !************************************
  !************************************  

  open(unit=10,status='REPLACE',file=data_output,action='WRITE')
  write(10,*) "#size of the gauge group: NMAT=",NMAT
  write(10,*) "#size of the deconfined sector: NDEC=",NDEC
  write(10,*) "#number of sites=",NSite
  write(10,*) "#Lattice Spacing=",LatticeSpacing
  write(10,*) "#temperature=",temperature
  !write(10,*) "#probe fermion mass=",mass
  write(10,*) "#'t Hooft coupling = infinity"
  write(10,*) "#Ntau=",Ntau
  write(10,*) "#Dtau for U=",Dtau_U
  write(10,*) "#Dtau for alpha=",Dtau_alpha
  write(10,*) "#P_dec is fixed to:",P_fix, "+/-",P_fix_width
  write(10,*) "#P_con is fixed below:",P_fix_width
  write(10,*) "#Coefficient for fixing P_dec and P_con =",g_coeff 
  write(10,*) "# sweep, perm_count, dH, |Pol_dec.|, |Pol_con.|, action, |Wilson|, acceptance rate"
  write(10,*)'#------------------------------------------------'

  open(unit=11,status='REPLACE',file=output_phase,action='WRITE')
  if(ncorrelator.eq.1)then
     open(unit=12,status='REPLACE',file=output_cor_dec,action='WRITE')
     open(unit=13,status='REPLACE',file=output_cor_con,action='WRITE')
     open(unit=14,status='REPLACE',file=output_cor_mix,action='WRITE')
  end if
  !************************************
  !************************************
  !     Take expectation value 
  !************************************
  !************************************       

  do l=1,iteration
  
     call MolecularEvolution(NMAT,NDIM,NSite,Ntau,Dtau_U,Dtau_alpha,U,alpha,&
          LatticeSpacing,acceptance,trial,g_coeff,P_fix,P_fix_width,dH,ndec,g_coeff_alpha)
     if(ndec.LT.nmat)then
        do iperm=1,10
           call permutation(nmat,ndim,nsite,U,alpha,ndec,g_coeff,P_fix,P_fix_width,perm_count)
        end do
     end if
     ! Calculate physical quantities 
     if(MOD(l,NoSkip).EQ.0)then
        NumCalObservable=NumCalObservable+1         
        call CalcPol_no_phase_2(NMAT,NSite,alpha,Pol_dec,Pol_con,ndec)
        call total_energy(tot_ene,U,alpha,NMAT,NDIM,NSite,LatticeSpacing)
        tot_ene=tot_ene/dble(nmat*nmat)
        call Calc_Wilson(U,NMAT,NDIM,NSite,wilson)

        call CalcPol(NMAT,alpha,Pol,phase)
        if(nphase.EQ.1)then
           write(11,*)phase
        end if
        if(ncorrelator.EQ.1)then
           call Calc_Pol_Correlator(NMAT,nsite,ndim,alpha,cor_dec,cor_con,cor_mix,U,ndec,nprod)
        
           write(12,*)cor_dec
           write(13,*)cor_con
           write(14,*)cor_mix
        end if
        write(10,*)l,perm_count,dH,abs(Pol_dec),abs(Pol_con),tot_ene,wilson,&
        &dble(acceptance)/dble(trial)
        write(*,*)l,perm_count,dH,abs(Pol_dec),abs(Pol_con),tot_ene,wilson,&
        &dble(acceptance)/dble(trial)

     end if

  end do
  !**************************************************
  !**************************************************
  !   End of iteration
  !**************************************************
  !**************************************************
  close(10)
  close(11)
  close(12)
  close(13)
  close(14)
  
  
  !call random_seed(get=seed)


  open(UNIT = 22, File = output_config, STATUS = "REPLACE", ACTION = "WRITE")
  !write(22,*) seed
  write(22,*) U
  write(22,*) alpha
  call mtsaveu(22)
  close(22)

  !write(*,*)"seed=", seed




end program Eguchi_Kawai

!***********************************************************
!***********************************************************
!     Box-Muller method for generating Gaussian random number

SUBROUTINE BoxMuller(p,q)  

  use mtmod !Mersenne twistor
  implicit none 

  doubleprecision p,q,r,s,Pi

  Pi=2d0*DASIN(1d0)
  !call random_number(r)
  !call random_number(s)
  r=grnd()
  s=grnd()
  p=dsqrt(-2d0*dlog(r))*DSIN(2d0*Pi*s)
  q=dsqrt(-2d0*dlog(r))*DCOS(2d0*Pi*s)

  return

END SUBROUTINE BoxMuller


!***********************************************************
!***********************************************************
!     the kinetic term; (-1/2a)*Tr(Vp*Upq*Vq^dag*Upq) + c.c. 
SUBROUTINE kinetic(kin_total,U,alpha,NMAT,NDIM,NSite,LatticeSpacing) 

  use omp_lib
  implicit none

  integer t,tp1(1:nsite),imat,jmat,idim,NMAT,NSite,NDIM
  doublecomplex U(1:NSite,1:NDIM,1:NMAT,1:NMAT)
  double precision alpha(1:nmat)
  double complex expialpha(1:nmat)
  
  doubleprecision kin(1:nsite),kin_total,LatticeSpacing
  

  do imat=1,nmat
     expialpha(imat)=dcmplx(dcos(alpha(imat)/dble(nsite)))+(0d0,1d0)*dcmplx(dsin(alpha(imat)/dble(nsite)))
  end do

  kin_total=0d0
  !$OMP PARALLEL DO REDUCTION(+:kin_total)
  do t=1,nsite
     !write(*,*) 'Thread ', omp_get_thread_num()
     if(t.NE.NSite)then
        tp1(t)=t+1
     else
        tp1(t)=1
     end if

     do idim=1,ndim
        kin(t)=0d0
        do imat=1,NMAT
           do jmat=1,NMAT
              kin(t)=kin(t)+dble(expialpha(imat)*U(tp1(t),idim,imat,jmat)&
                   *dconjg(expialpha(jmat))*dconjg(U(t,idim,imat,jmat)))
           end do
        end do
        kin(t)=dble(NMAT)-kin(t)
        kin(t)=kin(t)/LatticeSpacing*dble(NMAT)
        kin_total=kin_total+kin(t)
     end do
  end do
  !$OMP END PARALLEL DO
     
  return

END SUBROUTINE kinetic
!***********************************************************
!***********************************************************
!     the kinetic term; (-1/2a)*Tr(Vp*Upq*Vq^dag*Upq) + c.c. 

SUBROUTINE Calc_Wilson(U,NMAT,NDIM,NSite,wilson) 

  implicit none

  integer NMAT,NSite,NDIM
  integer t, imat,idim
  doublecomplex U(1:NSite,1:NDIM,1:NMAT,1:NMAT),temp(1:nsite)
  
  doubleprecision wilson

  wilson=0d0
  !$OMP PARALLEL DO REDUCTION(+:wilson)
  do t=1,nsite
     do idim=1,ndim
        temp(t)=(0d0,0d0)
        do imat=1,nmat
           temp(t)=temp(t)+U(t,idim,imat,imat)
        end do
        wilson=wilson+abs(temp(t))
     end do
  end do
  !$OMP END PARALLEL DO
  wilson=wilson/dble(nsite*ndim*nmat)

  return

END SUBROUTINE Calc_Wilson
!***********************************************************
!***********************************************************
!   Derivative of Hamiltonian w.r.t. V(t,p,i,j)
SUBROUTINE MakeDerivHamiltonian(U,alpha,NMAT,NDIM,NSite,DerivH_U,DerivH_alpha,&
     LatticeSpacing,g_coeff,P_fix,P_fix_width,ndec,g_coeff_alpha)


  implicit none 
  

  integer t,NSite,NMAT,NDIM,imat,jmat,idim,ndec,isite

  doublecomplex U(1:NSite,1:NDIM,1:NMAT,1:NMAT)
  doublecomplex V(1:NMAT)
  double precision alpha(1:nmat)
  doubleprecision LatticeSpacing,g_coeff,P_fix,P_fix_width,g_coeff_alpha
  doublecomplex DerivH_U(1:NSite,1:NDIM,1:NMAT,1:NMAT)
  doubleprecision DerivH_alpha(1:NMAT)

  do imat=1,nmat
     V(imat)=dcmplx(dcos(alpha(imat)/dble(nsite)))+(0d0,1d0)*dcmplx(dsin(alpha(imat)/dble(nsite)))
  end do
     

  call MakeDerivHamiltonian_alpha(U,alpha,NMAT,NDIM,NSite,DerivH_alpha,LatticeSpacing,&
       g_coeff,P_fix,P_fix_width,ndec,g_coeff_alpha)


  call MakeDerivHamiltonian_U(U,V,NMAT,NDIM,NSite,DerivH_U,LatticeSpacing)
  
  
  return
  
END SUBROUTINE MakeDerivHamiltonian


!***********************************************************
!***********************************************************
!     Molecular evolution (leap frog method)
SUBROUTINE MolecularEvolution(NMAT,NDIM,NSite,Ntau,Dtau_U,Dtau_alpha,U,alpha,&
     LatticeSpacing,acceptance,trial,g_coeff,P_fix,P_fix_width,dH,ndec,g_coeff_alpha)

  use mtmod !Mersenne twistor
  implicit none
  
  integer t,l,p,q,r,NSite,Ntau,NMAT,acceptance,trial,ite,NDIM,idim,ndec
  integer imat,jmat,kmat,isite
  doublecomplex BackUp_U(1:NSite,1:NDIM,1:NMAT,1:NMAT)
  double precision BackUp_alpha(1:NMAT)
  doublecomplex U(1:NSite,1:NDIM,1:NMAT,1:NMAT)
  double precision alpha(1:nmat),P_alpha(1:NMAT)
  double precision alpha_max,alpha_min,pi,g_coeff_alpha
  doubleprecision Dtau_U,Dtau_alpha,ransu1,ransu2,HamStart,HamEnd,&
       LatticeSpacing,kin,hantei,g_coeff,P_fix,P_fix_width,tot_kin,dH
  doublecomplex DerivH_U(1:NSite,1:NDIM,1:NMAT,1:NMAT)
  doubleprecision DerivH_alpha(1:NMAT)
  double complex P_U(1:NSite,1:NDIM,1:NMAT,1:NMAT),&
       MAT(1:NMAT,1:NMAT),MATiEXP(1:NMAT,1:NMAT),&
       U_Old(1:NMAT,1:NMAT),pol_dec,pol_con
  double complex MATiEXP_2(1:nmat,1:nmat,1:nsite),MAT_2(1:nmat,1:nmat,1:nsite)
  doublecomplex U_Old_2(1:NSite,1:NMAT,1:NMAT)
  
  BackUp_U=U
  BackUp_alpha=alpha

  pi=2d0*dasin(1d0)
  
  ! initial condition for P_U, P_alpha
  do imat=1,NMAT
     call BoxMuller(ransu1,ransu2)
     P_alpha(imat)=dcmplx(ransu1)
  end do
  do t=1,nsite
     do idim=1,ndim
        do imat=1,NMAT-1
           do jmat=imat+1,NMAT
              call BoxMuller(ransu1,ransu2)
              P_U(t,idim,imat,jmat)=dcmplx(ransu1/dsqrt(2d0))+dcmplx(ransu2/dsqrt(2d0))*(0D0,1D0)
              P_U(t,idim,jmat,imat)=dcmplx(ransu1/dsqrt(2d0))-dcmplx(ransu2/dsqrt(2d0))*(0D0,1D0)
           end do
        end do
        do imat=1,NMAT
           call BoxMuller(ransu1,ransu2)
           P_U(t,idim,imat,imat)=dcmplx(ransu1)
        end do
     end do
  end do
  ! Calculate the Hamiltonian
  !call total_kinetic_energy(tot_kin,U,V,NMAT,DIM,NumSite,LatticeSpacing)
  call total_energy(HamStart,U,alpha,NMAT,NDIM,NSite,LatticeSpacing)
  do imat=1,NMAT
     HamStart=HamStart+P_alpha(imat)*P_alpha(imat)*0.5d0
  end do

  do t=1,nsite
     do idim=1,ndim
        do imat=1,NMAT
           do jmat=1,NMAT
              HamStart=HamStart+dble(P_U(t,idim,imat,jmat)*dconjg(P_U(t,idim,imat,jmat)))*0.5d0
           end do
        end do
     end do
  end do
  !****************************************
  !*** Fix Pol; dec and con, separately ***
  !****************************************
  call CalcPol_no_phase_2(NMAT,NSite,alpha,Pol_dec,Pol_con,ndec)
  if(abs(pol_dec).GT.P_fix+P_fix_width)then
     HamStart=HamStart+0.5d0*g_coeff*(abs(pol_dec)-P_fix-P_fix_width)**2d0
  else if(abs(pol_dec).LT.P_fix-P_fix_width)then
     HamStart=HamStart+0.5d0*g_coeff*(abs(pol_dec)-P_fix+P_fix_width)**2d0
  end if
  !pol_con=0d0 if ndec=nmat.
  if(abs(pol_con).GT.P_fix_width)then
     HamStart=HamStart+0.5d0*g_coeff*(abs(pol_con)-P_fix_width)**2d0
  end if
  !***************
  !*** FP term ***
  !***************
  do imat=1,nmat-1
     do jmat=imat+1,nmat
        HamStart=HamStart-2d0*dlog(dabs(dsin(0.5d0*(alpha(imat)-alpha(jmat)))))
     end do
  end do
  !****************************  
  !*** constraint for alpha ***
  !****************************
  alpha_max=alpha(1)
  alpha_min=alpha(1)
  do imat=2,nmat
     if(alpha(imat).GT.alpha_max)then
        alpha_max=alpha(imat)
     else if(alpha(imat).LT.alpha_min)then
        alpha_min=alpha(imat)
     end if
  end do
  if(alpha_max-alpha_min.LT.2d0*pi)then
     hamstart=hamstart-dlog(2d0*pi-(alpha_max-alpha_min))
  end if
  !*********** leap frog : start *********************
  ! first step of leap frog
  !$OMP PARALLEL DO
  do t=1,NSite
     do idim=1,ndim
        do imat=1,nmat
           do jmat=1,nmat
              MAT_2(imat,jmat,t)=P_U(t,idim,imat,jmat)*dcmplx(0.5d0*Dtau_U)
           end do
        end do
        call MATRIX_iEXP_2(NMAT,MAT_2,MATiEXP_2,t,nsite)
        do imat=1,NMAT
           do jmat=1,NMAT
              U_Old_2(t,imat,jmat)=U(t,idim,imat,jmat)
              U(t,idim,imat,jmat)=(0d0,0d0)
           end do
        end do
        do imat=1,NMAT
           do jmat=1,NMAT
              do kmat=1,NMAT
                 U(t,idim,imat,jmat)=U(t,idim,imat,jmat)+MATiEXP_2(imat,kmat,t)*U_Old_2(t,kmat,jmat) 
              end do
           end do
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  alpha=alpha+P_alpha*dtau_alpha*0.5d0
  ! second,...,Ntau-th   
  do ite=2,Ntau
     call MakeDerivHamiltonian(U,alpha,NMAT,NDIM,NSite,DerivH_U,&
          DerivH_alpha,LatticeSpacing,g_coeff,P_fix,P_fix_width,ndec,g_coeff_alpha)
     P_U=P_U-DerivH_U*dcmplx(Dtau_U)
     P_alpha=P_alpha-DerivH_alpha*Dtau_alpha
     !$OMP PARALLEL DO 
     do t=1,NSite
        do idim=1,ndim
           do imat=1,nmat
              do jmat=1,nmat
                 MAT_2(imat,jmat,t)=P_U(t,idim,imat,jmat)*dcmplx(Dtau_U)
              end do
           end do
           call MATRIX_iEXP_2(NMAT,MAT_2,MATiEXP_2,t,nsite)
           do imat=1,NMAT
              do jmat=1,NMAT
                 U_Old_2(t,imat,jmat)=U(t,idim,imat,jmat)
                 U(t,idim,imat,jmat)=(0d0,0d0)
              end do
           end do
           do imat=1,NMAT
              do jmat=1,NMAT
                 do kmat=1,NMAT
                    U(t,idim,imat,jmat)=U(t,idim,imat,jmat)+MATiEXP_2(imat,kmat,t)*U_Old_2(t,kmat,jmat) 
                 end do
              end do
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     alpha=alpha+P_alpha*dtau_alpha
  end do
  ! last step
  call MakeDerivHamiltonian(U,alpha,NMAT,NDIM,NSite,DerivH_U,DerivH_alpha,&
       LatticeSpacing,g_coeff,P_fix,P_fix_width,ndec,g_coeff_alpha)
  
  P_U=P_U-DerivH_U*dcmplx(Dtau_U)
  P_alpha=P_alpha-DerivH_alpha*Dtau_alpha

  !$OMP PARALLEL DO 
  do t=1,NSite
     do idim=1,ndim
        do imat=1,nmat
           do jmat=1,nmat
              MAT_2(imat,jmat,t)=P_U(t,idim,imat,jmat)*dcmplx(0.5d0*Dtau_U)
           end do
        end do
        call MATRIX_iEXP_2(NMAT,MAT_2,MATiEXP_2,t,nsite)
        do imat=1,NMAT
           do jmat=1,NMAT
              U_Old_2(t,imat,jmat)=U(t,idim,imat,jmat)
              U(t,idim,imat,jmat)=(0d0,0d0)
           end do
        end do
        do imat=1,NMAT
           do jmat=1,NMAT
              do kmat=1,NMAT
                 U(t,idim,imat,jmat)=U(t,idim,imat,jmat)+MATiEXP_2(imat,kmat,t)*U_Old_2(t,kmat,jmat) 
              end do
           end do
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  alpha=alpha+P_alpha*dtau_alpha*0.5d0
  !*********** leap frog : end *********************
  ! Calculate the Hamiltonian
  call total_energy(HamEnd,U,alpha,NMAT,NDIM,NSite,LatticeSpacing)
  do imat=1,NMAT
     HamEnd=HamEnd+P_alpha(imat)*P_alpha(imat)*0.5d0
  end do
  do t=1,NSite
     do idim=1,ndim
        do imat=1,NMAT
           do jmat=1,NMAT
              HamEnd=HamEnd+dble(P_U(t,idim,imat,jmat)*dconjg(P_U(t,idim,imat,jmat)))*0.5d0
           end do
        end do
     end do
  end do
  !****************************************
  !*** Fix Pol; dec and con, separately ***
  !****************************************
  call CalcPol_no_phase_2(NMAT,NSite,alpha,Pol_dec,Pol_con,ndec)
  if(abs(pol_dec).GT.P_fix+P_fix_width)then
     HamEnd=HamEnd+0.5d0*g_coeff*(abs(pol_dec)-P_fix-P_fix_width)**2d0
  else if(abs(pol_dec).LT.P_fix-P_fix_width)then
     HamEnd=HamEnd+0.5d0*g_coeff*(abs(pol_dec)-P_fix+P_fix_width)**2d0
  end if
  !pol_con=0d0 if ndec=nmat.
  if(abs(pol_con).GT.P_fix_width)then
     HamEnd=HamEnd+0.5d0*g_coeff*(abs(pol_con)-P_fix_width)**2d0
  end if
  !***************
  !*** FP term ***
  !***************
  do imat=1,nmat-1
     do jmat=imat+1,nmat
        HamEnd=HamEnd-2d0*dlog(dabs(dsin(0.5d0*(alpha(imat)-alpha(jmat)))))
     end do
  end do
  !****************************  
  !*** constraint for alpha ***
  !****************************
  alpha_max=alpha(1)
  alpha_min=alpha(1)
  do imat=2,nmat
     if(alpha(imat).GT.alpha_max)then
        alpha_max=alpha(imat)
     else if(alpha(imat).LT.alpha_min)then
        alpha_min=alpha(imat)
     end if
  end do
  if(alpha_max-alpha_min.LT.2d0*pi)then
     hamend=hamend-dlog(2d0*pi-(alpha_max-alpha_min))
  end if

  
  trial=trial+1 
  ! metropolis test

  if(alpha_max-alpha_min.GT.2d0*pi)then
     !automatic reject
     U=BackUp_U
     alpha=BackUp_alpha
  else if(alpha_max-alpha_min.LT.2d0*pi)then
  
     !call random_number(hantei)
     hantei=grnd()
     if(dexp(HamStart-HamEnd) > hantei)THEN
        acceptance=acceptance+1
     else 
        U=BackUp_U
        alpha=BackUp_alpha
     end if
  end if

  dH=-HamStart+HamEnd

  return
  
END SUBROUTINE MolecularEvolution
!***********************************************************
!***********************************************************
!   Derivative of Hamiltonian w.r.t. V(t,p,i,j)

SUBROUTINE MakeDerivHamiltonian_alpha(U,alpha,NMAT,DIM,NumSite,DerivH,LatticeSpacing,g_coeff,P_fix,P_fix_width,ndec,g_coeff_alpha)

  implicit none 
  

  integer t,l,p,NumSite,NMAT,DIM,it,imat,jmat,kmat,ndec,imax,imin,tp1(1:numsite)

  doublecomplex U(1:NumSite,1:DIM,1:NMAT,1:NMAT)
  doubleprecision alpha(1:NMAT),alpha_max,alpha_min,g_coeff_alpha,pi
  
  doubleprecision LatticeSpacing,g_coeff,P_fix,P_fix_width
  doubleprecision DerivH(1:NMAT)
  doublecomplex temp(1:NMAT,1:NMAT)

  double complex MAT1(1:NMAT,1:NMAT),MAT2(1:NMAT,1:NMAT)

  double complex pol_dec,pol_con
  double complex expialpha(1:nmat),polyakov_line(1:nmat)

  DerivH=0d0

  pi=2d0*dasin(1d0)
  
  do imat=1,nmat
     expialpha(imat)=dcmplx(dcos(alpha(imat)/dble(numsite)))+(0d0,1d0)*dcmplx(dsin(alpha(imat)/dble(numsite)))
  end do

  temp=(0d0,0d0)
  do t=1,numsite
     if(t.NE.NumSite)then
        tp1(t)=t+1
     else
        tp1(t)=1
     end if
     
     do p=1,DIM
        do imat=1,NMAT
           do jmat=1,NMAT
              temp(imat,jmat)=temp(imat,jmat)&
                   +U(tp1(t),p,imat,jmat)*dconjg(U(t,p,imat,jmat))&
                   *expialpha(imat)*dconjg(expialpha(jmat))
              
           end do
        end do
     end do
  end do
  temp=temp*(-1d0)*dcmplx(nmat)/dcmplx(2d0*LatticeSpacing)/dble(numsite)
  do imat=1,nmat
     do jmat=1,nmat
        DerivH(imat)=DerivH(imat)+dble((0d0,1d0)*temp(imat,jmat)-(0d0,1d0)*temp(jmat,imat))*2d0
     end do
  end do
  
  !********************************
  !*** deriv of constraint term ***
  !********************************
  do imat=1,NMAT
     polyakov_line(imat)=dcmplx(dcos(alpha(imat)))+(0d0,1d0)*dcmplx(dsin(alpha(imat)))
  end do

  Pol_dec=(0d0,0d0)
  do imat=1,ndec
     Pol_dec=Pol_dec+polyakov_line(imat)
  end do
  Pol_dec=Pol_dec/dcmplx(ndec)

  pol_con=(0d0,0d0)
  if(ndec.LT.nmat)then
     do imat=ndec+1,nmat
        Pol_con=Pol_con+polyakov_line(imat)
     end do
     Pol_con=Pol_con/dcmplx(nmat-ndec)
  end if
  
  if(abs(pol_dec).GT.P_fix+P_fix_width)then
     do imat=1,ndec
        DerivH(imat)=DerivH(imat)+g_coeff*(abs(pol_dec)-P_fix-P_fix_width)/abs(pol_dec)/dble(ndec)&
             *dble(dconjg(pol_dec)*(0d0,1d0)*polyakov_line(imat))
     end do
  else if(abs(pol_dec).LT.P_fix-P_fix_width)then
     do imat=1,ndec
        DerivH(imat)=DerivH(imat)+g_coeff*(abs(pol_dec)-P_fix+P_fix_width)/abs(pol_dec)/dble(ndec)&
             *dble(dconjg(pol_dec)*(0d0,1d0)*polyakov_line(imat))
     end do
  end if
  !pol_con=0 if ndec=nmat
  if(abs(pol_con).GT.P_fix_width)then
     do imat=ndec+1,nmat
        DerivH(imat)=DerivH(imat)+g_coeff*(abs(pol_con)-P_fix_width)/abs(pol_con)/dble(nmat-ndec)&
             *dble(dconjg(pol_con)*(0d0,1d0)*polyakov_line(imat))
     end do
  end if
  

  !*************************
  !*** gauge-fixing term ***
  !*************************
  do imat=1,nmat
     do jmat=1,nmat
        if(jmat.NE.imat)then
           derivh(imat)=derivh(imat)& 
                -1d0/dtan(0.5d0*(alpha(imat)-alpha(jmat)))
        end if
     end do
  end do
  !****************************
  !****************************  
  !*** constraint for alpha ***
  !****************************
  !****************************
  imax=1
  imin=1
  alpha_max=alpha(1)
  alpha_min=alpha(1)
  do imat=2,nmat
     if(alpha(imat).GT.alpha_max)then
        imax=imat
        alpha_max=alpha(imat)
     else if(alpha(imat).LT.alpha_min)then
        imin=imat
        alpha_min=alpha(imat)
     end if
  end do
  if(alpha_max-alpha_min.LT.2d0*pi)then
     derivh(imax)=derivh(imax)&
          +1d0/(2d0*pi-(alpha_max-alpha_min)+1d0/g_coeff_alpha)
     derivh(imin)=derivh(imin)&
          -1d0/(2d0*pi-(alpha_max-alpha_min)+1d0/g_coeff_alpha)
  else if(alpha_max-alpha_min.GE.2d0*pi)then
     derivh(imax)=derivh(imax)+g_coeff_alpha
     derivh(imin)=derivh(imin)-g_coeff_alpha
  end if
  
  
  return
  
END SUBROUTINE MakeDerivHamiltonian_Alpha
!***********************************************************
!***********************************************************
!   Derivative of Hamiltonian w.r.t. U(t,p,q,i,j); p<q
SUBROUTINE MakeDerivHamiltonian_U(U,V,NMAT,NDIM,NSite,DerivH,LatticeSpacing)

  implicit none 
  
  integer t,tp1(1:nsite),tm1(1:nsite),l,q,r,NSite,NMAT,NDIM,idim,jdim,imat,jmat,kmat
  
  doublecomplex U(1:NSite,1:NDIM,1:NMAT,1:NMAT)
  doublecomplex V(1:NMAT)
  doubleprecision LatticeSpacing
  doublecomplex DerivH(1:nsite,1:ndim,1:NMAT,1:NMAT),DerivH_pre(1:nsite,1:ndim,1:NMAT,1:NMAT)
  doublecomplex UVd(1:NSITE,1:NMAT,1:NMAT)
  doublecomplex UV(1:NSITE,1:NMAT,1:NMAT)
  doublecomplex UdVd(1:NSITE,1:NMAT,1:NMAT)
  doublecomplex UdV(1:NSITE,1:NMAT,1:NMAT)
  !doublecomplex UU(1:NSITE,1:NMAT,1:NMAT)
  !doublecomplex UdUd(1:NSITE,1:NMAT,1:NMAT)
  !doublecomplex UUd(1:NSITE,1:NMAT,1:NMAT)
  !doublecomplex UdU(1:NSITE,1:NMAT,1:NMAT)
  !doublecomplex temp(1:NSITE,1:NMAT,1:NMAT)
  
  
  DerivH_pre=(0d0,0d0)
  !V is diagonal
  !*************************
  !***** kinetic term ******
  !*************************
  !$OMP PARALLEL DO 
  do t=1,nsite
     if(t.NE.NSite)then
        tp1(t)=t+1
     else
        tp1(t)=1
     end if
     
     if(t.NE.1)then
        tm1(t)=t-1
     else
        tm1(t)=NSite
     end if
     
     do idim=1,ndim
        do imat=1,nmat
           do jmat=1,nmat
              UVd(t,imat,jmat)=(0d0,0d0)
              UdV(t,imat,jmat)=(0d0,0d0)
           end do
        end do
        do imat=1,NMAT
           do jmat=1,NMAT
              UVd(t,imat,jmat)=UVd(t,imat,jmat)+U(t,idim,imat,jmat)*dconjg(V(jmat))
              UdV(t,imat,jmat)=UdV(t,imat,jmat)+dconjg(U(tm1(t),idim,jmat,imat))*V(jmat)
           end do
        end do
        
        do imat=1,NMAT
           do jmat=1,NMAT
              do kmat=1,NMAT
                 DerivH_pre(t,idim,imat,jmat)=DerivH_pre(t,idim,imat,jmat)+UVd(t,imat,kmat)*UdV(t,kmat,jmat)
              end do
           end do
        end do
        do imat=1,nmat
           do jmat=1,nmat
              UdVd(t,imat,jmat)=(0d0,0d0)
              UV(t,imat,jmat)=(0d0,0d0)
           end do
        end do
        do imat=1,NMAT
           do jmat=1,NMAT
              UdVd(t,imat,jmat)=UdVd(t,imat,jmat)+dconjg(U(tp1(t),idim,jmat,imat))*dconjg(V(jmat))
              UV(t,imat,jmat)=UV(t,imat,jmat)+U(t,idim,imat,jmat)*V(jmat)
           end do
        end do
        
        do imat=1,NMAT
           do jmat=1,NMAT
              do kmat=1,NMAT
                 DerivH_pre(t,idim,imat,jmat)=DerivH_pre(t,idim,imat,jmat)+UV(t,imat,kmat)*UdVd(t,kmat,jmat)
              end do
           end do
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  DerivH_pre=DerivH_pre*dcmplx(0.5d0/LatticeSpacing*dble(nmat))
  
  DerivH_pre=DerivH_pre*(0d0,-1d0)
  !add U-deriv and U^dag-deriv
  do t=1,nsite
     do idim=1,ndim
        do imat=1,NMAT
           do jmat=1,NMAT
              DerivH(t,idim,imat,jmat)=DerivH_pre(t,idim,imat,jmat)+dconjg(DerivH_pre(t,idim,jmat,imat))
           end do
        end do
     end do
  end do
  
  return

END SUBROUTINE MakeDerivHamiltonian_U

!***********************************************************
!***********************************************************
!   exp of a matrix : P -> Exp(i*P)
SUBROUTINE MATRIX_iEXP(MatSize,MAT,MATiEXP)

  implicit none 

  integer MatSize
  double complex MAT(1:MatSize,1:MatSize),MATiEXP(1:MatSize,1:MatSize)

  integer i,j,k
  double complex eig_iexp(1:MatSize)
  !lapack
  character jobvl,jobvr
  integer lda,ldvl,ldvr,info,lwork
  double complex w(1:MatSize),VL(1:MatSize),VR(1:MatSize,1:MATSIZE),&
       work(1:2*MatSize),rwork(1:2*MatSize)

  jobvl='N'
  jobvr='V'
  lda=MatSize
  ldvl=1
  ldvr=MatSize
  lwork=2*MatSize

  !MAT2=MAT

  call ZGEEV(jobvl,jobvr,MatSize,MAT,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork,info)
  ! in this case w(i) is real. 
  do i=1,MatSize
     !write(*,*)w(i)
     eig_iexp(i)=dcmplx(dcos(dble(w(i))))+dcmplx(dsin(dble(w(i))))*(0d0,1d0)
     !write(*,*)eig_log(i) 
 end do
  MATiEXP=(0d0,0d0)
  do i=1,MatSize
     do j=1,MatSize
        do k=1,MatSize
           MATiEXP(i,j)=MATiEXP(i,j)+VR(i,k)*eig_iexp(k)*dconjg(VR(j,k))
       end do
     end do
  end do

  return

END SUBROUTINE MATRIX_iEXP
!***********************************************************
!***********************************************************
!   exp of a matrix : P -> Exp(i*P)
SUBROUTINE MATRIX_iEXP_2(MatSize,MAT,MATiEXP,t,nsite)

  implicit none 

  integer MatSize, t,nsite
  double complex MAT(1:MatSize,1:MatSize,1:nsite),MATiEXP(1:MatSize,1:MatSize,1:nsite),&
       MAT_2(1:MatSize,1:MatSize)

  integer i,j,k
  double complex eig_iexp(1:MatSize)
  !lapack
  character jobvl,jobvr
  integer lda,ldvl,ldvr,info,lwork
  double complex w(1:MatSize),VL(1:MatSize),VR(1:MatSize,1:MATSIZE),&
       work(1:2*MatSize),rwork(1:2*MatSize)

  jobvl='N'
  jobvr='V'
  lda=MatSize
  ldvl=1
  ldvr=MatSize
  lwork=2*MatSize

  do i=1,MatSize
     do j=1,MatSize
        MAT_2(i,j)=MAT(i,j,t)
     end do
  end do
  call ZGEEV(jobvl,jobvr,MatSize,MAT_2,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork,info)
  ! in this case w(i) is real. 
  do i=1,MatSize
     !write(*,*)w(i)
     eig_iexp(i)=dcmplx(dcos(dble(w(i))))+dcmplx(dsin(dble(w(i))))*(0d0,1d0)
     !write(*,*)eig_log(i) 
  end do
  do i=1,MatSize
     do j=1,MatSize
        MATiEXP(i,j,t)=(0d0,0d0)
     end do
  end do
  do i=1,MatSize
     do j=1,MatSize
        do k=1,MatSize
           MATiEXP(i,j,t)=MATiEXP(i,j,t)+VR(i,k)*eig_iexp(k)*dconjg(VR(j,k))
       end do
     end do
  end do

  return

END SUBROUTINE MATRIX_iEXP_2
!***********************************************************
!***********************************************************
!    Polaakov loop 
SUBROUTINE CalcPol_no_phase(NMAT,NumSite,alpha,Pol)
  
  implicit none
  
  integer t,imat,jmat,kmat,p,NumSite,NMAT
  double precision alpha(1:NMAT)
  double complex Pol,MAT1(1:NMAT,1:NMAT),MAT2(1:NMAT,1:NMAT),eig(1:NMAT)
    double complex polyakov_line(1:NMAT)
 
  do imat=1,NMAT
     polyakov_line(imat)=dcmplx(dcos(alpha(imat)))+(0d0,1d0)*dcmplx(dsin(alpha(imat)))
  end do

  Pol=(0d0,0d0)
  do imat=1,nmat
     Pol=Pol+polyakov_line(imat)
  end do

  Pol=Pol/dcmplx(NMAT)
  
  return
  
END SUBROUTINE CalcPol_no_phase


!***********************************************************
!***********************************************************
!    Polaakov loop 
SUBROUTINE CalcPol_no_phase_2(NMAT,NumSite,alpha,Pol_dec,Pol_con,ndec)
  
  implicit none
  
  integer t,imat,jmat,kmat,p,NumSite,NMAT,ndec
  double precision alpha(1:nmat)
  double complex Pol_dec,pol_con,MAT1(1:NMAT,1:NMAT),MAT2(1:NMAT,1:NMAT),eig(1:NMAT)
  double complex polyakov_line(1:NMAT)
 
  do imat=1,NMAT
     polyakov_line(imat)=dcmplx(dcos(alpha(imat)))+(0d0,1d0)*dcmplx(dsin(alpha(imat)))
  end do

  Pol_dec=(0d0,0d0)
  do imat=1,ndec
     Pol_dec=Pol_dec+polyakov_line(imat)
  end do
  Pol_dec=Pol_dec/dcmplx(ndec)

  pol_con=(0d0,0d0)
  if(ndec.LT.nmat)then
     do imat=ndec+1,nmat
        Pol_con=Pol_con+polyakov_line(imat)
     end do
     Pol_con=Pol_con/dcmplx(nmat-ndec)
  end if
  !pol_con=0 if ndec=nmat
  
  return
  
END SUBROUTINE CalcPol_no_phase_2

!***********************************************************
!***********************************************************
!     total kinetic energy - (zero-point energy)

SUBROUTINE total_energy(tot_ene,U,alpha,NMAT,DIM,NumSite,LatticeSpacing) 

  implicit none

  integer t,tp1,i,j,k,p,q,r,NMAT,NumSite,DIM,idim,jdim,imat,jmat,kmat
  doublecomplex U(1:NumSite,1:DIM,1:NMAT,1:NMAT),temp(1:nmat,1:nmat),temp2(1:nmat,1:nmat)
  doubleprecision alpha(1:NMAT)

  
  doubleprecision kin,LatticeSpacing,tot_kin,tot_ene


  call kinetic(tot_kin,U,alpha,NMAT,DIM,NumSite,LatticeSpacing) 
 
  tot_ene=tot_kin
  !** no potential energy in this case **
  
  return

END SUBROUTINE total_energy

!##############################################
subroutine permutation(nmat,ndim,nsite,U,alpha,ndec,g_coeff,P_fix,P_fix_width,perm_count)

  use mtmod !Mersenne twistor
  implicit none

  integer nmat,ndim,nsite,ndec
  doublecomplex U(1:nsite,1:NDIM,1:NMAT,1:NMAT),temp_u(1:NMAT),polyakov_line_1,polyakov_line_2,pol_dec,pol_con
  doubleprecision alpha(1:NMAT),ransu1,ransu2,temp_alpha,g_coeff,P_fix,P_fix_width,hamstart,hamend,hantei
  
  integer idim
  integer imat,jmat,kmat
  integer isite,jsite
  integer perm_count
  
  !call random_number(ransu1)
  !call random_number(ransu2)
  ransu1=grnd()
  ransu2=grnd()
  imat=int(ransu1*dble(ndec))+1
  jmat=int(ransu2*dble(nmat-ndec))+1+ndec

  polyakov_line_1=dcmplx(dcos(alpha(imat)))+(0d0,1d0)*dcmplx(dsin(alpha(imat)))
  polyakov_line_2=dcmplx(dcos(alpha(jmat)))+(0d0,1d0)*dcmplx(dsin(alpha(jmat)))

  call CalcPol_no_phase_2(NMAT,NSite,alpha,Pol_dec,Pol_con,ndec)
  HamStart=0d0
  if(abs(pol_dec).GT.P_fix+P_fix_width)then
     HamStart=HamStart+0.5d0*g_coeff*(abs(pol_dec)-P_fix-P_fix_width)**2d0
  else if(abs(pol_dec).LT.P_fix-P_fix_width)then
     HamStart=HamStart+0.5d0*g_coeff*(abs(pol_dec)-P_fix+P_fix_width)**2d0
  end if
  if(abs(pol_con).GT.P_fix_width)then
     HamStart=HamStart+0.5d0*g_coeff*(abs(pol_con)-P_fix_width)**2d0
  end if

  Pol_dec=Pol_dec+(polyakov_line_2-polyakov_line_1)/dcmplx(ndec)
  Pol_con=Pol_con+(polyakov_line_1-polyakov_line_2)/dcmplx(nmat-ndec)

  HamEnd=0d0
  if(abs(pol_dec).GT.P_fix+P_fix_width)then
     HamEnd=HamEnd+0.5d0*g_coeff*(abs(pol_dec)-P_fix-P_fix_width)**2d0
  else if(abs(pol_dec).LT.P_fix-P_fix_width)then
     HamEnd=HamEnd+0.5d0*g_coeff*(abs(pol_dec)-P_fix+P_fix_width)**2d0
  end if
  if(abs(pol_con).GT.P_fix_width)then
     HamEnd=HamEnd+0.5d0*g_coeff*(abs(pol_con)-P_fix_width)**2d0
  end if


  !call random_number(hantei)
  hantei=grnd()
  if(dexp(HamStart-HamEnd) > hantei)THEN
     perm_count=perm_count+1
     !permute alpha(imat) and alpha(jmat)
     temp_alpha=alpha(imat)
     alpha(imat)=alpha(jmat)
     alpha(jmat)=temp_alpha
     !permute imat-th column and jmat-th column
     do idim=1,ndim
        do isite=1,nsite
           do kmat=1,nmat
              temp_u(kmat)=U(isite,idim,kmat,imat)
           end do
           do kmat=1,nmat
              u(isite,idim,kmat,imat)=u(isite,idim,kmat,jmat)
           end do
           do kmat=1,nmat
              u(isite,idim,kmat,jmat)=temp_u(kmat)
           end do
        end do
     end do
     !permute imat-th row and jmat-th row
     do idim=1,ndim
        do isite=1,nsite
           do kmat=1,nmat
              temp_u(kmat)=U(isite,idim,imat,kmat)
           end do
           do kmat=1,nmat
              u(isite,idim,imat,kmat)=u(isite,idim,jmat,kmat)
           end do
           do kmat=1,nmat
              u(isite,idim,jmat,kmat)=temp_u(kmat)
           end do
        end do
     end do
     
     
  end if
  
  
  return

END subroutine Permutation
!***********************************************************
!***********************************************************
!    Polyakov loop 
SUBROUTINE CalcPol(NMAT,alpha,Pol,phase)
  
  implicit none
  
  integer t,imat,NMAT
  double precision alpha(1:NMAT)
  double complex Pol
  double precision pi,phase(1:nmat)
  double complex overall_phase
  double complex polyakov_line(1:NMAT)
 
  do imat=1,NMAT
     polyakov_line(imat)=dcmplx(dcos(alpha(imat)))+(0d0,1d0)*dcmplx(dsin(alpha(imat)))
  end do
  
  pi=2d0*dasin(1d0)
  Pol=(0d0,0d0)
  do imat=1,NMAT
     pol=pol+polyakov_line(imat)
  end do

  Pol=Pol/dcmplx(NMAT)
  overall_phase=pol/dcmplx(abs(pol))
  !write(*,*)overall_phase
  polyakov_line=polyakov_line/overall_phase
  do imat=1,nmat
     if(dble(polyakov_line(imat)).GE.0d0)then
        phase(imat)=dasin(dble((0d0,-1d0)*polyakov_line(imat)))
     else if(dble((0d0,-1d0)*polyakov_line(imat)).GT.0d0)then
        phase(imat)=pi-dasin(dble((0d0,-1d0)*polyakov_line(imat)))
     else
        phase(imat)=-pi-dasin(dble((0d0,-1d0)*polyakov_line(imat)))
     end if
  end do

  
  return
  
END SUBROUTINE CalcPol


!***********************************************************
!***********************************************************
!    Wilsobn loop analogous to 2-pt Polyakov loop correlator
SUBROUTINE Calc_Pol_Correlator(NMAT,nsite,ndim,alpha,cor_dec,cor_con,cor_mix,U,ndec,nprod)
  
  implicit none
  
  integer t,imat,jmat,kmat,NMAT,ndec,ndim,idim,iprod,nprod,isite,nsite
  double precision alpha(1:NMAT)
  double precision cor_dec(1:nprod),cor_con(1:nprod),cor_mix(1:nprod)
  doublecomplex U(1:nsite,1:NDIM,1:NMAT,1:NMAT)
  double complex prod_U(1:nmat,1:nmat),temp(1:nmat,1:nmat)
  double complex Pol
  double complex overall_phase
  double complex polyakov_line(1:NMAT)
  
  do imat=1,NMAT
     polyakov_line(imat)=dcmplx(dcos(alpha(imat)))+(0d0,1d0)*dcmplx(dsin(alpha(imat)))
  end do
  
  cor_dec=0d0
  cor_con=0d0
  cor_mix=0d0
  do idim=1,ndim
     do isite=1,nsite
        do imat=1,nmat
           do jmat=1,nmat
              prod_U(imat,jmat)=U(isite,idim,imat,jmat)
           end do
        end do
        do iprod=1,nprod
           do imat=1,ndec
              do jmat=1,ndec
                 cor_dec(iprod)=cor_dec(iprod)+dble(polyakov_line(imat)*prod_u(imat,jmat)&
                      *dconjg(polyakov_line(jmat))*dconjg(prod_u(imat,jmat)))
              end do
           end do
           do imat=ndec+1,nmat
              do jmat=ndec+1,nmat
                 cor_con(iprod)=cor_con(iprod)+dble(polyakov_line(imat)*prod_u(imat,jmat)&
                      *dconjg(polyakov_line(jmat))*dconjg(prod_u(imat,jmat)))
              end do
           end do

           do imat=1,ndec
              do jmat=ndec+1,nmat
                 cor_mix(iprod)=cor_mix(iprod)+dble(polyakov_line(imat)*prod_u(imat,jmat)&
                      *dconjg(polyakov_line(jmat))*dconjg(prod_u(imat,jmat)))
              end do
           end do
           
           if(iprod.LT.nprod)then
              temp=prod_u
              prod_u=(0d0,0d0)
              do imat=1,nmat
                 do jmat=1,nmat
                    do kmat=1,nmat
                       prod_U(imat,jmat)=prod_U(imat,jmat)+temp(imat,kmat)*U(isite,idim,kmat,jmat)
                    end do
                 end do
              end do
           end if
        end do
     end do
  end do
  cor_dec=cor_dec/(ndec*nsite*ndim)
  cor_con=cor_con/((nmat-ndec)*nsite*ndim)
  cor_mix=cor_mix/(nsite*ndim)
  
  return
  
END SUBROUTINE Calc_Pol_Correlator
