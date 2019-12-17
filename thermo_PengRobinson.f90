
!!===========================================================================================================
!!
!! Module: mod_eosPR
!!
!! Thermodynamic calculations using Peng Robinson equation of state
!!
!! Use module name: mod_eosPR
!!
!!<> Notes:
!!    - all variable have SI units (e.g. Pressure [Pa], Density [kg/m3], Molar Enthalpy [J/kmol])
!!
!!============================================================================================================
module mod_eosPR      
  implicit none
  private
  
  ! <--- Do not modify the following parameters --->
  integer, parameter :: PREC    = 8 ! Precision of internals (routines, ODE solver etc) 
  real(PREC), parameter :: T_STD = 298.15_PREC, R_UNIVERSAL = 8314.4720_PREC  
  
  ! Provided interface
    public :: eosPR
    public :: eosPR_init,eosPR_nspecies
    public :: eosPR_pressure
    public :: eosPR_molarVolume, eosPR_molarVolume2
    public :: eosPR_enthalpy, eosPR_enthalpy2
    public :: eosPR_temperature_from_e
    public :: eosPR_temperature_from_EOS
    public :: eosPR_density, eosPR_density_from_EOS
    public :: eosPR_alpha_D
    public :: eosPR_IntCp
    public :: eosPR_energy
    public :: eosPR_mu, eosPR_mu_R
    public :: eosPR_Cp
    public :: eosPR_spSound
    public :: eosPR_heatFlux_BK,eosPR_massFlux_BK,eosPR_heatFlux_IK,eosPR_massFlux_IK
    public :: eosPR_molecularMass
    public :: eosPR_AmBm
    public :: eosPR_dAmdx
    public :: eosPR_K2
    public :: eosPR_mixMolecularMass
    public :: eosPR_moleFrac
    public :: eosPR_par_h
    public :: eosPR_par_molarVolume
    public :: eosPR_parp_parv_TX
    
  ! supercritical thermo Class
  type eosPR
     private
     ! Species specifications
     integer :: nspecies_
     real(PREC), pointer :: molecularMass_(:) !molecular mass [kg/kmol]
     real(PREC), pointer :: omega_(:),vc_(:),Tc_(:),pc_(:),zc_(:)
     real(PREC) :: D_L_,D_U_,T_L_,T_U_
     integer, pointer :: speciesIndex_(:)
     character(16), pointer :: speciesName_(:) 
  end type eosPR

CONTAINS

  !!============================================================================
  !!
  !! Sub: eosPR_init
  !!    Construct a thermo object 
  !!
  !!  init(eos, gas, chem, name)
  !!
  !!  spt [inout, type(eos)] - The object to construct.
  !=============================================================================
  subroutine eosPR_init(this) !eosinp)
    type(eosPR), intent(inout) :: this    
!     character(*), intent(in) :: cheminp
    !
    ! Temporary: Number of Species specifications: This part will be supplied by a reader in the future
    !    
    this%nspecies_ = 2    ! number of species
    !
    ! Allocations
    !    
    allocate( this%molecularMass_(this%nspecies_))
    allocate( this%omega_(this%nspecies_))
    allocate( this%vc_(this%nspecies_))
    allocate( this%Tc_(this%nspecies_))
    allocate( this%Zc_(this%nspecies_))      
    allocate( this%pc_(this%nspecies_))   
    allocate( this%speciesIndex_(this%nspecies_))
    allocate( this%speciesName_(this%nspecies_))
                                
    !
    ! Temporary: Species specifications. This part will be supplied by a reader in the future
    !
    this%speciesName_(1)="N2"
    this%speciesIndex_(1)=1
    this%omega_(1)=0.039_PREC
    this%vc_(1)=0.0898_PREC     !m3/kmol
    this%Tc_(1)=126.26_PREC     !K
    this%pc_(1)=33.55_PREC*101325._PREC      !Pa
    this%zc_(1)=0.29_PREC 
    this%molecularMass_(1)=28.013_PREC !Molar Mass  kg/kmol
    
    this%speciesName_(2)="C7H16"
    this%speciesIndex_(2)=2
    this%omega_(2)=0.349_PREC
    this%vc_(2)=0.432_PREC     ![m3/kmol]
    this%Tc_(2)=540.3_PREC     ![K]
    this%pc_(2)=27.04_PREC*101325._PREC
    this%zc_(2)=0.263_PREC  
    this%molecularMass_(2)=100.205_PREC !Molar Mass  kg/kmol  

    this%D_L_=8._PREC    ![kg/m3]
    this%D_U_=248._PREC  ![kg/m3]
    this%T_L_=600._PREC  ![K]
    this%T_U_=1100._PREC ![K]

  end subroutine eosPR_init
  
  !!============================================================================
  !!
  !! Sub: eosPR_nspecies
  !!    Returns number of species
  !!
  !! eosPR_nspecies(this)
  !!
  !!  this [type(eos), inout] - eosPR object
  !!  nsp [integer,out] - number of species
  !===========================================================================
  function eosPR_nspecies(this) result(nsp)
    implicit none
    type(eosPR), intent(inout) :: this
    integer :: nsp
    !
    nsp=this%nspecies_
  end function eosPR_nspecies
  
  !!============================================================================
  !!
  !! Sub: eosPR_pressure
  !!    Returns the pressure [Eq. (2.6)]
  !!
  !! eosPR_pressure(this,T,D,massFrac,p)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  T [real, in] - Temperature [K]
  !!  D [real, in] - Density [kg/m3]  
  !!  p [real, out] - Pressure [Pa]    
  !===========================================================================  
  subroutine eosPR_pressure(this,T,D,massFrac,p)
    implicit none
    type(eosPR), intent(inout) :: this 
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: T,D 
    real(PREC), intent(out) :: p
    real(PREC) :: molarVolume
    real(PREC) :: moleFrac(this%nspecies_)
    real(PREC) :: Bm,Am
    real(PREC) :: A(this%nspecies_,this%nspecies_)
    real(PREC) :: B(this%nspecies_) 
    !
    call eosPR_molarVolume(this,D,massFrac,molarVolume)
    call eosPR_moleFrac(this,massFrac,moleFrac)
    call eosPR_AmBm(this,T,moleFrac,Am,Bm,A,B)
    p=R_UNIVERSAL*T/(molarVolume-Bm)-Am/(molarVolume**2+2*molarVolume*Bm-Bm**2) !Pa

  end subroutine eosPR_pressure
 
  !!============================================================================
  !!
  !! Sub: eosPR_ABC
  !!    Returns A, B, C, Tc, vc, pc [Appendix A]
  !!
  !! eosPR_ABC(this,T,moleFrac,A,B,C,Tc,Vc,pc)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  moleFrac [real(1:nspecies), in] - mole fraction of all chemical species
  !!  T [real, in] - Temperature [K]
  !!  A, B, C, Tc, vc, pc
  !===========================================================================  
  subroutine eosPR_ABC(this,T,moleFrac,A,B,C,Tc,vc,pc)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: moleFrac(this%nspecies_)
    real(PREC), intent(in) :: T
    real(PREC), intent(out) :: A(this%nspecies_,this%nspecies_)
    real(PREC), intent(out) :: B(this%nspecies_)
    real(PREC), intent(out) :: C(this%nspecies_,this%nspecies_)
    real(PREC), intent(out) :: Tc(this%nspecies_,this%nspecies_)
    real(PREC), intent(out) :: vc(this%nspecies_,this%nspecies_)
    real(PREC), intent(out) :: pc(this%nspecies_,this%nspecies_)
    real(PREC) :: omega(this%nspecies_,this%nspecies_)
    real(PREC) :: zc(this%nspecies_,this%nspecies_)
    real(PREC) :: k(this%nspecies_,this%nspecies_)
    integer :: i,j
    ! 
    omega(1,1)=this%omega_(1)
    omega(2,2)=this%omega_(2)
    omega(1,2)=0.5_PREC*(this%omega_(1)+this%omega_(2)) 
    omega(2,1)=omega(1,2)

    zc(1,1)=this%zc_(1)
    zc(2,2)=this%zc_(2)
    zc(1,2)=0.5_PREC*(this%zc_(1)+this%zc_(2))
    zc(2,1)=zc(1,2)
    
    k(1,1)=0._PREC
    k(2,2)=0._PREC
    k(1,2)=0.1_PREC
    k(2,1)=0.1_PREC

    Tc(1,1)=this%Tc_(1)
    Tc(2,2)=this%Tc_(2)
    Tc(1,2)=sqrt(this%Tc_(1)*this%Tc_(2))*(1._PREC-k(1,2))
    Tc(2,1)=Tc(1,2)

    vc(1,1)=this%vc_(1)
    vc(2,2)=this%vc_(2)
    vc(1,2)=(this%vc_(1)**(1._PREC/3._PREC)+this%vc_(2)**(1._PREC/3._PREC))**3/8._PREC
    vc(2,1)=vc(1,2)

    pc(1,1)=this%pc_(1)
    pc(2,2)=this%pc_(2)
    pc(1,2)=zc(1,2)*R_UNIVERSAL*Tc(1,2)/vc(1,2)
    pc(2,1)=pc(1,2)

    do i=1,this%nspecies_
       do j=1,this%nspecies_
          C(i,j)=0.37464_PREC+1.54226_PREC*omega(i,j)-0.26992_PREC*omega(i,j)**2
       enddo
    enddo
    
    do i=1,this%nspecies_
       B(i)=0.077796_PREC*R_UNIVERSAL*this%Tc_(i)/this%pc_(i)
    enddo 

    do i=1,this%nspecies_
       do j=1,this%nspecies_
          A(i,j)=0.457236_PREC*(R_UNIVERSAL*Tc(i,j))**2*&
                 (1._PREC+C(i,j)*(1._PREC-sqrt(T/Tc(i,j))))**2/pc(i,j)
       enddo
    enddo

  end subroutine eosPR_ABC

  !!============================================================================
  !!
  !! Sub: eosPR_AmBm 
  !!    Returns Am, Bm, A, B [Eq. (2.7)]
  !!
  !! eosPR_AmBm(this,T,moleFrac,Am,Bm,A,B)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  moleFrac [real(1:nspecies), in] - mole fraction of all chemical species
  !!  T [real, in] - Temperature [K]
  !!  Am [real, out] - [(m3/kmol)2]
  !!  Bm [real, out] - [m3/kmol] 
  !!  A [real(1:nspecies,1:nspecies), out] - [(m3/kmol)2]
  !!  B [real(1:nspecies), out] - [m3/kmol]
  !===========================================================================  
  subroutine eosPR_AmBm(this,T,moleFrac,Am,Bm,A,B)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: moleFrac(this%nspecies_)
    real(PREC), intent(in) :: T
    real(PREC), intent(out) :: Am,Bm
    real(PREC), intent(out) :: A(this%nspecies_,this%nspecies_)
    real(PREC), intent(out) :: B(this%nspecies_)
    real(PREC) :: C(this%nspecies_,this%nspecies_)
    real(PREC) :: Tc(this%nspecies_,this%nspecies_)
    real(PREC) :: vc(this%nspecies_,this%nspecies_)
    real(PREC) :: pc(this%nspecies_,this%nspecies_)
    integer :: i,j
    ! 
    call eosPR_ABC(this,T,moleFrac,A,B,C,Tc,vc,pc)
    
    Am=0._PREC
    Bm=0._PREC

    do j=1,this%nspecies_
       Bm=moleFrac(j)*B(j)+Bm
    enddo

    do i=1,this%nspecies_
       do j=1,this%nspecies_
          Am=moleFrac(i)*moleFrac(j)*A(i,j)+Am
       enddo
    enddo 

  end subroutine eosPR_AmBm

  !!============================================================================
  !!
  !! Sub: eosPR_AmDerivatives
  !!    Returns derivatives of Am [Appendix B]
  !!
  !! eosPR_AmDerivatives(this,T,moleFrac,parAm_parXj,parAm_parT,
  !!                     par2Am_parXj_parT,par2Am_parT2)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  moleFrac [real(1:nspecies), in] - mole fraction of all chemical species
  !!  T [real, in] - Temperature [K]
  !!  moleFrac,parAm_parXj,parAm_parT,par2Am_parXj_parT,par2Am_parT2 [real, out]
  !===========================================================================  
  subroutine eosPR_AmDerivatives(this,T,moleFrac,parAm_parXj,parAm_parT,&
                                  par2Am_parXj_parT,par2Am_parT2)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: moleFrac(this%nspecies_)
    real(PREC), intent(in) :: T
    real(PREC), intent(out) :: parAm_parXj(this%nspecies_)
    real(PREC), intent(out) :: parAm_parT
    real(PREC), intent(out) :: par2Am_parXj_parT(this%nspecies_)
    real(PREC), intent(out) :: par2Am_parT2
    real(PREC) :: A(this%nspecies_,this%nspecies_)
    real(PREC) :: B(this%nspecies_)
    real(PREC) :: C(this%nspecies_,this%nspecies_)
    real(PREC) :: G(this%nspecies_,this%nspecies_)
    real(PREC) :: Tc(this%nspecies_,this%nspecies_)
    real(PREC) :: vc(this%nspecies_,this%nspecies_)
    real(PREC) :: pc(this%nspecies_,this%nspecies_)
    integer :: i,j
    ! 
    call eosPR_ABC(this,T,moleFrac,A,B,C,Tc,vc,pc)

    do i=1,this%nspecies_
       do j=1,this%nspecies_
          G(i,j)=C(i,j)*sqrt(T/Tc(i,j))/(1._PREC+C(i,j)*(1._PREC-sqrt(T/Tc(i,j))))
       enddo
    enddo
    
    parAm_parXj(1:this%nspecies_)=0._PREC
    do i=1,this%nspecies_
       do j=1,this%nspecies_
          parAm_parXj(i)=A(i,j)*moleFrac(j)+parAm_parXj(i)
       enddo
       parAm_parXj(i)=2._PREC*parAm_parXj(i)
    enddo

    parAm_parT=0._PREC
    do i=1,this%nspecies_
       do j=1,this%nspecies_
          parAm_parT=moleFrac(i)*moleFrac(j)*A(i,j)*G(i,j)+parAm_parT
       enddo
    enddo
    parAm_parT=-1._PREC*parAm_parT/T

    par2Am_parXj_parT(1:this%nspecies_)=0._PREC
    do i=1,this%nspecies_
       do j=1,this%nspecies_
          par2Am_parXj_parT(i)=moleFrac(j)*A(i,j)*G(i,j)+par2Am_parXj_parT(i)
       enddo
       par2Am_parXj_parT(i)=-2._PREC*par2Am_parXj_parT(i)/T
    enddo 

    par2Am_parT2=0._PREC
    do i=1,this%nspecies_
       do j=1,this%nspecies_
          par2Am_parT2=moleFrac(i)*moleFrac(j)*C(i,j)*(1+C(i,j))*Tc(i,j)/pc(i,j)*sqrt(Tc(i,j)/T)&
                      +par2Am_parT2
       enddo
    enddo
    par2Am_parT2=0.457236_PREC*R_UNIVERSAL**2*par2Am_parT2/(2._PREC*T)

  end subroutine eosPR_AmDerivatives

  !!============================================================================
  !!
  !! Sub: eosPR_parp_parv_TX
  !!    Returns (partial p)/(partial v) at constant T and X [Eq. (2.12)]
  !!
  !! eosPR_parp_parv_TX(this,T,D,massFrac,parp_parv_TX)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  T [real, in] - Temperature [K]
  !!  D [real, in] - Density [kg/m3]  
  !!  parp_parv_TX [real, out] - (partial p)/(partial v) at constant T and X
  !!                             [Pa.kmol/m3]    
  !===========================================================================  
  subroutine eosPR_parp_parv_TX(this,T,D,massFrac,parp_parv_TX)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: T,D
    real(PREC), intent(out) :: parp_parv_TX
    real(PREC) :: molarVolume
    real(PREC) :: moleFrac(this%nspecies_)
    real(PREC) :: Bm,Am
    real(PREC) :: A(this%nspecies_,this%nspecies_)
    real(PREC) :: B(this%nspecies_)
    !
    call eosPR_molarVolume(this,D,massFrac,molarVolume)
    call eosPR_moleFrac(this,massFrac,moleFrac)
    call eosPR_AmBm(this,T,moleFrac,Am,Bm,A,B)
    parp_parv_TX=-R_UNIVERSAL*T/(molarVolume-Bm)**2*(1._PREC-2._PREC*Am/&
                 (R_UNIVERSAL*T*(molarVolume+Bm)*(molarVolume/(molarVolume-Bm)+Bm/(molarVolume+Bm))**2))

  end subroutine eosPR_parp_parv_TX

  !!============================================================================
  !!
  !! Sub: eosPR_par_molarVolume
  !!    Returns (partial v)/(partial X_j) [Eq. (2.22)]
  !!
  !! eosPR_par_molarVolume(this,T,D,massFrac,par_molarVolume)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  T [real, in] - Temperature [K]
  !!  D [real, in] - Density [kg/m3]  
  !!  par_molarVolume [real(1:nspecies), out] - (partial v)/(partial X_j) [m3/kmol]    
  !===========================================================================  
  subroutine eosPR_par_molarVolume(this,T,D,massFrac,par_molarVolume)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: T,D
    real(PREC), intent(out) :: par_molarVolume(1:this%nspecies_)
    real(PREC) :: molarVolume
    real(PREC) :: moleFrac(this%nspecies_)
    real(PREC) :: Bm,Am
    real(PREC) :: A(this%nspecies_,this%nspecies_)
    real(PREC) :: B(this%nspecies_)
    real(PREC) :: parp_parv_TX
    real(PREC) :: aux
    integer :: i,j
    !
    call eosPR_molarVolume(this,D,massFrac,molarVolume)
    call eosPR_moleFrac(this,massFrac,moleFrac)
    call eosPR_AmBm(this,T,moleFrac,Am,Bm,A,B)
    call eosPR_parp_parv_TX(this,T,D,massFrac,parp_parv_TX)
    do j=1,this%nspecies_
       aux = 0._PREC
       do i=1,this%nspecies_
          aux=A(j,i)*moleFrac(i)+aux
       enddo
       par_molarVolume(j)=-1/parp_parv_TX*( &
                          R_UNIVERSAL*T/(molarVolume-Bm)+R_UNIVERSAL*T*B(j)/(molarVolume-Bm)**2 &
                        + 2._PREC*Am*(molarVolume-Bm)*B(j)/(molarVolume**2+2._PREC*molarVolume*Bm-Bm**2)**2 &
                        - 2._PREC*aux/(molarVolume**2+2._PREC*molarVolume*Bm-Bm**2))
    enddo

  end subroutine eosPR_par_molarVolume

  !!============================================================================
  !!
  !! Sub: eosPR_par_h 
  !!    Returns (partial h)/(partial X_j) [Eq. (2.23)]
  !!
  !! eosPR_par_h(this,T,D,massFrac,par_h)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  T [real, in] - Temperature [K]
  !!  D [real, in] - Density [kg/m3]  
  !!  par_h [real(1:nspecies), out] - (partial h)/(partial X_j) [J/kmol]    
  !===========================================================================  
  subroutine eosPR_par_h(this,T,D,massFrac,par_h)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: T,D
    real(PREC), intent(out) :: par_h(this%nspecies_)
    real(PREC) :: molarVolume
    real(PREC) :: moleFrac(this%nspecies_)
    real(PREC) :: Bm,Am
    real(PREC) :: A(this%nspecies_,this%nspecies_)
    real(PREC) :: B(this%nspecies_)
    real(PREC) :: parAm_parXj(this%nspecies_)
    real(PREC) :: parAm_parT
    real(PREC) :: par2Am_parXj_parT(this%nspecies_)
    real(PREC) :: par2Am_parT2
    real(PREC) :: par_molarVolume(this%nspecies_)
    real(PREC) :: h0(this%nspecies_)
    real(PREC) :: p
    real(PREC) :: K1
    integer :: i,j
    !
    call eosPR_molarVolume(this,D,massFrac,molarVolume)
    call eosPR_moleFrac(this,massFrac,moleFrac)
    call eosPR_AmBm(this,T,moleFrac,Am,Bm,A,B)
    call eosPR_AmDerivatives(this,T,moleFrac,parAm_parXj,parAm_parT,&
                                             par2Am_parXj_parT,par2Am_parT2)
    call eosPR_par_molarVolume(this,T,D,massFrac,par_molarVolume)
    call eosPR_pressure(this,T,D,massFrac,p)
    call eosPR_K1_2(this,molarVolume,Bm,K1)
    call eosPR_h0_species(this,T,h0)
    
    do i=1,this%nspecies_
       par_h(i)=h0(i)+p*par_molarVolume(i)-R_UNIVERSAL*T&
                         + (Am-T*parAm_parT)*(par_molarVolume(i)-molarVolume*B(i)/Bm)/&
                           (molarVolume**2+2._PREC*molarVolume*Bm-Bm**2)&
                         + K1*(parAm_parXj(i)-T*par2Am_parXj_parT(i)&
                         - (Am-T*parAm_parT)*B(i)/Bm) 
    enddo
 
  end subroutine eosPR_par_h

  !!============================================================================
  !! (ONLY FOR TEST)
  !! Sub: eosPR_par_h_2
  !!    Returns (partial h)/(partial X_j) [Eq. (2.23)]
  !!
  !! eosPR_par_h_2(this,T,D,massFrac,par_h)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  T [real, in] - Temperature [K]
  !!  D [real, in] - Density [kg/m3]  
  !!  par_h [real(1:nspecies), out] - (partial h)/(partial X_j) [J/kmol]    
  !=========================================================================== 
  subroutine eosPR_par_h_2(this,T,D,massFrac,par_h)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: T,D
    real(PREC), intent(out) :: par_h(this%nspecies_)
    real(PREC) :: molarVolume
    real(PREC) :: moleFrac(this%nspecies_)
    real(PREC) :: Bm,Am
    real(PREC) :: A(this%nspecies_,this%nspecies_)
    real(PREC) :: B(this%nspecies_)
    real(PREC) :: parAm_parXj(this%nspecies_)
    real(PREC) :: parAm_parT
    real(PREC) :: par2Am_parXj_parT(this%nspecies_)
    real(PREC) :: par2Am_parT2
    real(PREC) :: par_molarVolume(this%nspecies_)
    real(PREC) :: h0(this%nspecies_)
    real(PREC) :: p
    real(PREC) :: K1
    integer :: i,j
    !
    call eosPR_molarVolume(this,D,massFrac,molarVolume)
    call eosPR_moleFrac(this,massFrac,moleFrac)
    call eosPR_AmBm(this,T,moleFrac,Am,Bm,A,B)
    call eosPR_Am(this,T,moleFrac,Am)
    call eosPR_Bm(this,T,moleFrac,Bm)

    call eosPR_dAmdT(this,T,moleFrac,parAm_parT)
    call eosPR_dAmdX(this,T,moleFrac,parAm_parXj)
    call eosPR_d2AmdXdT(this,T,moleFrac,par2Am_parXj_parT)
    call eosPR_par_molarVolume(this,T,D,massFrac,par_molarVolume)
    call eosPR_pressure(this,T,D,massFrac,p)
    call eosPR_K1(this,T,D,massFrac,K1)
    call eosPR_h0_species(this,T,h0)

    do i=1,this%nspecies_
       par_h(i)=h0(i)+p*par_molarVolume(i)-R_UNIVERSAL*T&
                         + (Am-T*parAm_parT)*(par_molarVolume(i)-molarVolume*B(i)/Bm)/&
                           (molarVolume**2+2._PREC*molarVolume*Bm-Bm**2)&
                         + K1*(parAm_parXj(i)-T*par2Am_parXj_parT(i)&
                         - (Am-T*parAm_parT)*B(i)/Bm)
    enddo

  end subroutine eosPR_par_h_2
  !!============================================================================
  !!
  !! Sub: eosPR_K1_2
  !!    Returns K1 [Eq. (2.13)]
  !!
  !! eosPR_K1_2(this,molarVolume,Bm,K1)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  molarVolume [real, in] - molar specific volume [m3/kmol]
  !!  Bm [real, in]
  !!  K1 [real, out]  
  !===========================================================================  
  subroutine eosPR_K1_2(this,molarVolume,Bm,K1)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: molarVolume
    real(PREC), intent(in) :: Bm
    real(PREC), intent(out) :: K1 
    !
    K1=log((molarVolume+(1._PREC-sqrt(2._PREC))*Bm)&
          /(molarVolume+(1._PREC+sqrt(2._PREC))*Bm))
    K1=K1/(2._PREC*sqrt(2._PREC)*Bm)

  end subroutine eosPR_K1_2

  !!============================================================================
  !!
  !! Sub: eosPR_K2
  !!    Returns K2 [Eq. (2.26)]
  !!
  !! eosPR_K2(this,T,molarVolume,A,B,Bm,K2)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  T [real, in] - temperature [K]
  !!  molarVolume [real, in] - molar specific volume [m3/kmol]
  !!  A [real(1:nspecies,1:nspecies), in] - [(m3/kmol)2]
  !!  B [real(1:nspecies), in] - [m3/kmol]
  !!  Bm [real, in] - [m3/kmol]
  !!  K2 [real, out]
  !===========================================================================  
  subroutine eosPR_K2(this,T,molarVolume,A,B,Bm,K2)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: T,molarVolume
    real(PREC), intent(in) :: A(this%nspecies_,this%nspecies_)
    real(PREC), intent(in) :: B(this%nspecies_)
    real(PREC), intent(in) :: Bm
    real(PREC), intent(out) :: K2
    real(PREC) :: K1
    !
    call eosPR_K1_2(this,molarVolume,Bm,K1)  
    K2=((B(2)-B(1))/Bm)**2+2._PREC/(R_UNIVERSAL*T)*K1/Bm**2&
      *(A(2,2)*B(1)**2+A(1,1)*B(2)**2-2._PREC*A(1,2)*B(1)*B(2))

  end subroutine eosPR_K2

  !!============================================================================
  !!
  !! Sub: eosPR_alpha_D
  !!    Returns mass diffusion factor [Eq. (2.25)]
  !!
  !! eosPR_alpha_D(this,T,D,massFrac,alpha_D)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  T [real, in] - Temperature [K]
  !!  D [real, in] - Density [kg/m3]  
  !!  alpha_D [real, out] - mass diffusion factor
  !===========================================================================  
  subroutine eosPR_alpha_D(this,T,D,massFrac,alpha_D)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: T,D
    real(PREC), intent(out) :: alpha_D
    real(PREC) :: v
    real(PREC) :: moleFrac(this%nspecies_)
    real(PREC) :: Bm,Am
    real(PREC) :: A(this%nspecies_,this%nspecies_)
    real(PREC) :: B(this%nspecies_) 
    real(PREC) :: dAmdX(this%nspecies_)
    real(PREC) :: parp_parv_TX
    real(PREC) :: K2
    !
    call eosPR_molarVolume(this,D,massFrac,v)
    call eosPR_moleFrac(this,massFrac,moleFrac)
    call eosPR_AmBm(this,T,moleFrac,Am,Bm,A,B)
    call eosPR_dAmdX(this,T,moleFrac,dAmdX)
    call eosPR_parp_parv_TX(this,T,D,massFrac,parp_parv_TX)
    call eosPR_K2(this,T,v,A,B,Bm,K2)

    alpha_D=1._PREC+moleFrac(1)*moleFrac(2)*&
            (R_UNIVERSAL*T/(parp_parv_TX*Bm**2)*&
            ((B(2)-B(1))/(v-Bm)+&
            (B(1)*dAmdX(2)-B(2)*dAmdX(1))/&
            (R_UNIVERSAL*T*(v**2+2._PREC*v*Bm-Bm**2)))**2+K2)

  end subroutine eosPR_alpha_D

  !!============================================================================
  !!
  !! Sub: eosPR_mu
  !!    Returns mu [Eq. (2.27)]
  !!
  !! eosPR_mu(this,T0,mu_R,T,mu)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  T0 [real, in] - reference temperature for mixing layer simulations [K]
  !!  mu_R [real, in] - reference dynamic viscosity [kg/m.s]
  !!  T [real, in] - temperature [K]
  !!  mu [real, out] - dynamic viscosity [kg/m.s]
  !===========================================================================  
  subroutine eosPR_mu(this,T0,mu_R,T,mu)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: T0
    real(PREC), intent(in) :: mu_R
    real(PREC), intent(in) :: T
    real(PREC), intent(out) :: mu
    !
    mu=mu_R*(T/T0)**0.7_PREC

  end subroutine eosPR_mu

  !!============================================================================
  !!
  !! Sub: eosPR_mu_R
  !!    Returns reference viscosity [Eq. (3.6)]
  !!
  !! eosPR_mu(this,D1,D2,U1,U2,delta_w0,Re0,mu_R)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  D1,D2 [real, in] - reference densities [kg/m3]
  !!  U1,U2 [real, in] - reference velocities [m/s]
  !!  delta_w0 [real, in] - initial vorticity tickness [m]
  !!  Re0 [real, in] - reference Re numbe
  !!  mu_R [real, out] - reference viscosity [kg/m.s] 
  !===========================================================================  
  subroutine eosPR_mu_R(this,D1,D2,U1,U2,delta_w0,Re0,mu_R)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: D1,D2,U1,U2,delta_w0,Re0
    real(PREC), intent(out) :: mu_R
    real(PREC) :: Delta_U0
    !
    Delta_U0=abs(U1-U2)
    mu_R=0.5_PREC*(D1+D2)*Delta_U0*delta_w0/Re0

  end subroutine eosPR_mu_R

  !!============================================================================
  !!
  !! Sub: eosPR_Sc
  !!    Returns Schmidt number [Eq. (2.28)]
  !!
  !! eosPR_Sc(this,massFrac,Sc)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  Sc [real, out] - Schmidt number
  !===========================================================================  
  subroutine eosPR_Sc(this,massFrac,Sc)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(out) :: Sc
    !
    Sc=1.5_PREC-massFrac(2)

  end subroutine eosPR_Sc

  !!============================================================================
  !!
  !! Sub: eosPR_massD
  !!    Returns binary diffusion coefficient [Eq. (2.28)]
  !!
  !! eosPR_massD(this,mu,T,D,massFrac,massD)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  mu [real, in] - viscosity [kg/m.s] 
  !!  T [real, in] - Temperature [K]
  !!  D [real, in] - Density [kg/m3]  
  !!  massD [real, out] - binary diffusion coefficient [m2/s] 
  !===========================================================================  
  subroutine eosPR_massD(this,mu,T,D,massFrac,massD)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: mu
    real(PREC), intent(in) :: T,D
    real(PREC), intent(out) :: massD
    real(PREC) :: Sc,alpha_D
    !
    call eosPR_Sc(this,massFrac,Sc)
    call eosPR_alpha_D(this,T,D,massFrac,alpha_D)
    massD=mu/(Sc*D*alpha_D)

  end subroutine eosPR_massD
 
  !!============================================================================
  !!
  !! Sub: eosPR_Pr
  !!    Returns Prandtl number [Eq. (2.28)]
  !!
  !! eosPR_Pr(this,massFrac,Pr)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  Pr [real, out] - Prandtl number
  !===========================================================================  
  subroutine eosPR_Pr(this,massFrac,Pr)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(out) :: Pr
    real(PREC) :: Sc
    !
    call eosPR_Sc(this,massFrac,Sc)
    Pr=Sc/(2._PREC*exp(-3._PREC*massFrac(2)/2._PREC))

  end subroutine eosPR_Pr

  !!============================================================================
  !!
  !! Sub: eosPR_lambda
  !!    Returns lambda used in calculation of thermal conductivity [Eq. (2.28)]
  !!
  !! eosPR_lambda(this,mu,T,D,massFrac,lambda)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  mu [real, in] - viscosity [kg/m.s] 
  !!  T [real, in] - Temperature [K]
  !!  D [real, in] - Density [kg/m3]  
  !!  lambda [real, out] - lambda used in calculation of thermal conductivity [W/m.K]
  !===========================================================================  
  subroutine eosPR_lambda(this,mu,T,D,massFrac,lambda)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: mu
    real(PREC), intent(in) :: T,D
    real(PREC), intent(out) :: lambda
    real(PREC) :: Cp,Pr,mixMolecularMass
    !
    call eosPR_Cp(this,T,D,massFrac,Cp) 
    call eosPR_Pr(this,massFrac,Pr)
    call eosPR_mixMolecularMass(this,massFrac,mixMolecularMass)
    lambda=mu*Cp/(Pr*mixMolecularMass)

  end subroutine eosPR_lambda

  !!============================================================================
  !!
  !! Sub: eosPR_lambda_p
  !!    Returns thermal conductivity [Eq. (2.20)]
  !!
  !! eosPR_lambda_p(this,alpha_BK,mu,T,D,massFrac,lambda_p)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  alpha_BK [real, in] - Bearman-Kirkwood thermal diffusion factor
  !!  mu [real, in] - viscosity [kg/m.s] 
  !!  T [real, in] - Temperature [K]
  !!  D [real, in] - Density [kg/m3]  
  !!  lambda_p [real, out] - thermal conductivity [W/m.K] 
  !===========================================================================  
  subroutine eosPR_lambda_p(this,alpha_BK,mu,T,D,massFrac,lambda_p)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: alpha_BK
    real(PREC), intent(in) :: mu
    real(PREC), intent(in) :: T,D
    real(PREC), intent(out) :: lambda_p
    real(PREC) :: mixMolecularMass,lambda,alpha_IK,massD
    real(PREC) :: moleFrac(this%nspecies_)
    ! 
    call eosPR_mixMolecularMass(this,massFrac,mixMolecularMass) 
    call eosPR_moleFrac(this,massFrac,moleFrac)
    call eosPR_lambda(this,mu,T,D,massFrac,lambda)
    call eosPR_alpha_IK(this,T,D,massFrac,alpha_BK,alpha_IK)
    call eosPR_massD(this,mu,T,D,massFrac,massD)
   
    lambda_p=lambda+moleFrac(1)*moleFrac(2)*alpha_IK*alpha_BK*&
                    R_UNIVERSAL*D*massD/mixMolecularMass

  end subroutine eosPR_lambda_p
  !!============================================================================
  !!
  !! Sub: eosPR_alpha_IK
  !!    Returns Irwing-Kirkwood thermal diffusion factor [Eq. (2.21)]
  !!
  !! eosPR_alpha_IK(this,T,D,massFrac,alpha_BK,alpha_IK)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  T [real, in] - Temperature [K]
  !!  D [real, in] - Density [kg/m3]  
  !!  alpha_BK [real, in] - Bearman-Kirkwood thermal diffusion factor
  !!  alpha_IK [real, out] - Irwing-Kirkwood thermal diffusion factor
  !===========================================================================  
  subroutine eosPR_alpha_IK(this,T,D,massFrac,alpha_BK,alpha_IK)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_) 
    real(PREC), intent(in) :: T,D
    real(PREC), intent(in) :: alpha_BK
    real(PREC), intent(out) :: alpha_IK
    real(PREC) :: par_h(this%nspecies_)
    real(PREC) :: mixMolecularMass
    !
    call eosPR_par_h(this,T,D,massFrac,par_h)
    call eosPR_mixMolecularMass(this,massFrac,mixMolecularMass)
    alpha_IK=alpha_BK+(this%molecularMass_(1)*this%molecularMass_(2)/mixMolecularMass)*&
             (par_h(2)/this%molecularMass_(2)-par_h(1)/this%molecularMass_(1))/&
             (R_UNIVERSAL*T)

  end subroutine eosPR_alpha_IK

  !!============================================================================
  !!
  !! Sub: eosPR_alpha_BK
  !!    Returns Bearman-Kirkwood thermal diffusion factor [Eq. (2.21)]
  !!
  !! eosPR_alpha_BK(this,T,D,massFrac,alpha_IK,alpha_BK)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  T [real, in] - Temperature [K]
  !!  D [real, in] - Density [kg/m3]  
  !!  alpha_IK [real, in]  - Irwing-Kirkwood thermal diffusion factor
  !!  alpha_BK [real, out] - Bearman-Kirkwood thermal diffusion factor
  !===========================================================================  
  subroutine eosPR_alpha_BK(this,T,D,massFrac,alpha_IK,alpha_BK)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_) 
    real(PREC), intent(in) :: T,D
    real(PREC), intent(in) :: alpha_IK    
    real(PREC), intent(out) :: alpha_BK
    real(PREC) :: par_h(this%nspecies_)
    real(PREC) :: mixMolecularMass
    !
    call eosPR_par_h(this,T,D,massFrac,par_h)
    call eosPR_mixMolecularMass(this,massFrac,mixMolecularMass)
    alpha_BK=alpha_IK-(this%molecularMass_(1)*this%molecularMass_(2)/mixMolecularMass)*&
             (par_h(2)/this%molecularMass_(2)-par_h(1)/this%molecularMass_(1))/&
             (R_UNIVERSAL*T)

  end subroutine eosPR_alpha_BK


  !!============================================================================
  !!
  !! Sub: eosPR_massFlux_p
  !!    Returns massFlux_p used in calculation of mass flux [Eq. (2.19)]
  !!
  !! eosPR_massFlux_p(this,parYh_parXj,parp_parXj,
  !!                 mu,T,D,massFrac,massFlux_p)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  parYh_parxj [real, in] - hepthane mass fraction gradient in j direction [1/m]
  !!  parp_parxj [real, in] - pressure gradient in j direction [Pa/m]
  !!  mu [real, in] - viscosity [kg/m.s] 
  !!  T [real, in] - Temperature [K]
  !!  D [real, in] - Density [kg/m3]  
  !!  massFlux_p [real, out] - massFlux_p used in calculation of mass flux [kg/m2.s]
  !===========================================================================  
  subroutine eosPR_massFlux_p(this,parYh_parxj,parp_parxj,&
                             mu,T,D,massFrac,massFlux_p)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: parYh_parxj,parp_parxj
    real(PREC), intent(in) :: mu
    real(PREC), intent(in) :: T,D
    real(PREC), intent(out) :: massFlux_p
    real(PREC) :: mixMolecularMass
    real(PREC) :: alpha_D,massD
    real(PREC) :: par_molarVolume(this%nspecies_)
    !
    call eosPR_mixMolecularMass(this,massFrac,mixMolecularMass)
    call eosPR_alpha_D(this,T,D,massFrac,alpha_D)
    call eosPR_massD(this,mu,T,D,massFrac,massD)
    call eosPR_par_molarVolume(this,T,D,massFrac,par_molarVolume)
    massFlux_p=D*massD*(alpha_D*parYh_parxj+massFrac(1)*massFrac(2)/(R_UNIVERSAL*T)*&
          this%molecularMass_(1)*this%molecularMass_(2)/mixMolecularMass*&
          (par_molarVolume(2)/this%molecularMass_(2)-&
          par_molarVolume(1)/this%molecularMass_(1))*&
          parp_parxj)

  end subroutine eosPR_massFlux_p

  !!============================================================================
  !!
  !! Sub: eosPR_massFlux_BK
  !!    Returns mass flux [Eq. (2.18)] given Bearman-Kirkwood thermal diffusion factor
  !!
  !! eosPR_massFlux_BK(this,parT_parxj,parYh_parXj,parp_parXj,
  !!        alpha_BK,mu,T,D,massFrac,massFlux)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  parT_parxj [real, in] - temperature gradient in j direction [T/m]
  !!  parYh_parxj [real, in] - hepthane mass fraction gradient in j direction [1/m]
  !!  parp_parxj [real, in] - pressure gradient in j direction [Pa/m]
  !!  alpha_BK [real, in] - Bearman-Kirkwood thermal diffusion factor
  !!  mu [real, in] - viscosity [kg/m.s] 
  !!  T [real, in] - Temperature [K]
  !!  D [real, in] - Density [kg/m3]  
  !!  massFlux [real, out] - mass flux [kg/m2.s]
  !===========================================================================  
  subroutine eosPR_massFlux_BK(this,parT_parxj,parYh_parxj,parp_parxj,&
                    alpha_BK,mu,T,D,massFrac,massFlux)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: parT_parxj,parYh_parxj,parp_parxj
    real(PREC), intent(in) :: alpha_BK
    real(PREC), intent(in) :: mu
    real(PREC), intent(in) :: T,D
    real(PREC), intent(out) :: massFlux
    real(PREC) :: massD
    real(PREC) :: massFlux_p
    !
    call eosPR_massD(this,mu,T,D,massFrac,massD)
    call eosPR_massFlux_p(this,parYh_parxj,parp_parxj,mu,T,D,massFrac,massFlux_p)
    massFlux=-(massFlux_p+alpha_BK*massFrac(1)*massFrac(2)*D*massD*parT_parxj/T)

  end subroutine eosPR_massFlux_BK

  !!============================================================================
  !!
  !! Sub: eosPR_massFlux_IK
  !!    Returns mass flux [Eq. (2.18)] given Irwing-Kirkwood thermal diffusion factor
  !!
  !! eosPR_massFlux_IK(this,parT_parxj,parYh_parXj,parp_parXj,
  !!        alpha_IK,mu,T,D,massFrac,massFlux)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  parT_parxj [real, in] - temperature gradient in j direction [T/m]
  !!  parYh_parxj [real, in] - hepthane mass fraction gradient in j direction [1/m]
  !!  parp_parxj [real, in] - pressure gradient in j direction [Pa/m]
  !!  alpha_IK [real, in] - Irwing-Kirkwood thermal diffusion factor
  !!  mu [real, in] - viscosity [kg/m.s] 
  !!  T [real, in] - Temperature [K]
  !!  D [real, in] - Density [kg/m3]  
  !!  massFlux [real, out] - mass flux [kg/m2.s]
  !===========================================================================  
  subroutine eosPR_massFlux_IK(this,parT_parxj,parYh_parxj,parp_parxj,&
                    alpha_IK,mu,T,D,massFrac,massFlux)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: parT_parxj,parYh_parxj,parp_parxj
    real(PREC), intent(in) :: alpha_IK
    real(PREC), intent(in) :: mu
    real(PREC), intent(in) :: T,D
    real(PREC), intent(out) :: massFlux
    real(PREC) :: alpha_BK    
    real(PREC) :: massD
    real(PREC) :: massFlux_p    
    !
    call eosPR_massD(this,mu,T,D,massFrac,massD)
    call eosPR_alpha_BK(this,T,D,massFrac,alpha_IK,alpha_BK)       
    call eosPR_massFlux_p(this,parYh_parxj,parp_parxj,mu,T,D,massFrac,massFlux_p)
    massFlux=-(massFlux_p+alpha_BK*massFrac(1)*massFrac(2)*D*massD*parT_parxj/T)
  end subroutine eosPR_massFlux_IK


  !!============================================================================
  !!
  !! Sub: eosPR_heatFlux_BK
  !!    Returns heat flux [Eq. (2.17)] given Bearman-Kirkwood thermal diffusion factor
  !!
  !! eosPR_heatFlux_BK(this,parT_parxj,parYh_parXj,parp_parXj,
  !!        alpha_BK,mu,T,D,massFrac,heatFlux)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  parT_parxj [real, in] - temperature gradient in j direction [T/m]
  !!  parYh_parxj [real, in] - hepthane mass fraction gradient in j direction [1/m]
  !!  parp_parxj [real, in] - pressure gradient in j direction [Pa/m]
  !!  alpha_BK [real, in] - Bearman-Kirkwood thermal diffusion factor
  !!  mu [real, in] - viscosity [kg/m.s] 
  !!  T [real, in] - Temperature [K]
  !!  D [real, in] - Density [kg/m3]  
  !!  heatFlux [real, out] - heat flux [J/m2.s]
  !===========================================================================  
  subroutine eosPR_heatFlux_BK(this,parT_parxj,parYh_parxj,parp_parxj,&
                    alpha_BK,mu,T,D,massFrac,heatFlux)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: parT_parxj,parYh_parxj,parp_parxj
    real(PREC), intent(in) :: alpha_BK
    real(PREC), intent(in) :: mu
    real(PREC), intent(in) :: T,D
    real(PREC), intent(out) :: heatFlux
    real(PREC) :: mixMolecularMass
    real(PREC) :: massFlux_p
    real(PREC) :: alpha_IK
    real(PREC) :: lambda_p
    !
    call eosPR_mixMolecularMass(this,massFrac,mixMolecularMass)    
    call eosPR_massFlux_p(this,parYh_parxj,parp_parxj,mu,T,D,massFrac,massFlux_p)
    call eosPR_alpha_IK(this,T,D,massFrac,alpha_BK,alpha_IK)   
    call eosPR_lambda_p(this,alpha_BK,mu,T,D,massFrac,lambda_p)
    heatFlux=-(lambda_p*parT_parxj+alpha_IK*R_UNIVERSAL*T*&
              (mixMolecularMass/(this%molecularMass_(1)*this%molecularMass_(2)))*massFlux_p)

  end subroutine eosPR_heatFlux_BK 
  
  !!============================================================================
  !!
  !! Sub: eosPR_heatFlux_IK
  !!    Returns heat flux [Eq. (2.17)] given Irwing-Kirkwood thermal diffusion factor
  !!
  !! eosPR_heatFlux_IK(this,parT_parxj,parYh_parXj,parp_parXj,
  !!        alpha_IK,mu,T,D,massFrac,heatFlux)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  parT_parxj [real, in] - temperature gradient in j direction [T/m]
  !!  parYh_parxj [real, in] - hepthane mass fraction gradient in j direction [1/m]
  !!  parp_parxj [real, in] - pressure gradient in j direction [Pa/m]
  !!  alpha_IK [real, in] - Irwing-Kirkwood thermal diffusion factor
  !!  mu [real, in] - viscosity [kg/m.s] 
  !!  T [real, in] - Temperature [K]
  !!  D [real, in] - Density [kg/m3]  
  !!  heatFlux [real, out] - heat flux [J/m2.s]
  !===========================================================================  
  subroutine eosPR_heatFlux_IK(this,parT_parxj,parYh_parxj,parp_parxj,&
                    alpha_IK,mu,T,D,massFrac,heatFlux)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: parT_parxj,parYh_parxj,parp_parxj
    real(PREC), intent(in) :: alpha_IK
    real(PREC), intent(in) :: mu
    real(PREC), intent(in) :: T,D
    real(PREC), intent(out) :: heatFlux
    real(PREC) :: mixMolecularMass
    real(PREC) :: alpha_BK    
    real(PREC) :: massFlux_p
    real(PREC) :: lambda_p
    !
    call eosPR_mixMolecularMass(this,massFrac,mixMolecularMass)    
    call eosPR_massFlux_p(this,parYh_parxj,parp_parxj,mu,T,D,massFrac,massFlux_p)
    call eosPR_alpha_BK(this,T,D,massFrac,alpha_IK,alpha_BK)
    call eosPR_lambda_p(this,alpha_BK,mu,T,D,massFrac,lambda_p)
    heatFlux=-(lambda_p*parT_parxj+alpha_IK*R_UNIVERSAL*T*&
              (mixMolecularMass/(this%molecularMass_(1)*this%molecularMass_(2)))*massFlux_p)

  end subroutine eosPR_heatFlux_IK                                                                             

  !!============================================================================
  !!
  !! Sub: eosPR_molarVolume2
  !!    Returns molar specific volume partial molar volumes [Page 7]
  !!
  !! eosPR_molarVolume2(this,T,D,massFrac,molarVolume2)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  D [real, in] - density [kg/m3]
  !!  T [real, in] - temperature [K]
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  molarVolume2 [real, out] - molar specific volume [m3/kmol]  
  !=========================================================================== 
  subroutine eosPR_molarVolume2(this,T,D,massFrac,molarVolume2)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: T,D
    real(PREC), intent(out) :: molarVolume2
    real(PREC) :: par_molarVolume(this%nspecies_)
    real(PREC) :: moleFrac(this%nspecies_)
    integer :: i
    !
    call eosPR_par_molarVolume(this,T,D,massFrac,par_molarVolume)
    call eosPR_moleFrac(this,massFrac,moleFrac)
    molarVolume2=0._PREC
    do i=1,this%nspecies_
       molarVolume2=molarVolume2+moleFrac(i)*par_molarVolume(i)
    enddo
  end subroutine eosPR_molarVolume2

  !!============================================================================
  !!
  !! Sub: eosPR_enthalpy2
  !!    Returns molar specific enthalpy using partial molar specific enthalpies
  !!    [Page 7]
  !!
  !! eosPR_enthalpy2(this,T,D,massFrac,h)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  D [real, in] - density kg/m3
  !!  T [real, in] - temperature K
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  h [real, out] - molar specific enthalpy J/kmol  
  !=========================================================================== 
  subroutine eosPR_enthalpy2(this,T,D,massFrac,h)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: T,D
    real(PREC), intent(out) :: h
    real(PREC) :: par_h(this%nspecies_)
    real(PREC) :: moleFrac(this%nspecies_)
    integer :: i
    !
    call eosPR_par_h(this,T,D,massFrac,par_h)
    call eosPR_moleFrac(this,massFrac,moleFrac)
    h=0._PREC
    do i=1,this%nspecies_
       h=h+moleFrac(i)*par_h(i)
    enddo
  end subroutine eosPR_enthalpy2

  !!============================================================================
  !!
  !! Sub: eosPR_energy
  !!    Returns molar internal energy [e=h-pv]
  !!
  !! eosPR_energy(this,T,D,massFrac,e)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  T [real, in] - Temperature [K]
  !!  D [real, in] - Density [kg/m3] 
  !!  e [real, out] - Molar internal enthalpy [J/kmol]    
  !===========================================================================  
  subroutine eosPR_energy(this,T,D,massFrac,e)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: T,D
    real(PREC), intent(out) :: e
    real(PREC) :: v,p,h
    !
    call eosPR_molarVolume(this,D,massFrac,v)
    call eosPR_pressure(this,T,D,massFrac,p)
    call eosPR_enthalpy(this,T,D,massFrac,h)
    e=h-p*v

  end subroutine eosPR_energy

  !!============================================================================
  !!
  !! Sub: eosPR_temperature_from_e
  !!    Returns temperature using internal energy [Eq. (2.31)]
  !!
  !! eosPR_temperature_from_e(this,internal_e,D,massFrac,T)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  D [real, in] - density kg/m3
  !!  internal_e [real, in] - internal energy [J/kmol]
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  T [real, out] - temperature [K] 
  !=========================================================================== 
  subroutine eosPR_temperature_from_e(this,internal_e,D,massFrac,T)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: internal_e,D
    real(PREC), intent(out) :: T
    real(PREC) :: E(4,4) 
    real(PREC) :: ek(4)
    real(PREC) :: teta,e_L,e_U
    real(PREC) :: mixMolecularMass,internal_eMJ
    integer :: i
    !
    call eosPR_mixMolecularMass(this,massFrac,mixMolecularMass)
    internal_eMJ=internal_e/mixMolecularMass/1000000._PREC

    E(1,1)=0.44142_PREC
    E(1,2)=0.51813_PREC
    E(1,3)=-0.0018118_PREC
    E(1,4)=0.00000612810_PREC

    E(2,1)=0.41489_PREC 
    E(2,2)=0.47173_PREC
    E(2,3)=-0.044262_PREC
    E(2,4)=0.00225095_PREC

    E(3,1)=0.86090_PREC
    E(3,2)=1.78435_PREC
    E(3,3)=-0.0019375_PREC
    E(3,4)=0.00000478990_PREC

    E(4,1)=0.85262_PREC
    E(4,2)=1.76317_PREC
    E(4,3)=-0.049197_PREC
    E(4,4)=0.00259000_PREC 

    do i=1,4
       ek(i)=E(i,1)+E(i,2)*massFrac(2)+E(i,3)*massFrac(2)**2+E(i,4)*massFrac(2)**3
    enddo

    teta=1.14_PREC+0.667_PREC*massFrac(2)**0.676_PREC
    
    e_L=ek(1)+(D-this%D_L_)/(this%D_U_-this%D_L_)*(ek(2)-ek(1))
    e_U=ek(3)+(D-this%D_L_)/(this%D_U_-this%D_L_)*(ek(4)-ek(3))
    
   T=(this%T_L_**teta+(internal_eMJ-e_L)/(e_U-e_L)*(this%T_U_**teta-&
     this%T_L_**teta))**(1._PREC/teta)
   
  end subroutine eosPR_temperature_from_e

  !!============================================================================
  !!
  !! Sub: eosPR_temperature_from_EOS
  !!    Returns temperature using equation of state [Eq. (2.6)]
  !!
  !! eosPR_temperature_from_EOS(this,p,D,massFrac,T)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  D [real, in] - density kg/m3
  !!  p [real, in] - pressure [Pa]
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  T [real, out] - temperature [K] 
  !===========================================================================
  subroutine eosPR_temperature_from_EOS(this,p,D,massFrac,T)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: p,D
    real(PREC), intent(out) :: T
    real(PREC) :: molarVolume
    real(PREC) :: moleFrac(this%nspecies_)
    real(PREC) :: Bm,Am
    real(PREC) :: A(this%nspecies_,this%nspecies_)
    real(PREC) :: B(this%nspecies_)
    real(PREC) :: C(this%nspecies_,this%nspecies_)
    real(PREC) :: Tc(this%nspecies_,this%nspecies_)
    real(PREC) :: vc(this%nspecies_,this%nspecies_)
    real(PREC) :: pc(this%nspecies_,this%nspecies_)
    real(PREC) :: G,a2,a1,a0
    real(PREC) :: Q(this%nspecies_,this%nspecies_)
    real(PREC) :: x(2)
    integer :: i,j
    !
    call eosPR_molarVolume(this,D,massFrac,molarVolume)
    call eosPR_moleFrac(this,massFrac,moleFrac)
    call eosPR_Bm(this,T,moleFrac,Bm)
    call eosPR_ABC(this,T,moleFrac,A,B,C,Tc,vc,pc)
   
    G=R_UNIVERSAL/(molarVolume-Bm)
    do i=1,this%nspecies_
       do j=1,this%nspecies_
          Q(i,j)=moleFrac(i)*moleFrac(j)*0.457236_PREC*(R_UNIVERSAL*Tc(i,j))**2/&
                              (pc(i,j)*(molarVolume**2+2*molarVolume*Bm-Bm**2))
       enddo
    enddo
    
    a2=G
    do i=1,this%nspecies_
       do j=1,this%nspecies_
          a2=a2-Q(i,j)*C(i,j)**2/Tc(i,j)
       enddo
    enddo
    a1=0._PREC
    do i=1,this%nspecies_
       do j=1,this%nspecies_
          a1=a1+2._PREC*Q(i,j)*C(i,j)*(1._PREC+C(i,j))/sqrt(Tc(i,j))
       enddo
    enddo
    a0=0._PREC
    do i=1,this%nspecies_
       do j=1,this%nspecies_
          a0=a0-Q(i,j)*(1._PREC+C(i,j))**2
       enddo
    enddo

   x(1)=(-a1+sqrt(a1**2-4._PREC*a2*(a0-p)))/(2._PREC*a2)
   x(2)=(-a1-sqrt(a1**2-4._PREC*a2*(a0-p)))/(2._PREC*a2)
    
   do i=1,2
      if (x(i)>0) then
         if (x(i)**2<4000.) then
             T=x(i)**2 
         endif
      endif
   enddo

  end subroutine eosPR_temperature_from_EOS

  !!============================================================================
  !!
  !! Sub: eosPR_h0_species
  !!    Returns h0 - Reference molar specific enthalpies [Eq. (2.30)]
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  T [real, in] - temperature [K]
  !!  h0 [real(1:nspecies), out] - [J/kmol] 
  !===========================================================================  
  subroutine eosPR_h0_species(this,T,h0)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: T
    real(PREC), intent(out) :: h0(this%nspecies_)

    h0(1)=this%molecularMass_(1)*656.72_PREC*T**(1.071_PREC)
    h0(2)=this%molecularMass_(2)*27.877_PREC*T**(1.6414_PREC)

  end subroutine eosPR_h0_species

  !!============================================================================
  !!
  !! Sub: eosPR_density_from_EOS
  !!    Returns density using pressure [Eq. (2.6)]
  !!
  !! eosPR_density_from_EOS(this,T,p,massFrac,D)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  T [real, in] - Temperature [K]
  !!  D [real, out] - Density [kg/m3]  
  !!  p [real, in] - Pressure [Pa]    
  !===========================================================================  
  subroutine eosPR_density_from_EOS(this,T,p,massFrac,D)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: T,p
    real(PREC), intent(out) :: D
    real(PREC) :: molarVolume
    real(PREC) :: moleFrac(this%nspecies_)
    real(PREC) :: Bm,Am
    real(PREC) :: A(this%nspecies_,this%nspecies_)
    real(PREC) :: B(this%nspecies_)
    real(PREC) :: a3,a2,a1,a0
    complex(PREC) ::x(3)
    real(PREC) :: aux
    real(PREC) :: D_L,D_U
    integer :: i,n_rp
    !
    D_L=0.01_PREC  ![kg/m3]
    D_U=1000._PREC ![kg/m3]
    call eosPR_moleFrac(this,massFrac,moleFrac)
    call eosPR_AmBm(this,T,moleFrac,Am,Bm,A,B)
    a3=p
    a2=Bm*p-R_UNIVERSAL*T
    a1=Am-2*R_UNIVERSAL*T*Bm-3*Bm**2*p
    a0=p*Bm**3+(R_UNIVERSAL*T*Bm-Am)*Bm

    call cubic_root(a3,a2,a1,a0,x)
 
    D=0._PREC
    n_RP=0
    
    !x(1) = 0.51927
    !x(2) = complex(0.10815,0.23736)
    !x(3) = complex(0.10815,-0.23736)
    
    !print*,a3,a2,a1,a0
    do i=1,3
    !print*,x(i)
    !print*,Real(x(i)), AIMAG(x(i))
    !print*,Real(Real(x(i))),Real(AIMAG(x(i)))
    !print*,abs(Real(AIMAG(x(i))))

    if ((abs(Real(AIMAG(x(i))))<10.**(-14)) .AND. (Real(Real(x(i)))>0)) then
       print*,i
       n_rp=n_rp+1
       molarVolume=Real(Real(x(i)),PREC)
       call eosPR_density(this,molarVolume,massFrac,aux)  
       if ((aux.ge.D_L) .AND. (aux.le.D_U)) then
       !if ((aux.ge.this%D_L_) .AND. (aux.le.this%D_U_)) then
          D=aux
       endif
    endif
    enddo

    if (D==0) then
       print*,"Error: density is out of the range"
    endif

  end subroutine eosPR_density_from_EOS

  !!============================================================================
  !!
  !! Sub: cubic_root
  !!    Returns roots of cubic function [https://en.wikipedia.org/wiki/Cubic_function]
  !!
  !! cubic_root(a,b,c,d,aux)
  !!
  !! a,b,c,d [real, in] 
  !! aux [complex, out] 
  !!============================================================================
  subroutine cubic_root(a,b,c,d,aux)
  implicit none
  real(PREC), intent(in) :: a,b,c,d
  complex(PREC),intent(out) :: aux(3)
  real(PREC) :: x(3)
  complex(PREC) :: i
  real(PREC) :: Delta
  complex(PREC) :: Delta_0,Delta_1,C1,C2,C3
  
  Delta=18*a*b*c*d-4*b**3*d+b**2*c**2-4*a*c**3-27*a**2*d**2
  if (Delta>0) then 
     !print*,'the equation has three distinct real roots.'
  elseif (Delta==0.) then 
     !print*,'the equation has a multiple root and all of its roots are real.'
  else
     !print*,'the equation has one real root and two non-real complex conjugate roots.'
  endif
  Delta_0=dcmplx(b**2-3*a*c,0._PREC)
  Delta_1=dcmplx(2*b**3-9*a*b*c+27*a**2*d,0._PREC)
  !print*,'a',a
  !print*,'b',b
  !print*,'c',c
  !print*,'d',d
  !print*,'Delta',Delta
  !print*,'Delta_0',Delta_0
  !print*,'Delta_1',Delta_1
  if ((Delta==0.) .AND. (Delta_0==0.)) then
     !print*,'1'
     aux(1)=-b/(3._PREC*a)
     aux(2)=-b/(3._PREC*a)
     aux(3)=-b/(3._PREC*a)
  elseif ((Delta==0.) .AND. (Delta_0/=0.)) then
     !print*,'2'
     aux(1)=(9._PREC*a*d-b*c)/(2*Delta_0)
     aux(2)=(9._PREC*a*d-b*c)/(2*Delta_0)
     aux(3)=(4._PREC*a*b*c-9._PREC*a**2*d-b**3)/(a*Delta_0)
  else
     C1=(0.5_PREC*(Delta_1+sqrt(Delta_1**2-4*Delta_0**3)))**(1._PREC/3._PREC)
     if (C1==0.) then
        C1=(0.5_PREC*(Delta_1-sqrt(Delta_1**2-4*Delta_0**3)))**(1._PREC/3._PREC)
     endif
     i=(0._PREC,1._PREC)
     aux(1)=-1/(3._PREC*a)*(b+C1+Delta_0/C1)
     C2=(-0.5_PREC+0.5_PREC*sqrt(3._PREC)*i)*C1
     C3=(-0.5_PREC-0.5_PREC*sqrt(3._PREC)*i)*C1
     aux(2)=-1/(3._PREC*a)*(b+C2+Delta_0/C2)
     aux(3)=-1/(3._PREC*a)*(b+C3+Delta_0/C3)
  endif

end subroutine cubic_root

  !!============================================================================
  !!
  !! Sub: eosPR_IntCp
  !!    Returns change in enthalpy by integrating Cp 
  !!    Integration based on Simpson's rule 
  !!    [https://en.wikipedia.org/wiki/Simpson%27s_rule]
  !!
  !! eosPR_IntCp(this,T1,T2,D,massFrac,Delta_h)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  T1 and T2 [real, in] - Initial and final temperatures [K]
  !!  D [real, in] - Density [kg/m3] 
  !!  Delta_h [real, out] - Molar enthalpy change [J/kmol]    
  !===========================================================================  
  subroutine eosPR_IntCp(this,T1,T2,D,massFrac,Delta_h)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: T1,T2,D
    real(PREC), intent(out) :: Delta_h
    real(PREC) :: T_a,T_b
    real(PREC) :: Cp_a,Cp_b,Cp_ab2
    integer :: n,i
    !
    n=int(T2-T1)*100 
    Delta_h=0._PREC
    do i=1,n
       T_a=T1+(i-1)*((T2-T1)/n)
       T_b=T1+(i)*((T2-T1)/n)
       call eosPR_Cp(this,T_a,D,massFrac,Cp_a)
       call eosPR_Cp(this,T_b,D,massFrac,Cp_b)
       call eosPR_Cp(this,((T_a+T_b)/2._PREC),D,massFrac,Cp_ab2)
       !--- TEST ---!
       !Cp_a=T_a-1._PREC/T_a+T_a**2
       !Cp_b=T_b-1._PREC/T_b+T_b**2
       !Cp_ab2=(T_a+T_b)/2-2._PREC/(T_a+T_b)+ &
       !     ((T_a+T_b)/2._PREC)**2
       !-------------
       Delta_h=(T_b-T_a)/6._PREC*(Cp_a+4._PREC*Cp_ab2+Cp_b)+Delta_h
     enddo

  end subroutine eosPR_IntCp
 
  !!============================================================================
  !!
  !! Sub: eosPR_enthalpy
  !!    Returns the molar specific enthalpy [Eq. (2.9)]
  !!
  !! eosPR_enthalpy(this,T,D,massFrac,h)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  T [real, in] - Temperature [K]
  !!  D [real, in] - Density [kg/m3] 
  !!  h [real, out] - Molar specific enthalpy [J/kmol]    
  !===========================================================================  
  subroutine eosPR_enthalpy(this,T,D,massFrac,h)
    implicit none
    type(eosPR), intent(inout) :: this 
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: T,D
    real(PREC), intent(out) :: h
    real(PREC) :: molarVolume,p
    real(PREC) :: moleFrac(this%nspecies_)
    real(PREC) :: Am,K1,dAmdT,h0

    !
    call eosPR_molarVolume(this,D,massFrac,molarVolume)
    call eosPR_moleFrac(this,massFrac,moleFrac)
    call eosPR_Am(this,T,moleFrac,Am)
    call eosPR_pressure(this,T,D,massFrac,p)
    call eosPR_K1(this,T,D,massFrac,K1)
    call eosPR_dAmdT(this,T,moleFrac,dAmdT)
    call eosPR_h0(this,T,massFrac,h0)
    h=h0+p*molarVolume-R_UNIVERSAL*T+K1*(Am-T*dAmdT) ![J/kmol]

  end subroutine eosPR_enthalpy
  
  !!============================================================================
  !!
  !! Sub: eosPR_Cp
  !!    Returns the molar specific heat at constant pressure (Cp) [Eq. (2.10)]
  !!
  !! eosPR_Cp(this,T,D,massFrac,Cp)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  T [real, in] - Temperature [K]
  !!  D [real, in] - Density [kg/m3] 
  !!  Cp [real, out] - Molar specific heat at constant pressure [J/kmol.K]    
  !===========================================================================  
  subroutine eosPR_Cp(this,T,D,massFrac,Cp)
    implicit none
    type(eosPR), intent(inout) :: this 
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: T,D
    real(PREC), intent(out) :: Cp
    real(PREC) :: molarVolume
    real(PREC) :: moleFrac(this%nspecies_)
    real(PREC) :: Am,K1,d2AmdT2,dpdT,dpdv,Cp0
    !
    call eosPR_molarVolume(this,D,massFrac,molarVolume)
    call eosPR_moleFrac(this,massFrac,moleFrac)
    call eosPR_Am(this,T,moleFrac,Am)
    call eosPR_dpdT(this,T,D,massFrac,dpdT)
    call eosPR_dpdv(this,T,D,massFrac,dpdv)
    call eosPR_K1(this,T,D,massFrac,K1)
    call eosPR_d2AmdT2(this,T,massFrac,d2AmdT2)
    call eosPR_Cp0(this,T,massFrac,Cp0)
    call eosPR_parp_parv_TX(this,T,D,massFrac,dpdv)
    Cp=Cp0-T*(dpdT**2)/(dpdv)-R_UNIVERSAL-T*d2AmdT2*K1 ![J/kmol]
    
  end subroutine eosPR_Cp
  
  !!============================================================================
  !!
  !! Sub: eosPR_spSound
  !!    Returns the speed of the sound [Eq. (2.14)]
  !!
  !! eosPR_spSound(this,T,D,massFrac,as)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  T [real, in] - Temperature [K]
  !!  D [real, in] - Density [kg/m3] 
  !!  as [real, out] - The speed of sound [m/s]    
  !===========================================================================  
  subroutine eosPR_spSound(this,T,D,massFrac,as)
    implicit none
    type(eosPR), intent(inout) :: this 
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: T,D
    real(PREC), intent(out) :: as
    real(PREC) :: molarVolume,dpdT,dpdv,Cp
    real(PREC) :: ks,alpha,kT
    !
    call eosPR_molarVolume(this,D,massFrac,molarVolume)
    call eosPR_dpdT(this,T,D,massFrac,dpdT)
    call eosPR_dpdv(this,T,D,massFrac,dpdv)
    call eosPR_Cp(this,T,D,massFrac,Cp)
    alpha=-1._PREC*dpdT/(molarVolume*dpdv)
    kT=-1._PREC/(molarVolume*dpdv)
    ks=kT-molarVolume*T*alpha**2/Cp
    as=sqrt(1._PREC/(D*ks)) ![J/kmol]
    
  end subroutine eosPR_spSound

  !!============================================================================
  !!
  !! Sub: eosPR_Para
  !!    Returns A,B,C,G,pc,Tc [Appendix A]
  !!
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  T [real, in] - Temperature [K]
  !!  Am [real, out] - [(m3/kmol)2]
  !!  Bm [real, out] - [m3/kmol] 
  !===========================================================================  
  subroutine eosPR_Para(this,T,A,B,C,G,pc,Tc)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: T
    real(PREC), intent(out) :: A(this%nspecies_,this%nspecies_)
    real(PREC), intent(out) :: B(this%nspecies_)
    real(PREC), intent(out) :: C(this%nspecies_,this%nspecies_)
    real(PREC), intent(out) :: G(this%nspecies_,this%nspecies_)
    real(PREC), intent(out) :: Tc(this%nspecies_,this%nspecies_)
    real(PREC), intent(out) :: pc(this%nspecies_,this%nspecies_)
    real(PREC) :: omega(this%nspecies_,this%nspecies_)
    real(PREC) :: zc(this%nspecies_,this%nspecies_)
    real(PREC) :: k(this%nspecies_,this%nspecies_)
    real(PREC) :: vc(this%nspecies_,this%nspecies_)
    integer :: i,j
    ! 
    omega(1,1)=this%omega_(1)
    omega(2,2)=this%omega_(2)
    omega(1,2)=0.5_PREC*(this%omega_(1)+this%omega_(2)) 
    omega(2,1)=omega(1,2)

    zc(1,1)=this%zc_(1)
    zc(2,2)=this%zc_(2)
    zc(1,2)=0.5_PREC*(this%zc_(1)+this%zc_(2))
    zc(2,1)=zc(1,2)
    
    k(1,1)=0._PREC
    k(2,2)=0._PREC
    k(1,2)=0.1_PREC
    k(2,1)=0.1_PREC

    Tc(1,1)=this%Tc_(1)
    Tc(2,2)=this%Tc_(2)
    Tc(1,2)=sqrt(this%Tc_(1)*this%Tc_(2))*(1._PREC-k(1,2))
    Tc(2,1)=Tc(1,2)

    vc(1,1)=this%vc_(1)
    vc(2,2)=this%vc_(2)
    vc(1,2)=(this%vc_(1)**(1._PREC/3._PREC)+this%vc_(2)**(1._PREC/3._PREC))**3/8._PREC
    vc(2,1)=vc(1,2)

    pc(1,1)=this%pc_(1)
    pc(2,2)=this%pc_(2)
    pc(1,2)=zc(1,2)*R_UNIVERSAL*Tc(1,2)/vc(1,2)
    pc(2,1)=pc(1,2)

    do i=1,this%nspecies_
       do j=1,this%nspecies_
          C(i,j)=0.37464_PREC+1.54226_PREC*omega(i,j)-0.26992_PREC*omega(i,j)**2
          G(i,j)=(C(i,j)*sqrt(T/Tc(i,j)))/(1._PREC+C(i,j)*(1._PREC-sqrt(T/Tc(i,j))))
       enddo
    enddo

    B(1)=0.077796_PREC*R_UNIVERSAL*this%Tc_(1)/this%pc_(1)
    B(2)=0.077796_PREC*R_UNIVERSAL*this%Tc_(2)/this%pc_(2)

    do i=1,this%nspecies_
       do j=1,this%nspecies_
          A(i,j)=0.457236_PREC*(R_UNIVERSAL*Tc(i,j))**2* &
                   ((1._PREC+C(i,j)*(1._PREC-sqrt(T/Tc(i,j))))**2)/pc(i,j)
       enddo
    enddo
 
  end subroutine eosPR_Para
  
  !!============================================================================
  !!
  !! Sub: eosPR_Am
  !!    Returns Am [Eq. (2.7)]
  !!
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  
  !!  Am [real, out] - [(m3/kmol)2]
  !===========================================================================  
  subroutine eosPR_Am(this,T,moleFrac,Am)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: moleFrac(this%nspecies_)
    real(PREC), intent(in) :: T
    real(PREC), intent(out) :: Am
    real(PREC)  :: A(this%nspecies_,this%nspecies_)
    real(PREC)  :: B(this%nspecies_)
    real(PREC)  :: C(this%nspecies_,this%nspecies_)
    real(PREC)  :: G(this%nspecies_,this%nspecies_)
    real(PREC)  :: Tc(this%nspecies_,this%nspecies_)
    real(PREC)  :: pc(this%nspecies_,this%nspecies_)
    
    integer :: i,j
    
    call eosPR_Para(this,T,A,B,C,G,pc,Tc)
 
    Am=0._PREC

    do i=1,this%nspecies_
       do j=1,this%nspecies_
          Am=moleFrac(i)*moleFrac(j)*A(i,j)+Am
       enddo
    enddo 

  end subroutine eosPR_Am
  
  !!============================================================================
  !!
  !! Sub: eosPR_Bm
  !!    Returns Bm [Eq. (2.7)]
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  
  !!  Bm [real, out] - [(m3/kmol)2]
  !===========================================================================  
  subroutine eosPR_Bm(this,T,moleFrac,Bm)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: moleFrac(this%nspecies_)
    real(PREC), intent(in) :: T
    real(PREC), intent(out) :: Bm
    real(PREC)  :: A(this%nspecies_,this%nspecies_)
    real(PREC)  :: B(this%nspecies_)
    real(PREC)  :: C(this%nspecies_,this%nspecies_)
    real(PREC)  :: G(this%nspecies_,this%nspecies_)
    real(PREC)  :: Tc(this%nspecies_,this%nspecies_)
    real(PREC)  :: pc(this%nspecies_,this%nspecies_)
    
    integer :: i,j
    
    call eosPR_Para(this,T,A,B,C,G,pc,Tc)
    
    Bm=0._PREC

     do j=1,this%nspecies_
        Bm=moleFrac(j)*B(j)+Bm
     enddo

  end subroutine eosPR_Bm
  
  !!============================================================================
  !!
  !! Sub: eosPR_dAmdT
  !!    Returns dAmdT [Appendix B]
  !!
  !!  eosPR_dAmdT(this,T,moleFrac,dAmdT)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!
  !!  dAmdT [real, out] 
  !===========================================================================  
  subroutine eosPR_dAmdT(this,T,moleFrac,dAmdT)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: moleFrac(this%nspecies_)
    real(PREC), intent(in) :: T
    real(PREC), intent(out) :: dAmdT
    real(PREC)  :: A(this%nspecies_,this%nspecies_)
    real(PREC)  :: B(this%nspecies_)
    real(PREC)  :: C(this%nspecies_,this%nspecies_)
    real(PREC)  :: G(this%nspecies_,this%nspecies_)
    real(PREC)  :: Tc(this%nspecies_,this%nspecies_)
    real(PREC)  :: pc(this%nspecies_,this%nspecies_)
    
    integer :: i,j
    
    call eosPR_Para(this,T,A,B,C,G,pc,Tc)
    
    dAmdT=0._PREC

    do i=1,this%nspecies_
       do j=1,this%nspecies_
          dAmdT=moleFrac(i)*moleFrac(j)*A(i,j)*G(i,j)+dAmdT
       enddo
    enddo 
    
    dAmdT=-1._PREC*dAmdT/T

  end subroutine eosPR_dAmdT
  
  !!============================================================================
  !!
  !! Sub: eosPR_d2AmdT2
  !!    Returns d2AmdT2 [Appendix B]
  !!
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!
  !!  d2AmdT2 [real, out] 
  !===========================================================================  
  subroutine eosPR_d2AmdT2(this,T,massFrac,d2AmdT2)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: T
    real(PREC), intent(out) :: d2AmdT2
    real(PREC) :: moleFrac(this%nspecies_)
    real(PREC)  :: A(this%nspecies_,this%nspecies_)
    real(PREC)  :: B(this%nspecies_)
    real(PREC)  :: C(this%nspecies_,this%nspecies_)
    real(PREC)  :: G(this%nspecies_,this%nspecies_)
    real(PREC)  :: Tc(this%nspecies_,this%nspecies_)
    real(PREC)  :: pc(this%nspecies_,this%nspecies_)
    
    integer :: i,j
    
    call eosPR_moleFrac(this,massFrac,moleFrac)
    call eosPR_Para(this,T,A,B,C,G,pc,Tc)
    
    d2AmdT2=0._PREC

    do i=1,this%nspecies_
       do j=1,this%nspecies_
          d2AmdT2=moleFrac(i)*moleFrac(j)*C(i,j)*(1+C(i,j))*(Tc(i,j)/pc(i,j))* &
                     sqrt(Tc(i,j)/T)+d2AmdT2
       enddo
    enddo 
    
    d2AmdT2=(0.457236_PREC*R_UNIVERSAL**2)*d2AmdT2/(2._PREC*T)

  end subroutine eosPR_d2AmdT2
  
  !!============================================================================
  !!
  !! Sub: eosPR_dAmdX
  !!    Returns dAmdX [Appendix B]
  !!
  !! eosPR_dAmdx(this,massFrac)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  
  !!  dAmdX [real, out] - [(m3/kmol)2]
  !===========================================================================  
  subroutine eosPR_dAmdX(this,T,moleFrac,dAmdX)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: moleFrac(this%nspecies_)
    real(PREC), intent(in) :: T
    real(PREC), intent(out) :: dAmdX(this%nspecies_)
    real(PREC)  :: A(this%nspecies_,this%nspecies_)
    real(PREC)  :: B(this%nspecies_)
    real(PREC)  :: C(this%nspecies_,this%nspecies_)
    real(PREC)  :: G(this%nspecies_,this%nspecies_)
    real(PREC)  :: Tc(this%nspecies_,this%nspecies_)
    real(PREC)  :: pc(this%nspecies_,this%nspecies_)
    
    integer :: i,j
    
    call eosPR_Para(this,T,A,B,C,G,pc,Tc)
    
    dAmdX(1:this%nspecies_)=0._PREC
    
    do i=1,this%nspecies_
      do j=1,this%nspecies_
        dAmdX(i)=moleFrac(j)*A(i,j)+dAmdX(i)
      enddo
      dAmdX(i)=2._PREC*dAmdX(i)
    enddo

  end subroutine eosPR_dAmdX
  
  !!============================================================================
  !!
  !! Sub: eosPR_d2AmdXdT
  !!    Returns d2AmdXdT [Appendix B]
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  
  !!  d2AmdXdT [vec, out] - 
  !===========================================================================  
  subroutine eosPR_d2AmdXdT(this,T,moleFrac,d2AmdXdT)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: moleFrac(this%nspecies_)
    real(PREC), intent(in) :: T
    real(PREC), intent(out) :: d2AmdXdT(this%nspecies_)
    real(PREC)  :: A(this%nspecies_,this%nspecies_)
    real(PREC)  :: B(this%nspecies_)
    real(PREC)  :: C(this%nspecies_,this%nspecies_)
    real(PREC)  :: G(this%nspecies_,this%nspecies_)
    real(PREC)  :: Tc(this%nspecies_,this%nspecies_)
    real(PREC)  :: pc(this%nspecies_,this%nspecies_)
    
    integer :: i,j
    
    call eosPR_Para(this,T,A,B,C,G,pc,Tc)
    
    d2AmdXdT(1:this%nspecies_)=0._PREC
    
    do i=1,this%nspecies_
      do j=1,this%nspecies_
        d2AmdXdT(i)=moleFrac(j)*A(i,j)*G(i,j)+d2AmdXdT(i)
      enddo
      d2AmdXdT(i)=-2._PREC*d2AmdXdT(i)/T
    enddo
  
  end subroutine eosPR_d2AmdXdT
  
  !!============================================================================
  !!
  !! Sub: eosPR_K1
  !!    Returns the K1 parameter [Eq. (2.13)]
  !!
  !! eosPR_K1(this,T,D,massFrac,K1)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  T [real, in] - Temperature [K]
  !!  D [real, in] - Density [kg/m3]  
  !!  K1 [real, out] - K1 parameter    
  !===========================================================================  
  subroutine eosPR_K1(this,T,D,massFrac,K1)
    implicit none
    type(eosPR), intent(inout) :: this 
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: T,D 
    real(PREC), intent(out) :: K1
    real(PREC) :: molarVolume
    real(PREC) :: moleFrac(this%nspecies_)
    real(PREC) :: Bm,Am
    !
    call eosPR_molarVolume(this,D,massFrac,molarVolume)
    call eosPR_moleFrac(this,massFrac,moleFrac)
    call eosPR_Bm(this,T,moleFrac,Bm)
    !p=R_UNIVERSAL*T/(molarVolume-Bm)-Am/(molarVolume**2+2*molarVolume*Bm-Bm**2)
    K1=log((molarVolume+(1._PREC-sqrt(2._PREC))*Bm)/  &
           (molarVolume+(1._PREC+sqrt(2._PREC))*Bm))
    K1=K1/(2._PREC*sqrt(2._PREC)*Bm)
    
  end subroutine eosPR_K1
  
  !!============================================================================
  !!
  !! Sub: eosPR_dpdT
  !!    Returns the dp/dT [Eq. (2.11)]
  !!
  !! eosPR_pressure(this,massFrac)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  T [real, in] - Temperature [K]
  !!  D [real, in] - Density [kg/m3]  
  !!  dpdT [real, out] - [pa/K]    
  !===========================================================================  
  subroutine eosPR_dpdT(this,T,D,massFrac,dpdT)
    implicit none
    type(eosPR), intent(inout) :: this 
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: T,D 
    real(PREC), intent(out) :: dpdT
    real(PREC) :: molarVolume
    real(PREC) :: moleFrac(this%nspecies_)
    real(PREC) :: Bm,dAmdT
    !
    call eosPR_molarVolume(this,D,massFrac,molarVolume)
    call eosPR_moleFrac(this,massFrac,moleFrac)
    call eosPR_dAmdT(this,T,moleFrac,dAmdT)
    call eosPR_Bm(this,T,moleFrac,Bm)
    dpdT=R_UNIVERSAL/(molarVolume-Bm)-dAmdT/(molarVolume**2+2._PREC*molarVolume*Bm-Bm**2)
   
  end subroutine eosPR_dpdT
   
  !!============================================================================
  !!
  !! Sub: eosPR_dpdv
  !!    Returns the dp/dv [Eq. (2.12)]
  !!
  !! eosPR_pressure(this,massFrac)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  T [real, in] - Temperature [K]
  !!  D [real, in] - Density [kg/m3]  
  !!  dpdv [real, out] - [pa.kmol/m3]    
  !===========================================================================  
  subroutine eosPR_dpdv(this,T,D,massFrac,dpdv)
    implicit none
    type(eosPR), intent(inout) :: this 
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: T,D 
    real(PREC), intent(out) :: dpdv
    real(PREC) :: molarVolume
    real(PREC) :: moleFrac(this%nspecies_)
    real(PREC) :: Bm,Am
    !
    call eosPR_molarVolume(this,D,massFrac,molarVolume)
    call eosPR_moleFrac(this,massFrac,moleFrac)
    call eosPR_Am(this,T,moleFrac,Am)
    call eosPR_Bm(this,T,moleFrac,Bm)
    dpdv=(-R_UNIVERSAL*T/((molarVolume-Bm)**2))* &
          (1._PREC-2._PREC*Am/(R_UNIVERSAL*T*(molarVolume+Bm)* &
            ((molarVolume/(molarVolume-Bm)+Bm/(molarVolume+Bm)))**2))
  
  end subroutine eosPR_dpdv
   
  !!============================================================================
  !!
  !! Sub: eosPR_h0
  !!    Returns h0 - Reference molar specific enthalpy [h0=sum h0_j*X_j]
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  
  !!  h0 [real, out] - [J/kmol] 
  !===========================================================================  
  subroutine eosPR_h0(this,T,massFrac,h0)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: T
    real(PREC), intent(out) :: h0
    real(PREC) :: h0_species(this%nspecies_)
    real(PREC) :: moleFrac(this%nspecies_)
    
    integer :: i
    
    call eosPR_moleFrac(this,massFrac,moleFrac)    
    call eosPR_h0_species(this,T,h0_species)

    h0=0._PREC
    do i=1,this%nspecies_
          h0=moleFrac(i)*h0_species(i)+h0
          !h0=moleFrac(i)*this%DELGF_(i)+h0 
    enddo
    

  end subroutine eosPR_h0
  
  !!============================================================================
  !!
  !! Sub: eosPR_Cp0
  !!    Returns Cp0 - Reference molar heat capacity at constant pressure [Eq. (2.30)]
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  
  !!  Cp0 [real, out] - [J/kmol.K] 
  !===========================================================================  
  subroutine eosPR_Cp0(this,T,massFrac,Cp0)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(in) :: T
    real(PREC), intent(out) :: Cp0
    real(PREC) :: molarEnthalpy(this%nspecies_)
    real(PREC) :: moleFrac(this%nspecies_)
    real(PREC) :: Cp0_species(this%nspecies_)
    
    integer :: i
    
    call eosPR_moleFrac(this,massFrac,moleFrac)

    Cp0_species(1)=this%molecularMass_(1)*656.72_PREC*1.071_PREC*T**(0.071_PREC)
    Cp0_species(2)=this%molecularMass_(2)*27.877_PREC*1.6414_PREC*T**(0.6414_PREC)
    
    Cp0=0._PREC
    do i=1,this%nspecies_
          Cp0=moleFrac(i)*Cp0_species(i)+Cp0 
    enddo

  end subroutine eosPR_Cp0
  
  !!============================================================================
  !!
  !! Sub: eosPR_moleFrac
  !!    Returns the mole fractions using mass fractions [X_j=m/m_j*Y_j]
  !!
  !! eosPR_moleFrac(this,massFrac,moleFrac)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  moleFrac [real(1:nspecies), out] - mole fraction of all chemical species  
  !===========================================================================  
  subroutine eosPR_moleFrac(this,massFrac,moleFrac)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(out) :: moleFrac(this%nspecies_)
    real(PREC) :: mixMolecularMass
    !
    call eosPR_mixMolecularMass(this,massFrac,mixMolecularMass)
    moleFrac(1:this%nspecies_)=mixMolecularMass/this%molecularMass_(1:this%nspecies_)*&
                                                           massFrac(1:this%nspecies_)

  end subroutine eosPR_moleFrac
  
  !!============================================================================
  !!
  !! Sub: eosPR_mixMolecularMass
  !!    Returns mixture molecular mass using mass fractions [m=1/sum(Y_j/m_j)]
  !!
  !! eosPR_mixMolecularMass(this,massFrac,mixMolecularMass)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  mixMolecularMass [real, out] - mixture molecular mass [kg/kmol]
  !===========================================================================  
  subroutine eosPR_mixMolecularMass(this,massFrac,mixMolecularMass)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(out) :: mixMolecularMass
    integer :: j
    !
    mixMolecularMass=0._PREC
    do j=1,this%nspecies_
       mixMolecularMass=massFrac(j)/this%molecularMass_(j)+mixMolecularMass
    enddo
    mixMolecularMass=1._PREC/mixMolecularMass
    
  end subroutine eosPR_mixMolecularMass
  
  !!============================================================================
  !!
  !! Sub: eosPR_molecularMass
  !!    Returns mixture array of molecular masses
  !!
  !! eosPR_molecularMass(this,molecularMass)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  molecularMass [real(1:nspecies), out] - molecular mass of all chemical species
  !===========================================================================  
  subroutine eosPR_molecularMass(this,molecularMass)
    implicit none
    type(eosPR), intent(in) :: this
    real(PREC), intent(out) :: molecularMass(this%nspecies_)
    !
    molecularMass=this%molecularMass_
  end subroutine eosPR_molecularMass

  !!============================================================================
  !!
  !! Sub: eosPR_molarVolume
  !!    Returns molar specific volume  using density and mass fractions 
  !!    [v=m/density]
  !!
  !! eosPR_molarVolume(this,D,massFrac,molarVolume)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  D [real, in] - density [kg/m3]
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species
  !!  molarVolume [real, out] - molar specific volume [m3/kmol]  
  !===========================================================================  
  subroutine eosPR_molarVolume(this,D,massFrac,molarVolume)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: D
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(out) :: molarVolume
    real(PREC) :: mixMolecularMass
    !
    call eosPR_mixMolecularMass(this,massFrac,mixMolecularMass)
    molarVolume=mixMolecularMass/D 
  end subroutine eosPR_molarVolume

  !!============================================================================
  !!
  !! Sub: eosPR_density
  !!    Returns density  using molar specific volume and mass fractions
  !!    [density=m/v]
  !!
  !! eosPR_density(this,molarVolume,massFrac,D)
  !!
  !!  this [type(eosPR), inout] - eosPR object
  !!  molarVolume [real, out] - molar specific volume [m3/kmol]
  !!  massFrac [real(1:nspecies), in] - mass fraction of all chemical species  
  !!  D [real, out] - density [kg/m3]
  !===========================================================================  
  subroutine eosPR_density(this,molarVolume,massFrac,D)
    implicit none
    type(eosPR), intent(inout) :: this
    real(PREC), intent(in) :: molarVolume
    real(PREC), intent(in) :: massFrac(this%nspecies_)
    real(PREC), intent(out) :: D
    real(PREC) :: mixMolecularMass
    !
    call eosPR_mixMolecularMass(this,massFrac,mixMolecularMass)
    D=mixMolecularMass/molarVolume
  end subroutine eosPR_density

end module mod_eosPR      
