module computer
    use config
    use utilities    
    use integral
    implicit none
    
    type :: GenericComputer
        !*****************************************************************************************************************
        !"""
        !Mother class of Forecaster and Analyser. Implements methods and members that are common to both.
        !"""
        !*****************************************************************************************************************
        class(ComputationConfig), allocatable :: cfg
        class(legendre_polys_type), allocatable :: algo_legendre_polys
        integer :: algo_nb_realisations, algo_seed
        
    contains
        procedure :: init_GenericComputer
        procedure :: compute_Ab
    end type
    
    
contains    
!==========================================================================================================================  
    subroutine init_GenericComputer(self, config, legendre_polys, nb_realisations, seed)
    !*****************************************************************************************************************
    !"""
    !Constructor of GenericComputer. Sets the args as members.
    !
    !:param algo: Algorithm object
    !:type algo: Algo
    !"""
    !*****************************************************************************************************************
        class(ComputationConfig), intent(in) :: config
        class(legendre_polys_type), intent(in) :: legendre_polys
        class(GenericComputer), intent(inout) :: self  
        integer :: nb_realisations, seed
        
        self.cfg = config
        self.algo_legendre_polys = legendre_polys
        self.algo_nb_realisations = nb_realisations
        self.algo_seed = seed
    end subroutine    
!==========================================================================================================================
    
!==========================================================================================================================  
    subroutine compute_Ab(self, input_core_state, AbT)
    !*****************************************************************************************************************
    !"""
    !Compute A(b) the operator DivH(uBr) in the spectral domain.
    !
    !:param input_core_state: 1D CoreState storing a minima the magnetic field coefficients.
    !:type input_core_state: 1D CoreState
    !:return: The computed A(b)
    !:rtype: 2D numpy.ndarray (dim: Nsv x Nu2)
    !"""
    !*****************************************************************************************************************
        class(input_core_state_type), intent(in) :: input_core_state
        class(GenericComputer), intent(in) :: self
        real(kind=8), allocatable, intent(out) :: AbT(:,:)
        real(kind=8), allocatable :: gauss_th(:), gauss_weights(:)
        integer :: tmax, LLb, LLu, LLsv
        real(kind=8), allocatable :: lpb(:,:), dlpb(:,:), d2lpb(:,:)
        real(kind=8), allocatable :: lpu(:,:), dlpu(:,:), d2lpu(:,:)
        real(kind=8), allocatable :: lpsv(:,:), dlpsv(:,:), d2lpsv(:,:)
        real(kind=8), allocatable :: Ab(:,:)
        
        allocate(gauss_th, source=self.algo_legendre_polys.thetas)
        allocate(gauss_weights, source=self.algo_legendre_polys.weights)
        tmax = SIZE(gauss_th)
        
        allocate(lpb, source=self.algo_legendre_polys.MF(1, :, :))
        allocate(dlpb, source=self.algo_legendre_polys.MF(2, :, :))
        allocate(d2lpb, source=self.algo_legendre_polys.MF(3, :, :))
        LLb = SIZE(lpb, 1)
        
        allocate(lpu, source=self.algo_legendre_polys.U(1, :, :))
        allocate(dlpu, source=self.algo_legendre_polys.U(2, :, :))
        allocate(d2lpu, source=self.algo_legendre_polys.U(3, :, :))
        LLu = SIZE(lpu, 1)
        
        allocate(lpsv, source=self.algo_legendre_polys.SV(1, :, :))
        allocate(dlpsv, source=self.algo_legendre_polys.SV(2, :, :))
        allocate(d2lpsv, source=self.algo_legendre_polys.SV(3, :, :))
        LLsv = SIZE(lpsv, 1)
        
        !# Function radmats computes A(b)^T
        allocate(Ab(input_core_state.Nu2, input_core_state.Nsv), source=0.0d0)
        allocate(AbT(input_core_state.Nsv, input_core_state.Nu2), source=0.0d0)
        call radmats(gauss_th, gauss_weights, lpsv, dlpsv, d2lpsv, &
                     lpu, dlpu, d2lpu, lpb, dlpb, d2lpb, &
                     Ab, input_core_state.B, &
                     input_core_state.Lsv, input_core_state.Lu, input_core_state.Lb, &
                     input_core_state.Nsv, input_core_state.Nu2 / 2, input_core_state.Nb, &
                     LLsv, LLu, LLb, tmax)
        AbT = TRANSPOSE(Ab)
    end subroutine
!==========================================================================================================================
    
end module