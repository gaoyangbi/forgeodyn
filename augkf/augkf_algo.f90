module augkf_algo
    use utilities
    use generic_algo
    use reads
    use priors
    use pca
    use common
    use mpi
    use forecaster
    use analyser
    use corestate
    implicit none
    
    type, extends(GenericAlgo) :: AugkfAlgo
        integer, allocatable :: attributed_models(:)
        integer :: seed
        class(NormedPCAOperator), allocatable :: pcaU_operator
        class(set_prior_type), allocatable :: avg_prior
        class(cov_prior_type), allocatable :: cov_prior
        class(legendre_polys_type), allocatable :: legendre_polys
        class(AugkfForecasterAR1), allocatable :: forecaster_1
        class(AugkfForecasterAR3), allocatable :: forecaster_3
        class(AugkfAnalyserAR1), allocatable :: analyser_1
        class(AugkfAnalyserAR3), allocatable :: analyser_3
    contains
        procedure :: init_AugkfAlgo, check_PCA, create_forecaster, create_analyser
        procedure :: init_corestates, analysis_step_algo, forecast_step_algo
        procedure :: is_equ
        procedure :: extract_prior_and_covariances
        procedure :: check_coef_size
        procedure :: set_prior_size, U_ER_to_Z2
    
    
    end type
    
contains
    
    subroutine create_augkf_algo(config, nb_realisations, seed, attributed_models, algo)
    !"""
    !Factory function that returns the AugKF algo
    !
    !:param config: Configuration of the algo
    !:type config: ComputationConfig
    !:param nb_realisations: number of realisations to consider in the algo
    !:type nb_realisations: int
    !:return: Algorithm object
    !:rtype: AugkfAlgo
    !"""
        class(ComputationConfig), intent(in) :: config
        integer, intent(in) :: nb_realisations, seed
        integer, intent(in) :: attributed_models(:)
        class(AugkfAlgo), intent(out) :: algo
        
        call algo.init_AugkfAlgo(config, nb_realisations, seed, attributed_models)
    end subroutine
    
!==========================================================================================================================    
    function is_equ(self, other) result(are_all_equal)
    !*****************************************************************************************************************
    !"""
    !Implementation of equality for two different AugkfAlgos
    !
    !:param: other
    !"""
    !*****************************************************************************************************************
        class(AugkfAlgo), intent(in) :: self, other
        logical :: are_all_equal
        
        
        are_all_equal = .true.
        write(10,'(A)') 'testing equalities of algo'
        !# loop on all items of the two instantiation of the classes
        
    contains
        subroutine print_equal_message(eq, key1, key2)
            logical, intent(in) :: eq
            character(len=*), intent(in) :: key1, key2
            if (eq) then
                ! do nothing
            else
                write(10,'(A)') "Not equal: "//trim(key1)//" is not equal in both algos. "
            end if
        
        end subroutine        
    end function  
!==========================================================================================================================
    
!==========================================================================================================================    
    subroutine init_AugkfAlgo(self, cfg, nb_realisations, seed, attributed_models)
    !*****************************************************************************************************************
    !"""
    !Base class that defines the interface of Algo for it to be used in run.py
    !"""
    !*****************************************************************************************************************
        class(AugkfAlgo), intent(inout) :: self
        class(ComputationConfig), intent(in) :: cfg
        integer, intent(in) :: nb_realisations, seed
        integer, intent(in) :: attributed_models(:)
        
        
        call self.init_GenericAlgo(cfg, nb_realisations)
        self.attributed_models = attributed_models
        allocate(self.avg_prior, self.cov_prior)
        call self.extract_prior_and_covariances(self.avg_prior, self.cov_prior)
        allocate(self.legendre_polys)
        call compute_legendre_polys(cfg.Nth_legendre, cfg.Lb, cfg.Lu, cfg.Lsv, self.legendre_polys)
        call self.create_forecaster()
        self.seed = seed
        call self.create_analyser()   
    end subroutine    
!==========================================================================================================================
    
!==========================================================================================================================
    function check_PCA(self) result(log)
        class(AugkfAlgo) :: self
        logical :: log
        
        if (self.config.pca == 1) then
            log = .true.
        else
            log = .false.
        end if
    end function
!==========================================================================================================================
    
!==========================================================================================================================    
    subroutine create_forecaster(self)
    !*****************************************************************************************************************
    !"""
    !Factory method to create the forecaster.
    !
    !:return: AugkfForecaster
    !"""
    !*****************************************************************************************************************
        class(AugkfAlgo), intent(inout) :: self
        
        if (trim(self.config.AR_type) == 'AR3') then
            allocate(self.forecaster_3)
            call self.forecaster_3.init_AugkfForecasterAR(self.config, self.legendre_polys)
        else
            allocate(self.forecaster_1)
            call self.forecaster_1.init_AugkfForecasterAR(self.config, self.legendre_polys)
        end if   
    end subroutine    
!==========================================================================================================================
    
!==========================================================================================================================    
    subroutine create_analyser(self)
    !*****************************************************************************************************************
    !"""
    !Factory method to create the analyser.
    !
    !:return: AugkfAnalyser
    !"""
    !*****************************************************************************************************************
        class(AugkfAlgo), intent(inout) :: self
        
        if (trim(self.config.AR_type) == 'AR3') then
            allocate(self.analyser_3)
            call self.analyser_3.init_AugkfAnalyserAR(self.config, self.legendre_polys, self.nb_realisations, self.seed)
        else
            allocate(self.analyser_1)
            call self.analyser_1.init_AugkfAnalyserAR(self.config, self.legendre_polys, self.nb_realisations, self.seed)
        end if   
    end subroutine    
!==========================================================================================================================

!==========================================================================================================================    
    subroutine init_corestates(self, random_state, computed_states, forecast_states, analysed_states, misfits, Z_AR3)
    !*****************************************************************************************************************
    !"""
    !Sets up the corestates needed for the AugKF algorithm.
    !Returns CoreStates of the adequate form and initialisation to perform the AugKF.
    !
    !:param random_state: Random state to use for normal draw (only used when init from noised priors)
    !:type random_state: np.random.RandomState
    !:return: computed_states, forecast_states, analysed_states, misfits, Z_AR3
    !:rtype: CoreState, CoreState, CoreState, CoreState, 3D numpy array (N_real x 3 x Ncoef) or None (if not AR3)
    !"""
    !*****************************************************************************************************************
        class(AugkfAlgo), intent(inout) :: self
        real(kind=8), intent(in) :: random_state
        type(corestate_measures_type) :: corestate_measures
        class(key_measures), allocatable :: computed_states_data(:), analysed_states_data(:), misfits_data(:)
        type(CoreState_type), intent(out) :: computed_states, forecast_states, analysed_states, misfits
        REAL(kind=8), allocatable, intent(out) :: Z_AR3(:,:,:)
        character(len=100) :: str_init
        
        !# Define the measures used for the computation and their number of coeffs
        corestate_measures.MF = self.config.Nb()
        corestate_measures.U = self.config.Nu2()
        corestate_measures.SV = self.config.Nsv()
        corestate_measures.ER = self.config.Nsv()
        corestate_measures.Z = self.config.Nz()
        
        !# Add derivatives needed to AR3 to measures
        if (TRIM(self.config.AR_type) == "AR3") then
            corestate_measures.dUdt = self.config.Nu2()
            corestate_measures.d2Udt2 = self.config.Nu2()
            corestate_measures.dERdt = self.config.Nsv()
            corestate_measures.d2ERdt2 = self.config.Nsv()
        end if
        
        !# Add shear measure ('S') if do_shear = 1
        if (self.config.do_shear == 1) then
            corestate_measures.S = self.config.Nu2()
        end if
        
        if ((TRIM(self.config.AR_type) == "AR3") .and.(self.config.do_shear == 1)) then
            allocate(computed_states_data(10))
            allocate(analysed_states_data(SIZE(computed_states_data)))
            
            computed_states_data(1).key = 'MF'
            analysed_states_data(1).key = 'MF'
            allocate(computed_states_data(1).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.MF), source = 0.0d0)
            allocate(analysed_states_data(1).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.MF), source = 0.0d0)
            computed_states_data(2).key = 'U'
            analysed_states_data(2).key = 'U'
            allocate(computed_states_data(2).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.U), source = 0.0d0)
            allocate(analysed_states_data(2).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.U), source = 0.0d0)
            computed_states_data(3).key = 'SV'
            analysed_states_data(3).key = 'SV'
            allocate(computed_states_data(3).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.SV), source = 0.0d0)
            allocate(analysed_states_data(3).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.SV), source = 0.0d0)
            computed_states_data(4).key = 'ER'
            analysed_states_data(4).key = 'ER'
            allocate(computed_states_data(4).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.ER), source = 0.0d0)
            allocate(analysed_states_data(4).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.ER), source = 0.0d0)
            computed_states_data(5).key = 'Z'
            analysed_states_data(5).key = 'Z'
            allocate(computed_states_data(5).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.Z), source = 0.0d0)
            allocate(analysed_states_data(5).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.Z), source = 0.0d0)
            computed_states_data(6).key = 'dUdt'
            analysed_states_data(6).key = 'dUdt'
            allocate(computed_states_data(6).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.dUdt), source = 0.0d0)
            allocate(analysed_states_data(6).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.dUdt), source = 0.0d0)
            computed_states_data(7).key = 'd2Udt2'
            analysed_states_data(7).key = 'd2Udt2'
            allocate(computed_states_data(7).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.d2Udt2), source = 0.0d0)
            allocate(analysed_states_data(7).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.d2Udt2), source = 0.0d0)
            computed_states_data(8).key = 'dERdt'
            analysed_states_data(8).key = 'dERdt'
            allocate(computed_states_data(8).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.dERdt), source = 0.0d0)
            allocate(analysed_states_data(8).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.dERdt), source = 0.0d0)
            computed_states_data(9).key = 'd2ERdt2'
            analysed_states_data(9).key = 'd2ERdt2'
            allocate(computed_states_data(9).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.d2ERdt2), source = 0.0d0)
            allocate(analysed_states_data(9).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.d2ERdt2), source = 0.0d0)
            computed_states_data(10).key = 'S'
            analysed_states_data(10).key = 'S'
            allocate(computed_states_data(10).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.S), source = 0.0d0)
            allocate(analysed_states_data(10).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.S), source = 0.0d0)
        elseif ((TRIM(self.config.AR_type) == "AR3") .and.(self.config.do_shear == 0)) then
            allocate(computed_states_data(9))
            allocate(analysed_states_data(SIZE(computed_states_data)))
            
            computed_states_data(1).key = 'MF'
            analysed_states_data(1).key = 'MF'
            allocate(computed_states_data(1).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.MF), source = 0.0d0)
            allocate(analysed_states_data(1).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.MF), source = 0.0d0)
            computed_states_data(2).key = 'U'
            analysed_states_data(2).key = 'U'
            allocate(computed_states_data(2).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.U), source = 0.0d0)
            allocate(analysed_states_data(2).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.U), source = 0.0d0)
            computed_states_data(3).key = 'SV'
            analysed_states_data(3).key = 'SV'
            allocate(computed_states_data(3).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.SV), source = 0.0d0)
            allocate(analysed_states_data(3).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.SV), source = 0.0d0)
            computed_states_data(4).key = 'ER'
            analysed_states_data(4).key = 'ER'
            allocate(computed_states_data(4).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.ER), source = 0.0d0)
            allocate(analysed_states_data(4).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.ER), source = 0.0d0)
            computed_states_data(5).key = 'Z'
            analysed_states_data(5).key = 'Z'
            allocate(computed_states_data(5).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.Z), source = 0.0d0)
            allocate(analysed_states_data(5).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.Z), source = 0.0d0)
            computed_states_data(6).key = 'dUdt'
            analysed_states_data(6).key = 'dUdt'
            allocate(computed_states_data(6).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.dUdt), source = 0.0d0)
            allocate(analysed_states_data(6).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.dUdt), source = 0.0d0)
            computed_states_data(7).key = 'd2Udt2'
            analysed_states_data(7).key = 'd2Udt2'
            allocate(computed_states_data(7).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.d2Udt2), source = 0.0d0)
            allocate(analysed_states_data(7).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.d2Udt2), source = 0.0d0)
            computed_states_data(8).key = 'dERdt'
            analysed_states_data(8).key = 'dERdt'
            allocate(computed_states_data(8).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.dERdt), source = 0.0d0)
            allocate(analysed_states_data(8).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.dERdt), source = 0.0d0)
            computed_states_data(9).key = 'd2ERdt2'
            analysed_states_data(9).key = 'd2ERdt2'
            allocate(computed_states_data(9).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.d2ERdt2), source = 0.0d0)
            allocate(analysed_states_data(9).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.d2ERdt2), source = 0.0d0)
        elseif ((TRIM(self.config.AR_type) .ne. "AR3") .and.(self.config.do_shear == 1)) then
            allocate(computed_states_data(6))
            allocate(analysed_states_data(SIZE(computed_states_data)))
            
            computed_states_data(1).key = 'MF'
            analysed_states_data(1).key = 'MF'
            allocate(computed_states_data(1).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.MF), source = 0.0d0)
            allocate(analysed_states_data(1).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.MF), source = 0.0d0)
            computed_states_data(2).key = 'U'
            analysed_states_data(2).key = 'U'
            allocate(computed_states_data(2).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.U), source = 0.0d0)
            allocate(analysed_states_data(2).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.U), source = 0.0d0)
            computed_states_data(3).key = 'SV'
            analysed_states_data(3).key = 'SV'
            allocate(computed_states_data(3).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.SV), source = 0.0d0)
            allocate(analysed_states_data(3).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.SV), source = 0.0d0)
            computed_states_data(4).key = 'ER'
            analysed_states_data(4).key = 'ER'
            allocate(computed_states_data(4).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.ER), source = 0.0d0)
            allocate(analysed_states_data(4).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.ER), source = 0.0d0)
            computed_states_data(5).key = 'Z'
            analysed_states_data(5).key = 'Z'
            allocate(computed_states_data(5).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.Z), source = 0.0d0)
            allocate(analysed_states_data(5).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.Z), source = 0.0d0)
            computed_states_data(6).key = 'S'
            analysed_states_data(6).key = 'S'
            allocate(computed_states_data(6).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.S), source = 0.0d0)
            allocate(analysed_states_data(6).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.S), source = 0.0d0)
        elseif ((TRIM(self.config.AR_type) .ne. "AR3") .and.(self.config.do_shear == 0)) then
            allocate(computed_states_data(5))
            allocate(analysed_states_data(SIZE(computed_states_data)))
            
            computed_states_data(1).key = 'MF'
            analysed_states_data(1).key = 'MF'
            allocate(computed_states_data(1).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.MF), source = 0.0d0)
            allocate(analysed_states_data(1).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.MF), source = 0.0d0)
            computed_states_data(2).key = 'U'
            analysed_states_data(2).key = 'U'
            allocate(computed_states_data(2).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.U), source = 0.0d0)
            allocate(analysed_states_data(2).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.U), source = 0.0d0)
            computed_states_data(3).key = 'SV'
            analysed_states_data(3).key = 'SV'
            allocate(computed_states_data(3).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.SV), source = 0.0d0)
            allocate(analysed_states_data(3).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.SV), source = 0.0d0)
            computed_states_data(4).key = 'ER'
            analysed_states_data(4).key = 'ER'
            allocate(computed_states_data(4).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.ER), source = 0.0d0)
            allocate(analysed_states_data(4).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.ER), source = 0.0d0)
            computed_states_data(5).key = 'Z'
            analysed_states_data(5).key = 'Z'
            allocate(computed_states_data(5).measure_data(size(self.attributed_models), self.config.nb_forecasts(), corestate_measures.Z), source = 0.0d0)
            allocate(analysed_states_data(5).measure_data(size(self.attributed_models), self.config.nb_analyses(), corestate_measures.Z), source = 0.0d0)
        end if
        
        !# Build the array of computed states and analysed_states
        call computed_states.init_CoreState(computed_states_data)
        call analysed_states.init_CoreState(analysed_states_data)
        
        !# Initialize the core state at t=0
        if (trim(self.config.core_state_init) == 'from_file') then
            call computed_states.initialise_from_file(self.config, self.attributed_models, Z_AR3)
            str_init = self.config.init_file
        else
            if (trim(self.config.AR_type) == "AR3") then
                call computed_states.initialise_from_noised_priors(random_state, self.config, self.attributed_models, self.avg_prior, self.cov_prior, self.analyser_3, self.pcaU_operator, Z_AR3)
            else
                call computed_states.initialise_from_noised_priors(random_state, self.config, self.attributed_models, self.avg_prior, self.cov_prior, self.analyser_1, self.pcaU_operator, Z_AR3)
            end if    
            str_init = 'normal draw around average priors'
        end if
        
        write(10, '(A,A)') 'Computed states initialised from ', TRIM(str_init)
        write(*, '(A,A)') 'Computed states initialised from ', TRIM(str_init)
        
        !# Create the array storing the result of only forecasts (also copies the value at t=0)
        forecast_states = computed_states
        
        !# Create the CoreState (1 realisation and 1 coef (max_degree forced to 0)) storing the misfits of analyses
        allocate(misfits_data(2))
        
        misfits_data(1).key = 'MF'
        allocate(misfits_data(1).measure_data(1, self.config.nb_analyses(), 1), source = 0.0d0)
        misfits_data(2).key = 'SV'
        allocate(misfits_data(2).measure_data(1, self.config.nb_analyses(), 1), source = 0.0d0)
        call misfits.init_CoreState(misfits_data)        
        misfits.max_degrees_(1).max_d = 0
        misfits.max_degrees_(2).max_d = 0
        write(10, '(A)') "AugKF CoreStates ready !"
        write(*, '(A)') "AugKF CoreStates ready !"
    end subroutine    
!==========================================================================================================================
    
!==========================================================================================================================    
    subroutine analysis_step_algo(self, input_core_state, analysis_time)
    !*****************************************************************************************************************
    !"""
    !Sets up the corestates needed for the AugKF algorithm.
    !Returns CoreStates of the adequate form and initialisation to perform the AugKF.
    !
    !:param random_state: Random state to use for normal draw (only used when init from noised priors)
    !:type random_state: np.random.RandomState
    !:return: computed_states, forecast_states, analysed_states, misfits, Z_AR3
    !:rtype: CoreState, CoreState, CoreState, CoreState, 3D numpy array (N_real x 3 x Ncoef) or None (if not AR3)
    !"""
    !*****************************************************************************************************************
        class(AugkfAlgo), intent(inout) :: self
        real(8), intent(in) :: input_core_state
        real(8), intent(in) :: analysis_time
        
        
    end subroutine    
!==========================================================================================================================
    
!==========================================================================================================================    
    subroutine forecast_step_algo(self, input_core_state, random_state)
    !*****************************************************************************************************************
    !"""
    !Sets up the corestates needed for the AugKF algorithm.
    !Returns CoreStates of the adequate form and initialisation to perform the AugKF.
    !
    !:param random_state: Random state to use for normal draw (only used when init from noised priors)
    !:type random_state: np.random.RandomState
    !:return: computed_states, forecast_states, analysed_states, misfits, Z_AR3
    !:rtype: CoreState, CoreState, CoreState, CoreState, 3D numpy array (N_real x 3 x Ncoef) or None (if not AR3)
    !"""
    !*****************************************************************************************************************
        class(AugkfAlgo), intent(inout) :: self
        real(8), intent(in) :: input_core_state
        integer, intent(in), optional :: random_state
        
        
    end subroutine    
!==========================================================================================================================
    
!==========================================================================================================================    
    subroutine extract_prior_and_covariances(self, avg_prior, cov_prior)
    !*****************************************************************************************************************
    !"""
    !Extracts the priors from the files in the config prior directory. Also sets
    !the covariance matrices by computation from priors or by reading them.
    !
    !:return: average priors and covariances matrices as dictionaries.
    !:rtype: dict, dict
    !"""
    !*****************************************************************************************************************
        class(AugkfAlgo), intent(inout) :: self
        real(8) :: dt_f
        integer :: Nz, Nb, Nsv, Nuz
        character(len=50) :: AR_type
        type(prior_data), allocatable :: prior_data_obj(:)
        character(len=50) :: measures(4)
        character(len=50), allocatable :: tag_list(:)
        character(len=50) :: t
        integer :: i, j, k, rank, ierr        
        type(cov_prior_type), intent(out) :: cov_prior
        type(set_prior_type), allocatable :: set_prior_obj(:)
        type(set_prior_type), intent(out) :: avg_prior
        real(kind=8) :: dt_sampling, dt_prior, blackman_w
        integer :: Nt
        real(kind=8), allocatable :: matrix_1(:), matrix_2(:,:), matrix_3(:,:)
        real(kind=8), allocatable :: matrix_temp(:,:), matrix_temp2(:,:)
        real(kind=8), allocatable :: Z(:,:), dZ(:,:), d2Z(:,:), d3Z(:,:)
        real(kind=8), allocatable :: u(:,:), du(:,:), d2u(:,:), d3u(:,:)
        real(kind=8), allocatable :: e(:,:), de(:,:), d2e(:,:), d3e(:,:)
        real(kind=8), allocatable :: b(:,:), db(:,:), d2b(:,:), d3b(:,:)
        real(kind=8), allocatable :: u_e(:,:), du_de(:,:), d2u_d2e(:,:), d3u_d3e(:,:)
        real(kind=8), allocatable :: Uz(:,:), dUz(:,:), d2Uz(:,:), d3Uz(:,:)
        real(kind=8), allocatable :: ERz(:,:), dERz(:,:), d2ERz(:,:), d3ERz(:,:)
        type(container_type), allocatable :: container(:), container2(:)
        real(kind=8), allocatable :: A(:,:), Chol(:,:), A2(:,:), Chol2(:,:), B_(:,:), C_(:,:), B_2(:,:), C_2(:,:)
               
        
        AR_type = self.config.AR_type
        Nz = self.config.Nz()
        Nb = self.config.Nb()
        Nsv = self.config.Nsv()
        Nuz = self.config.Nuz()
        dt_f = self.config.dt_f
        measures = ['times', 'B', 'U', 'ER']
        
        if (TRIM(self.config.prior_type) == '100path') then
            allocate(prior_data_obj(10))
            allocate(set_prior_obj(10))
            allocate(tag_list(2))
        else
            allocate(prior_data_obj(1))
            allocate(set_prior_obj(1))
            allocate(tag_list(1))
        end if
        
        call extract_realisations(self.config.prior_dir, self.config.prior_type, self.config.dt_smoothing, measures, prior_data_obj)
        
        
        ! # Applying unique function on array to get list of tags
        if (TRIM(self.config.prior_type) == '100path') then
            tag_list(1) = prior_data_obj(1).tag
            tag_list(2) = prior_data_obj(2).tag
        else
            tag_list(1) = prior_data_obj(1).tag
        end if
        
        !# Store covariance
        allocate(cov_prior.Z_Z(Nz,Nz))
        allocate(cov_prior.B_B(Nb,Nb))
        allocate(cov_prior.U_U(Nuz,Nuz))
        allocate(cov_prior.ER_ER(Nsv,Nsv))
        
        if (trim(AR_type) == 'AR3') then
            allocate(cov_prior.dZ_dZ(Nz,Nz))
            allocate(cov_prior.d2Z_d2Z(Nz,Nz))
        end if
        
        !# Store average
        !avg_prior = {}
        
        !# Init container
        allocate(container(size(tag_list)))
        !# If U ER not combined then a second container is needed
        if (self.config.combined_U_ER_forecast == 0) then
            allocate(container2(size(tag_list)))
        end if
        
        !# Loop over tag list 
        !# 100path MUST NOT BE THE FIRST IN TAG LIST THIS WILL RESULT IN AN ERROR
        !# BECAUSE PCA FIT (AND AVERAGES) MUST BE DONE ON A GEODYNAMO WITH LONG TIME SERIES (LIKE 71PATH)
        do j = 1, size(tag_list)
            
            t = tag_list(j)
            !# Init matrices to store iteratively if many time series in tag 
            !# (case for 100path which combines 100p and 71p priors) 
            !# Thus we can mix geodynamo priors
            if (ALLOCATED(Z)) deallocate(Z)
            if (ALLOCATED(dZ)) deallocate(dZ)
            if (ALLOCATED(d2Z)) deallocate(d2Z)
            if (ALLOCATED(d3Z)) deallocate(d3Z)
            allocate(Z(0,Nz), dZ(0,Nz), d2Z(0,Nz), d3Z(0,Nz))
            
            do i = 1, SIZE(prior_data_obj)
                if (trim(t) == prior_data_obj(i).tag) then
                    !# Set prior size 
                    call self.set_prior_size(prior_data_obj(i).U, self.config.Lu, 'U', set_prior_obj(i).U)
                    call self.set_prior_size(prior_data_obj(i).MF, self.config.Lb, 'MF', set_prior_obj(i).MF)
                    call self.set_prior_size(prior_data_obj(i).ER, self.config.Lsv, 'ER', set_prior_obj(i).ER)
                    
                    dt_sampling = prior_data_obj(i).dt_samp  !# time sampling
                    
                    !# times
                    allocate(set_prior_obj(i).times(SIZE(prior_data_obj(i).times, 1)))
                    set_prior_obj(i).times = prior_data_obj(i).times
                    
                    if (AR_type == "AR3") then
                        dt_prior = set_prior_obj(i).times(2)-set_prior_obj(i).times(1) !# geodynamo sampling
                        ! Compute blackman smoothing window length as the rounded ratio dt_sampling/dt_prior
                        Nt = NINT(dt_sampling / dt_prior)
                        ! # Compute blackman smoothing window total weight by summing all blackman window coeffs
                        blackman_w = sum_blackman(Nt)
                    else
                        ! No smoothing
                        Nt = 1
                        blackman_w = 1.0d0
                        ! # Perform a subsampling of U MF ER times according to self.config.dt_sampling
                        call sample_timed_quantity(prior_data_obj(i).times, set_prior_obj(i).U, self.config.dt_sampling, matrix_1, matrix_2)
                        deallocate(set_prior_obj(i).U)
                        allocate(set_prior_obj(i).U(SIZE(matrix_2, 1), SIZE(matrix_2, 2)))
                        set_prior_obj(i).U = matrix_2
                        
                        call sample_timed_quantity(prior_data_obj(i).times, set_prior_obj(i).MF, self.config.dt_sampling, matrix_1, matrix_2)
                        deallocate(set_prior_obj(i).MF)
                        allocate(set_prior_obj(i).MF(SIZE(matrix_2, 1), SIZE(matrix_2, 2)))
                        set_prior_obj(i).MF = matrix_2
                        
                        call sample_timed_quantity(prior_data_obj(i).times, set_prior_obj(i).ER, self.config.dt_sampling, matrix_1, matrix_2)
                        deallocate(set_prior_obj(i).ER)
                        allocate(set_prior_obj(i).ER(SIZE(matrix_2, 1), SIZE(matrix_2, 2)))
                        set_prior_obj(i).ER = matrix_2
                        
                        deallocate(set_prior_obj(i).times)
                        allocate(set_prior_obj(i).times(SIZE(matrix_1, 1)))
                        set_prior_obj(i).times = matrix_1
                        
                        dt_prior = set_prior_obj(i).times(2)-set_prior_obj(i).times(1) !# geodynamo sampling
                    end if
                    
                    if (TRIM(t) .ne. '100path') then  !#We consider average of 70path for 100path because longer time series
                        !# Compute U, B, ER averages
                        if (ALLOCATED(avg_prior.U)) deallocate(avg_prior.U)
                        if (ALLOCATED(avg_prior.MF)) deallocate(avg_prior.MF)
                        if (ALLOCATED(avg_prior.ER)) deallocate(avg_prior.ER)
                        if (ALLOCATED(avg_prior.times)) deallocate(avg_prior.times)
                        
                        allocate(avg_prior.U(1, SIZE(set_prior_obj(i).U, 2)))
                        allocate(avg_prior.MF(1, SIZE(set_prior_obj(i).MF, 2)))
                        allocate(avg_prior.ER(1, SIZE(set_prior_obj(i).ER, 2)))
                        
                        avg_prior.U = reshape(sum(set_prior_obj(i).U, dim=1) / SIZE(set_prior_obj(i).U, 1), [1, SIZE(avg_prior.U, 2)])
                        avg_prior.MF = reshape(sum(set_prior_obj(i).MF, dim=1) / SIZE(set_prior_obj(i).MF, 1), [1, SIZE(avg_prior.MF, 2)])
                        avg_prior.ER = reshape(sum(set_prior_obj(i).ER, dim=1) / SIZE(set_prior_obj(i).ER, 1), [1, SIZE(avg_prior.ER, 2)])

                    end if
                    
                    !# Center data
                    if (self.check_PCA()) then
                        if (TRIM(t) .ne. '100path') then ! # PCA fit over 71path because longer time series
                            if (.not. (ALLOCATED(self.pcaU_operator))) allocate(self.pcaU_operator)
                            call self.pcaU_operator.init_PCAOperator(self.config)
                            call self.pcaU_operator.fit(set_prior_obj(i).U)
                        end if
                        !# PCA transform of U
                        call self.pcaU_operator.transform(matrix_temp, set_prior_obj(i).U)
                        deallocate(set_prior_obj(i).U)
                        allocate(set_prior_obj(i).U(SIZE(matrix_temp, 1), SIZE(matrix_temp, 2)))
                        set_prior_obj(i).U = matrix_temp
                    else
                        !# Remove mean U
                        do concurrent (k = 1: SIZE(set_prior_obj(i).U, 1))
                             set_prior_obj(i).U(k,:) = set_prior_obj(i).U(k,:) - avg_prior.U(1,:)
                        end do
                    end if
                    
                    !# Remove mean ER
                    do concurrent (k = 1: SIZE(set_prior_obj(i).ER, 1))
                        set_prior_obj(i).ER(k,:) = set_prior_obj(i).ER(k,:) - avg_prior.ER(1,:)
                    end do
                    
                    if (trim(AR_type) /= "diag") then
                        call prep_AR_matrix(set_prior_obj(i).U, dt_prior, Nt, u, du, d2u, d3u)                        
                        call prep_AR_matrix(set_prior_obj(i).ER, dt_prior, Nt, e, de, d2e, d3e)
                        call prep_AR_matrix(set_prior_obj(i).MF, dt_prior, Nt, b, db, d2b, d3b)
                        
                        call self.U_ER_to_Z2(u, e, u_e)
                        call self.U_ER_to_Z2(du, de, du_de)
                        call self.U_ER_to_Z2(d2u, d2e, d2u_d2e)
                        call self.U_ER_to_Z2(d3u, d3e, d3u_d3e)
                        
                        call np_concatenate0(Z, u_e, matrix_temp)
                        deallocate(Z)
                        allocate(Z(SIZE(matrix_temp,1), SIZE(matrix_temp,2)))
                        Z = matrix_temp   
                        
                        call np_concatenate0(dZ, du_de, matrix_temp)
                        deallocate(dZ)
                        allocate(dZ(SIZE(matrix_temp,1), SIZE(matrix_temp,2)))
                        dZ = matrix_temp
                        
                        call np_concatenate0(d2Z, d2u_d2e, matrix_temp)
                        deallocate(d2Z)
                        allocate(d2Z(SIZE(matrix_temp,1), SIZE(matrix_temp,2)))
                        d2Z = matrix_temp
                        
                        call np_concatenate0(d3Z, d3u_d3e, matrix_temp)
                        deallocate(d3Z)
                        allocate(d3Z(SIZE(matrix_temp,1), SIZE(matrix_temp,2)))
                        d3Z = matrix_temp
                    else
                        if (ALLOCATED(u)) deallocate(u)
                        if (ALLOCATED(e)) deallocate(e)
                        if (ALLOCATED(b)) deallocate(b)
                        allocate(u(SIZE(set_prior_obj(i).U,1), SIZE(set_prior_obj(i).U,2)))
                        allocate(e(SIZE(set_prior_obj(i).ER,1), SIZE(set_prior_obj(i).ER,2)))
                        allocate(b(SIZE(set_prior_obj(i).MF,1), SIZE(set_prior_obj(i).MF,2)))
                        u = set_prior_obj(i).U
                        e = set_prior_obj(i).ER
                        b = set_prior_obj(i).MF                        
                    end if                    
                    
                    if (TRIM(t) /= "100path") then!#We consider covariance of 70path for 100path because longer time series
                        !# Compute U, B, ER covariance matrices
                        call cov(u / blackman_w, cov_prior.U_U)
                        call cov(b / blackman_w, cov_prior.B_B)
                        call cov(e / blackman_w, cov_prior.ER_ER)
                        
                        if (trim(AR_type) == "diag") then
                            ! # Diag is independant for forecast so diag block U,U and ER,ER
                            call block_diag(cov_prior.U_U, cov_prior.ER_ER, cov_prior.Z_Z)
                            
                        else
                            !# If U and ER dependant for forecast.
                            if (self.config.combined_U_ER_forecast == 1) then
                                call cov(Z / blackman_w, cov_prior.Z_Z)
                                call cov(d2Z / blackman_w, cov_prior.d2Z_d2Z)
                                call cov(dZ / blackman_w, cov_prior.dZ_dZ)
                            else
                                !# Concatenate
                                call block_diag(cov_prior.U_U, cov_prior.ER_ER, cov_prior.Z_Z)
                                
                                call cov(du / blackman_w, matrix_temp)
                                call cov(de / blackman_w, matrix_temp2)
                                call block_diag(matrix_temp, matrix_temp2, cov_prior.dZ_dZ)
                                
                                call cov(d2u / blackman_w, matrix_temp)
                                call cov(d2e / blackman_w, matrix_temp2)
                                call block_diag(matrix_temp, matrix_temp2, cov_prior.d2Z_d2Z)
                            end if                            
                        end if                       
                    end if                
                end if                
            end do
            
            if (ALLOCATED(Uz)) deallocate(Uz)
            if (ALLOCATED(dUz)) deallocate(dUz)
            if (ALLOCATED(d2Uz)) deallocate(d2Uz)
            if (ALLOCATED(d3Uz)) deallocate(d3Uz)
            if (ALLOCATED(ERz)) deallocate(ERz)
            if (ALLOCATED(dERz)) deallocate(dERz)
            if (ALLOCATED(d2ERz)) deallocate(d2ERz)
            if (ALLOCATED(d3ERz)) deallocate(d3ERz)
            
            allocate(Uz(SIZE(Z, 1), Nuz), dUz(SIZE(Z, 1), Nuz), d2Uz(SIZE(Z, 1), Nuz), d3Uz(SIZE(Z, 1), Nuz))
            allocate(ERz(SIZE(Z, 1), SIZE(Z, 2)-Nuz), dERz(SIZE(Z, 1), SIZE(Z, 2)-Nuz), d2ERz(SIZE(Z, 1), SIZE(Z, 2)-Nuz), d3ERz(SIZE(Z, 1), SIZE(Z, 2)-Nuz))
            
            if (trim(AR_type) /= "diag") then
                !# If U and ER independant for forecast
                if (self.config.combined_U_ER_forecast == 0)then 
                    ! # Split U ER parts of Z
                    Uz = Z(:,1:Nuz)
                    dUz = dZ(:,1:Nuz)
                    d2Uz = d2Z(:,1:Nuz)
                    d3Uz = d3Z(:,1:Nuz)
                    ERz = Z(:,Nuz+1:SIZE(Z, 2))
                    dERz = dZ(:,Nuz+1:SIZE(Z, 2))
                    d2ERz = d2Z(:,Nuz+1:SIZE(Z, 2))
                    d3ERz = d3Z(:,Nuz+1:SIZE(Z, 2))
                end if
                
                if (trim(AR_type) == "AR1") then
                    !# If U and ER dependant for forecast
                    if (self.config.combined_U_ER_forecast == 1) then
                        ! # Append container for Z
                        call init_container(container(j), TRANSPOSE(Z), TRANSPOSE(dZ)*dt_prior, dt_prior, Nt)
                    else
                        ! # Append container for U
                        call init_container(container(j), TRANSPOSE(Uz), TRANSPOSE(dUz)*dt_prior, dt_prior, Nt)
                        ! # Append container 2 for ER
                        call init_container(container2(j), TRANSPOSE(ERz), TRANSPOSE(dERz)*dt_prior, dt_prior, Nt)
                    end if
                else if (trim(AR_type) == "AR3") then
                    ! # If U and ER dependant for forecast
                    if (self.config.combined_U_ER_forecast == 1) then
                        ! # Append container for Z
                        call np_concatenate0(TRANSPOSE(Z), TRANSPOSE(dZ), matrix_temp)
                        call np_concatenate0(matrix_temp, TRANSPOSE(d2Z), matrix_3) ! matrix_1= np.concatenate((Z.T,dZ.T,d2Z.T),axis=0)
                        call np_concatenate0(TRANSPOSE(dZ), TRANSPOSE(d2Z), matrix_temp)
                        call np_concatenate0(matrix_temp, TRANSPOSE(d3Z), matrix_2) ! matrix_2= np.concatenate((dZ.T,d2Z.T,d3Z.T),axis=0)
                        call init_container(container(j), matrix_3, matrix_2*dt_prior, dt_prior, Nt)
                    else
                        ! # Append container for U
                        call np_concatenate0(TRANSPOSE(Uz), TRANSPOSE(dUz), matrix_temp)
                        call np_concatenate0(matrix_temp, TRANSPOSE(d2Uz), matrix_3) ! matrix_1= np.concatenate((Uz.T,dUz.T,d2Uz.T),axis=0)
                        call np_concatenate0(TRANSPOSE(dUz), TRANSPOSE(d2Uz), matrix_temp)
                        call np_concatenate0(matrix_temp, TRANSPOSE(d3Uz), matrix_2) ! matrix_2= np.concatenate((dUz.T,d2Uz.T,d3Uz.T),axis=0)
                        call init_container(container(j), matrix_3, matrix_2*dt_prior, dt_prior, Nt)
                        
                        ! # Append container 2 for ER
                        call np_concatenate0(TRANSPOSE(ERz), TRANSPOSE(dERz), matrix_temp)
                        call np_concatenate0(matrix_temp, TRANSPOSE(d2ERz), matrix_3) ! matrix_1= np.concatenate((ERz.T,dERz.T,d2ERz.T),axis=0)
                        call np_concatenate0(TRANSPOSE(dERz), TRANSPOSE(d2ERz), matrix_temp)
                        call np_concatenate0(matrix_temp, TRANSPOSE(d3ERz), matrix_2) ! matrix_2= np.concatenate((dERz.T,d2ERz.T,d3ERz.T),axis=0) 
                        call init_container(container2(j), matrix_3, matrix_2*dt_prior, dt_prior, Nt)
                    end if                   
                end if                
            end if            
        end do
        
        if (trim(AR_type) == "AR1") then
            ! # compute A, Chol of U or Z
            call compute_AR_coefs_avg(container, AR_type, A, B_, C_, Chol)
            if (self.config.combined_U_ER_forecast == 0)then
                !# compute A, Chol of ER
                call compute_AR_coefs_avg(container2, AR_type, A2, B_, C_, Chol2)
                !# diag block 
                if (ALLOCATED(matrix_temp)) deallocate(matrix_temp)
                if (ALLOCATED(matrix_temp2)) deallocate(matrix_temp2)
                allocate(matrix_temp, source=A)
                allocate(matrix_temp2, source=Chol)
                call block_diag(matrix_temp, A2, A)
                call block_diag(matrix_temp2, Chol2, Chol)
                !# compute A, Chol for forecast
                call compute_AR1_coefs_forecast(A, Chol, dt_f, Nz, cov_prior.A, cov_prior.Chol)
            end if
        else if (trim(AR_type) == "diag") then
            !# compute A, Chol
            call compute_diag_AR1_coefs(cov_prior.U_U, cov_prior.ER_ER, self.config.TauU, self.config.TauE, A, Chol)
            call compute_AR1_coefs_forecast(A, Chol, dt_f, Nz, cov_prior.A, cov_prior.Chol)
        else if (trim(AR_type) == "AR3") then
            !# compute A, B, C, Chol of U or Z
            call compute_AR_coefs_avg(container, AR_type, A, B_, C_, Chol)
            if (self.config.combined_U_ER_forecast == 0) then
                !# compute A, B, C, Chol of ER
                call compute_AR_coefs_avg(container2, AR_type, A2, B_2, C_2, Chol2)
                !# diag block 
                if (ALLOCATED(matrix_temp)) deallocate(matrix_temp)
                if (ALLOCATED(matrix_temp2)) deallocate(matrix_temp2)
                if (ALLOCATED(matrix_2)) deallocate(matrix_2)
                if (ALLOCATED(matrix_3)) deallocate(matrix_3)
                allocate(matrix_temp, source=A)
                allocate(matrix_temp2, source=B_)
                allocate(matrix_2, source=C_)
                allocate(matrix_3, source=Chol)
                call block_diag(matrix_temp, A2, A)
                call block_diag(matrix_temp2, B_2, B_)
                call block_diag(matrix_2, C_2, C_)
                call block_diag(matrix_3, Chol2, Chol)
            end if
            !# compute A, B, C, Chol forecast
            call compute_AR3_coefs_forecast(A, B_, C_, Chol, dt_f, Nz, cov_prior.A, cov_prior.B, cov_prior.C, cov_prior.Chol)
        end if
        write (10,'(A)') '================================================'
        write (10,'(A)') 'Reading and computation of priors/covariances of AugkfAlgo OK !'        
    end subroutine    
!==========================================================================================================================
    
!==========================================================================================================================
    subroutine set_prior_size(self, X, asked_max_degree, measure, result_)
        class(AugkfAlgo), intent(in) :: self
        real(kind=8), intent(in) :: X(:,:)
        integer, intent(in) :: asked_max_degree
        integer :: asked_coeffs, read_Nu
        character(len=*) :: measure
        real(kind=8), intent(out), allocatable :: result_(:,:)
        
        if (ALLOCATED(result_)) deallocate(result_)
        asked_coeffs = asked_max_degree * (asked_max_degree + 2)
        if (TRIM(measure) == 'MF' .OR. TRIM(measure) == 'SV' .OR. TRIM(measure) == 'ER') then
            call self.check_coef_size(X, asked_max_degree, asked_coeffs)
            allocate(result_(SIZE(X, 1), asked_coeffs))
            result_ = X(:, 1:asked_coeffs)
        else
            call self.check_coef_size(X, asked_max_degree, asked_coeffs)
            read_Nu = SIZE(X, 2) / 2
            allocate(result_(SIZE(X, 1), asked_coeffs*2))
            result_(:,1:asked_coeffs) = X(:, 1:asked_coeffs)
            result_(:,asked_coeffs+1:asked_coeffs*2) = X(:, read_Nu+1:read_Nu+asked_coeffs)
        endif    
    end subroutine    
    !function set_prior_size(self, X, asked_max_degree, measure) result(result_)
    !    class(AugkfAlgo) :: self
    !    real(kind=8) :: X(:,:)
    !    integer :: asked_max_degree, asked_coeffs, read_Nu
    !    character(len=*) :: measure
    !    real(kind=8), allocatable :: result_(:,:)
    !    
    !    asked_coeffs = asked_max_degree * (asked_max_degree + 2)
    !    if (TRIM(measure) == 'MF' .OR. TRIM(measure) == 'SV' .OR. TRIM(measure) == 'ER') then
    !        call self.check_coef_size(X, asked_max_degree, asked_coeffs)
    !        allocate(result_(SIZE(X, 1), asked_coeffs))
    !        result_ = X(:, 1:asked_coeffs)
    !    else
    !        call self.check_coef_size(X, asked_max_degree, asked_coeffs)
    !        read_Nu = SIZE(X, 2) / 2
    !        allocate(result_(SIZE(X, 1), asked_coeffs*2))
    !        result_(:,1:asked_coeffs) = X(:, 1:asked_coeffs)
    !        result_(:,asked_coeffs+1:asked_coeffs*2) = X(:, read_Nu+1:read_Nu+asked_coeffs)
    !    endif    
    !end function    
!==========================================================================================================================
    
!==========================================================================================================================
    subroutine check_coef_size(self, X, asked_max_degree, asked_coeffs)
    !*****************************************************************************************************************
    !"""
    !# Take only data corresponding to the asked coefs (and check that they are not bigger than the data...)
    !"""
    !*****************************************************************************************************************
        class(AugkfAlgo), intent(in) :: self
        real(kind=8), intent(in) :: X(:,:)
        integer, intent(in) :: asked_max_degree, asked_coeffs
        
        if (SIZE(X, 2) < asked_coeffs) then
            write(10,'(A, i, A, A, A)') 'Asked max degree', asked_max_degree ,'is too big to be handled by the prior data in ', TRIM(self.config.prior_dir), '. Please retry with a lower max degree.'
            stop
        end if
    end subroutine  
!==========================================================================================================================
    
!==========================================================================================================================
    subroutine U_ER_to_Z2(self, U, ER, result_)
    !*****************************************************************************************************************
    !"""
    !Compute Z from U ER
    !
    !:param U: U state
    !:type U: numpy array, dim Ncoef (1D) or Ntimes x Ncoef (2D) or Nreal x Ntimes x Ncoef (3D)
    !:return: U and ER states [U  ER]
    !:rtype: 1D or 2D or 3D Z states
    !"""
    !*****************************************************************************************************************
    class(AugkfAlgo), intent(in) :: self
    real(kind=8), intent(in) :: U(:,:), ER(:,:)
    real(kind=8), allocatable, intent(out) :: result_(:,:)
    
    allocate(result_(SIZE(U, 1), SIZE(U, 2) + SIZE(ER, 2)))
    result_(:,1:SIZE(U, 2)) = U
    result_(:,SIZE(U, 2)+1:SIZE(U, 2)+SIZE(ER, 2)) = ER
    end subroutine    
!==========================================================================================================================
    
end module
    
    