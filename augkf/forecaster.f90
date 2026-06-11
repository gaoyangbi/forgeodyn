module forecaster
    use mpi
    use computer
    use corestate
    use config
    use common
    use pca
    implicit none
    
    type, extends(GenericComputer),public :: AugkfForecasterAR1
        logical :: Cholesky_AR_check
    contains
        procedure :: init_AugkfForecasterAR
        procedure :: forecast_step_AR1
        procedure :: parallel_forecast_step_AR1, forecast_Z
    end type AugkfForecasterAR1
    
    type, extends(AugkfForecasterAR1),public :: AugkfForecasterAR3
        
    contains
        !procedure :: forecast_step
    end type AugkfForecasterAR3
    
contains    
!==========================================================================================================================
    subroutine init_AugkfForecasterAR(self, config, legendre_polys)
    !*****************************************************************************************************************
    !"""
    !Class that implements the forecasts using AugKF (Augmented state Kalman Filter) algorithm with DIFF treated as a contribution to ER.
    !"""
    !*****************************************************************************************************************
        class(AugkfForecasterAR1), intent(inout) :: self
        class(ComputationConfig), intent(in) :: config
        class(legendre_polys_type), intent(in) :: legendre_polys
        call self.init_GenericComputer(config, legendre_polys, 0, 0)
        
        !# Bool to deactivate checks on AR processes
        self.Cholesky_AR_check = .false.
    end subroutine init_AugkfForecasterAR
!==========================================================================================================================
    
!==========================================================================================================================
    subroutine forecast_step_AR1(self, input_core_state, algo_nb_realisations, algo_config, algo_pcaU_operator, algo_avg_prior, algo_cov_prior, Z_AR, seed, i_real, i_t, next_core_state)
    !*****************************************************************************************************************
    !"""
    !Forecasts the input_core_state using AR processes for Z, computation of SV and Euler scheme for B.
    !
    !:param input_core_state: core_state of a single realisation at a single date 
    !:type input_core_state: corestates.CoreState
    !:param Z_AR: forecast state AR
    !:type Z_AR: np.array 1D if AR1 (Ncoef) or 2D if AR3 (3 x Ncoef)
    !:param seed: random seed
    !:type seed: int
    !:param i_real: model realisation index
    !:type i_real: int
    !:param i_t: time index
    !:type i_t: int
    !:return: CoreState containing the result from the forecast
    !:rtype: corestates.CoreState
    !"""
    !*****************************************************************************************************************
        class(AugkfForecasterAR1), intent(in) :: self
        class(CoreState_type), intent(in) :: input_core_state
        integer, intent(in) :: algo_nb_realisations
        class(ComputationConfig), intent(in) :: algo_config
        class(NormedPCAOperator), intent(in) :: algo_pcaU_operator
        class(set_prior_type), intent(in) :: algo_avg_prior
        class(cov_prior_type), intent(in) :: algo_cov_prior
        real(kind=8), intent(in) :: Z_AR(:)
        integer, intent(in) :: seed, i_real, i_t
        class(CoreState_type), allocatable, intent(out) :: next_core_state
        integer :: N, rstate
        real(kind=8), allocatable :: Z_AR1_forecast(:), Ab(:,:)
        type(input_core_state_type) :: CoreState_temp
        
        !# copy input core state
        allocate(next_core_state, source=input_core_state)
        
        !# set random state
        N = algo_nb_realisations
        rstate = seed + i_real + N * i_t
        
        !# Compute Z(t+1)
        call self.forecast_Z(Z_AR, algo_cov_prior, rstate, Z_AR1_forecast)
        next_core_state.measures_(5).measure_data = RESHAPE(Z_AR1_forecast, SHAPE(next_core_state.measures_(5).measure_data))
        call Z_to_U_ER1(algo_config, algo_avg_prior, algo_pcaU_operator, Z_AR1_forecast, next_core_state.measures_(2).measure_data(1,1,:), next_core_state.measures_(4).measure_data(1,1,:))
        
        !# Compute A(b)
        CoreState_temp.Lsv = next_core_state.cs_Lsv()
        CoreState_temp.Lu = next_core_state.cs_Lu()
        CoreState_temp.Lb = next_core_state.cs_Lb()
        CoreState_temp.Nsv = next_core_state.cs_Nsv()
        CoreState_temp.Nu2 = next_core_state.cs_Nu2()
        CoreState_temp.Nb = next_core_state.cs_Nb()
        allocate(CoreState_temp.B, source = next_core_state.measures_(1).measure_data(1, 1, :))
        call self.compute_Ab(CoreState_temp, Ab)
        
        !# Compute SV(t+1) = A(b)U(t+1) + E(t+1)
        next_core_state.measures_(3).measure_data(1,1,:) = MATMUL(Ab, next_core_state.measures_(2).measure_data(1,1,:)) + next_core_state.measures_(4).measure_data(1,1,:)
        
        !# Compute B(t+1) = B(t) + dt*SV(t+1) (Euler scheme as SV=dB/dt)
        next_core_state.measures_(1).measure_data(1,1,:) = input_core_state.measures_(1).measure_data(1,1,:) + algo_config.dt_f * next_core_state.measures_(3).measure_data(1,1,:)
        
    end subroutine forecast_step_AR1
!==========================================================================================================================
    
!==========================================================================================================================
    subroutine forecast_Z(self, Z_AR1, alogo_cov_prior, rstate, Z_AR1_forecast)
    !*****************************************************************************************************************
    !"""
    !Forecast Z state with AR-1 process.
    !
    !:param Z_AR1: AR-1 forecast state
    !:type Z_AR1: np.array (Ncoef)
    !:param rstate: Random state to use for the AR-1 process
    !:type rstate: np.random.RandomState
    !:return: forecasted Z state
    !:rtype: np.array(Ncoef)
    !"""
    !*****************************************************************************************************************
        class(AugkfForecasterAR1), intent(in) :: self
        real(kind=8), intent(in) :: Z_AR1(:)
        class(cov_prior_type), intent(in) :: alogo_cov_prior
        integer, intent(in) :: rstate
        real(kind=8), intent(out), allocatable :: Z_AR1_forecast(:)
        
        !# AR1 process for Z
        call ar1_process(Z_AR1, alogo_cov_prior.A, alogo_cov_prior.Chol, rstate, .True., Z_AR1_forecast)        
    end subroutine forecast_Z
!==========================================================================================================================

!==========================================================================================================================
    subroutine parallel_forecast_step_AR1(self, algo_config, algo_nb_realisations, algo_attributed_models, algo_pcaU_operator, algo_avg_prior, algo_cov_prior, input_core_state, seed, i_t, forecast_at_t)
    !*****************************************************************************************************************
    !"""
    !parallelize the AR1 forecast step
    !
    !:param input_core_state: input_core_state at time t
    !:type input_core_states: corestates.Corestate
    !:param seed: random seed
    !:type seed: int
    !:param i_t: time index
    !:type i_t: int
    !:return: CoreState containing the result from the forecast
    !:rtype: corestates.CoreState
    !"""
    !*****************************************************************************************************************
        class(AugkfForecasterAR1), intent(inout) :: self
        class(ComputationConfig), intent(in) :: algo_config
        integer, intent(in) :: algo_nb_realisations
        integer, intent(in) :: algo_attributed_models(:)
        class(NormedPCAOperator), intent(in) :: algo_pcaU_operator
        class(set_prior_type), intent(in) :: algo_avg_prior
        class(cov_prior_type), intent(in) :: algo_cov_prior
        class(CoreState_type), intent(in) :: input_core_state
        class(CoreState_type), allocatable, intent(out) :: forecast_at_t
        integer, intent(in) :: seed, i_t
        real(kind=8) :: t
        type(CoreState_type) :: core_state_slice
        integer :: i_idx, j, i, num, comm, rank, ierr
        class(CoreState_type), allocatable :: next_core_state
        character(len=1000) :: to_print
        
        
        t = algo_config.t_forecasts(i_t+1)
        
        !# copy input core state
        allocate(forecast_at_t, source = input_core_state)
        !# set all measures to 0
        do i = 1, size(forecast_at_t.measures_)
            forecast_at_t.measures_(i).measure_data = 0.0d0
        end do
        
        !# Each process computes its attributed models
        allocate(core_state_slice.measures_, source=input_core_state.measures_)
        allocate(core_state_slice.max_degrees_, source=input_core_state.max_degrees_)
        
        !# MPI rank
        comm = MPI_COMM_WORLD
        call MPI_Comm_rank(comm, rank, ierr)
        
        do i_idx = 1, size(algo_attributed_models)
            ! slice the core state for the model i
            do j = 1, size(input_core_state.measures_)
                deallocate(core_state_slice.measures_(j).measure_data)
                allocate(core_state_slice.measures_(j).measure_data(1, 1, SIZE(input_core_state.measures_(j).measure_data, 3)))
                core_state_slice.measures_(j).measure_data = input_core_state.measures_(j).measure_data(i_idx:i_idx, :, :)
            end do
            num = SIZE(core_state_slice.measures_(5).measure_data)
            call self.forecast_step_AR1(core_state_slice, algo_nb_realisations, algo_config, algo_pcaU_operator, algo_avg_prior, algo_cov_prior, RESHAPE(core_state_slice.measures_(5).measure_data, [num]), seed, i_idx, i_t, next_core_state)
            
            do j = 1, size(input_core_state.measures_)
                forecast_at_t.measures_(j).measure_data(i_idx:i_idx,1:1,:) = next_core_state.measures_(j).measure_data
            end do
            call coef_print(next_core_state, 1, to_print)
            write(10, '(a,i2,a,a)') 'Process-'//'rank:', rank, '  ', trim(to_print)
            write(*, '(a,i2,a,a)') 'Process-'//'rank:', rank, '  ', trim(to_print)
        end do        
    end subroutine parallel_forecast_step_AR1
!==========================================================================================================================
    
end module