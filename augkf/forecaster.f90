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
        procedure :: forecast_Z_AR3
        procedure :: forecast_step_AR3
        procedure :: update_Z_AR3
        procedure :: parallel_forecast_step_AR3
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
    subroutine forecast_step_AR1(self, input_core_state, algo_nb_realisations, algo_config, algo_pcaU_operator, algo_avg_prior, algo_cov_prior, Z_AR, global_seed, global_i_real, i_t, next_core_state)
    !*****************************************************************************************************************
    !"""
    !Forecasts the input_core_state using AR processes for Z, computation of SV and Euler scheme for B.
    !
    !:param input_core_state: core_state of a single realisation at a single date 
    !:type input_core_state: corestates.CoreState
    !:param Z_AR: forecast state AR
    !:type Z_AR: np.array 1D if AR1 (Ncoef) or 2D if AR3 (3 x Ncoef)
    !:param global_seed: global random seed
    !:type global_seed: int
    !:param global_i_real: global model realisation index
    !:type global_i_real: int
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
        integer, intent(in) :: global_seed, global_i_real, i_t
        class(CoreState_type), allocatable, intent(out) :: next_core_state
        real(kind=8), allocatable :: Z_AR1_forecast(:), Ab(:,:)
        type(input_core_state_type) :: CoreState_temp
        
        !# copy input core state
        allocate(next_core_state, source=input_core_state)
        
        !# Compute Z(t+1)
        call self.forecast_Z(Z_AR, algo_cov_prior, global_seed, global_i_real, i_t, Z_AR1_forecast)
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
    subroutine forecast_step_AR3(self, input_core_state, algo_nb_realisations, algo_config, algo_pcaU_operator, algo_avg_prior, algo_cov_prior, Z_AR3, global_seed, global_i_real, i_t, next_core_state)
    !*****************************************************************************************************************
    !"""
    !Forecasts the input_core_state using an AR-3 process for Z,
    !computation of SV and Euler scheme for B.
    !
    !:param input_core_state: core state of a single realisation at a single date
    !:type input_core_state: corestates.CoreState
    !:param Z_AR3: AR-3 forecast state
    !:type Z_AR3: np.array (3 x Ncoef)
    !:param global_seed: global random seed
    !:type global_seed: int
    !:param global_i_real: global model realisation index
    !:type global_i_real: int
    !:param i_t: time index
    !:type i_t: int
    !:return: CoreState containing the result from the forecast
    !:rtype: corestates.CoreState
    !"""
    !*****************************************************************************************************************
        class(AugkfForecasterAR3), intent(in) :: self
        class(CoreState_type), intent(in) :: input_core_state
        integer, intent(in) :: algo_nb_realisations
        class(ComputationConfig), intent(in) :: algo_config
        class(NormedPCAOperator), intent(in) :: algo_pcaU_operator
        class(set_prior_type), intent(in) :: algo_avg_prior
        class(cov_prior_type), intent(in) :: algo_cov_prior
        real(kind=8), intent(in) :: Z_AR3(:,:)
        integer, intent(in) :: global_seed, global_i_real, i_t
        class(CoreState_type), allocatable, intent(out) :: next_core_state

        real(kind=8), allocatable :: Z_AR3_forecast(:), Ab(:,:)
        type(input_core_state_type) :: CoreState_temp

        !# Copy input core state
        allocate(next_core_state, source=input_core_state)

        !# Compute Z(t+1)
        call self.forecast_Z_AR3( &
            Z_AR3,            &
            algo_cov_prior,   &
            global_seed,      &
            global_i_real,    &
            i_t,              &
            Z_AR3_forecast    &
        )

        next_core_state.measures_(5).measure_data = &
            RESHAPE(Z_AR3_forecast, SHAPE(next_core_state.measures_(5).measure_data))

        !# Compute U(t+1) and ER(t+1) from Z(t+1)
        call Z_to_U_ER1( &
            algo_config,                                      &
            algo_avg_prior,                                   &
            algo_pcaU_operator,                               &
            Z_AR3_forecast,                                   &
            next_core_state.measures_(2).measure_data(1,1,:), &
            next_core_state.measures_(4).measure_data(1,1,:)  &
        )

        !# Compute A(b)
        CoreState_temp.Lsv = next_core_state.cs_Lsv()
        CoreState_temp.Lu = next_core_state.cs_Lu()
        CoreState_temp.Lb = next_core_state.cs_Lb()
        CoreState_temp.Nsv = next_core_state.cs_Nsv()
        CoreState_temp.Nu2 = next_core_state.cs_Nu2()
        CoreState_temp.Nb = next_core_state.cs_Nb()
        allocate( &
            CoreState_temp.B, &
            source=next_core_state.measures_(1).measure_data(1,1,:) &
        )
        call self.compute_Ab(CoreState_temp, Ab)

        !# Compute SV(t+1) = A(b)U(t+1) + ER(t+1)
        next_core_state.measures_(3).measure_data(1,1,:) = &
            MATMUL(Ab, next_core_state.measures_(2).measure_data(1,1,:)) &
            + next_core_state.measures_(4).measure_data(1,1,:)

        !# Compute B(t+1) = B(t) + dt*SV(t+1)
        next_core_state.measures_(1).measure_data(1,1,:) = &
            input_core_state.measures_(1).measure_data(1,1,:) &
            + algo_config.dt_f * next_core_state.measures_(3).measure_data(1,1,:)

    end subroutine forecast_step_AR3
!==========================================================================================================================

!==========================================================================================================================
    subroutine forecast_Z(self, Z_AR1, alogo_cov_prior, global_seed, global_i_real, i_t, Z_AR1_forecast)
    !*****************************************************************************************************************
    !"""
    !Forecast Z state with AR-1 process.
    !
    !:param Z_AR1: AR-1 forecast state
    !:type Z_AR1: np.array (Ncoef)
    !:param global_seed: global random seed
    !:type global_seed: int
    !:param global_i_real: global model realisation index
    !:type global_i_real: int
    !:param i_t: time index
    !:type i_t: int
    !:return: forecasted Z state
    !:rtype: np.array(Ncoef)
    !"""
    !*****************************************************************************************************************
        class(AugkfForecasterAR1), intent(in) :: self
        real(kind=8), intent(in) :: Z_AR1(:)
        class(cov_prior_type), intent(in) :: alogo_cov_prior
        integer, intent(in) :: global_seed, global_i_real, i_t
        real(kind=8), intent(out), allocatable :: Z_AR1_forecast(:)
        
        !# AR1 process for Z
        call ar1_process(Z_AR1, alogo_cov_prior.A, alogo_cov_prior.Chol, global_seed, global_i_real, i_t, .True., Z_AR1_forecast)
    end subroutine forecast_Z
!==========================================================================================================================

!==========================================================================================================================
    subroutine forecast_Z_AR3(self, Z_AR3, algo_cov_prior, global_seed, global_i_real, i_t, Z_AR3_forecast)
    !*****************************************************************************************************************
    !"""
    !Forecast Z state with AR-3 process.
    !
    !:param Z_AR3: AR-3 forecast state
    !:type Z_AR3: np.array (3 x Ncoef)
    !:param global_seed: global random seed
    !:type global_seed: int
    !:param global_i_real: global model realisation index
    !:type global_i_real: int
    !:param i_t: time index
    !:type i_t: int
    !:return: forecasted Z state
    !:rtype: np.array(Ncoef)
    !"""
    !*****************************************************************************************************************
        class(AugkfForecasterAR3), intent(in) :: self
        real(kind=8), intent(in) :: Z_AR3(:,:)
        class(cov_prior_type), intent(in) :: algo_cov_prior
        integer, intent(in) :: global_seed, global_i_real, i_t
        real(kind=8), allocatable, intent(out) :: Z_AR3_forecast(:)

        !# AR3 process for Z
        call ar3_process( &
            Z_AR3,               &
            algo_cov_prior.A,    &
            algo_cov_prior.B,    &
            algo_cov_prior.C,    &
            algo_cov_prior.Chol, &
            global_seed,         &
            global_i_real,       &
            i_t,                 &
            .True.,              &
            Z_AR3_forecast       &
        )
    end subroutine forecast_Z_AR3
!==========================================================================================================================

!==========================================================================================================================
    subroutine update_Z_AR3(self, Z, Z_AR3)
    !*****************************************************************************************************************
    !"""
    !Update the three-state forecast history used by the AR-3 process.
    !
    !The history is shifted from:
    ![Z(n-2), Z(n-1), Z(n)]
    !to:
    ![Z(n-1), Z(n), Z(n+1)]
    !
    !:param Z: newly forecasted Z states
    !:type Z: np.array (Nreal x Ncoef)
    !:param Z_AR3: AR-3 forecast history
    !:type Z_AR3: np.array (Nreal x 3 x Ncoef)
    !"""
    !*****************************************************************************************************************
        class(AugkfForecasterAR3), intent(in) :: self
        real(kind=8), intent(in) :: Z(:,:)
        real(kind=8), intent(inout) :: Z_AR3(:,:,:)

        if (SIZE(Z_AR3, 1) .ne. SIZE(Z, 1)) then
            write (10,'(A)') &
                'Z and Z_AR3 must contain the same number of realisations.'
            write (*,'(A)') &
                'Z and Z_AR3 must contain the same number of realisations.'
            stop
        end if

        if (SIZE(Z_AR3, 2) .ne. 3) then
            write (10,'(A, i4, A)') &
                'Second dimension of Z_AR3 should be 3. Got ', &
                SIZE(Z_AR3, 2), ' instead.'
            write (*,'(A, i4, A)') &
                'Second dimension of Z_AR3 should be 3. Got ', &
                SIZE(Z_AR3, 2), ' instead.'
            stop
        end if

        if (SIZE(Z_AR3, 3) .ne. SIZE(Z, 2)) then
            write (10,'(A)') &
                'Z and Z_AR3 must contain the same number of coefficients.'
            write (*,'(A)') &
                'Z and Z_AR3 must contain the same number of coefficients.'
            stop
        end if

        !# Keep the three most recent Z states
        Z_AR3(:,1,:) = Z_AR3(:,2,:)
        Z_AR3(:,2,:) = Z_AR3(:,3,:)
        Z_AR3(:,3,:) = Z

    end subroutine update_Z_AR3
!==========================================================================================================================

!==========================================================================================================================
    subroutine parallel_forecast_step_AR1(self, algo_config, algo_nb_realisations, algo_attributed_models, algo_pcaU_operator, algo_avg_prior, algo_cov_prior, input_core_state, global_seed, i_t, forecast_at_t)
    !*****************************************************************************************************************
    !"""
    !parallelize the AR1 forecast step
    !
    !:param input_core_state: input_core_state at time t
    !:type input_core_states: corestates.Corestate
    !:param global_seed: global random seed
    !:type global_seed: int
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
        integer, intent(in) :: global_seed, i_t
        real(kind=8) :: t
        type(CoreState_type) :: core_state_slice
        integer :: i_idx, global_i_real, j, i, num, comm, rank, ierr
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
            global_i_real = algo_attributed_models(i_idx)
            do j = 1, size(input_core_state.measures_)
                deallocate(core_state_slice.measures_(j).measure_data)
                allocate(core_state_slice.measures_(j).measure_data(1, 1, SIZE(input_core_state.measures_(j).measure_data, 3)))
                core_state_slice.measures_(j).measure_data = input_core_state.measures_(j).measure_data(i_idx:i_idx, :, :)
            end do
            num = SIZE(core_state_slice.measures_(5).measure_data)
            call self.forecast_step_AR1(core_state_slice, algo_nb_realisations, algo_config, algo_pcaU_operator, algo_avg_prior, algo_cov_prior, RESHAPE(core_state_slice.measures_(5).measure_data, [num]), global_seed, global_i_real, i_t, next_core_state)
            
            do j = 1, size(input_core_state.measures_)
                forecast_at_t.measures_(j).measure_data(i_idx:i_idx,1:1,:) = next_core_state.measures_(j).measure_data
            end do
            call coef_print(next_core_state, 1, to_print)
            write(10, '(a,i2,a,a)') 'Process-'//'rank:', rank, '  ', trim(to_print)
            write(*, '(a,i2,a,a)') 'Process-'//'rank:', rank, '  ', trim(to_print)
        end do        
    end subroutine parallel_forecast_step_AR1
!==========================================================================================================================

!==========================================================================================================================
    subroutine parallel_forecast_step_AR3(self, algo_config, algo_nb_realisations, algo_attributed_models, algo_pcaU_operator, algo_avg_prior, algo_cov_prior, input_core_state, Z_AR3, global_seed, i_t, forecast_at_t)
    !*****************************************************************************************************************
    !"""
    !Parallelize the AR-3 forecast step over the realisations
    !attributed to the current MPI process.
    !
    !:param input_core_state: input core state at time t
    !:type input_core_state: corestates.CoreState
    !:param Z_AR3: AR-3 forecast history
    !:type Z_AR3: np.array (Nlocal_real x 3 x Ncoef)
    !:param global_seed: global random seed
    !:type global_seed: int
    !:param i_t: time index
    !:type i_t: int
    !:return: CoreState containing the forecast result
    !:rtype: corestates.CoreState
    !"""
    !*****************************************************************************************************************
        class(AugkfForecasterAR3), intent(inout) :: self
        class(ComputationConfig), intent(in) :: algo_config
        integer, intent(in) :: algo_nb_realisations
        integer, intent(in) :: algo_attributed_models(:)
        class(NormedPCAOperator), intent(in) :: algo_pcaU_operator
        class(set_prior_type), intent(in) :: algo_avg_prior
        class(cov_prior_type), intent(in) :: algo_cov_prior
        class(CoreState_type), intent(in) :: input_core_state
        real(kind=8), intent(inout) :: Z_AR3(:,:,:)
        integer, intent(in) :: global_seed, i_t
        class(CoreState_type), allocatable, intent(out) :: forecast_at_t

        real(kind=8) :: t
        type(CoreState_type) :: core_state_slice
        integer :: i_idx, global_i_real, j, i, comm, rank, ierr
        class(CoreState_type), allocatable :: next_core_state
        character(len=1000) :: to_print

        t = algo_config.t_forecasts(i_t+1)

        if (SIZE(Z_AR3, 1) .ne. SIZE(algo_attributed_models)) then
            write (10,'(A)') &
                'Z_AR3 must contain one history for every attributed model.'
            write (*,'(A)') &
                'Z_AR3 must contain one history for every attributed model.'
            stop
        end if

        !# Copy input core state
        allocate(forecast_at_t, source=input_core_state)

        !# Set all measures to zero
        do i = 1, SIZE(forecast_at_t.measures_)
            forecast_at_t.measures_(i).measure_data = 0.0d0
        end do

        !# Allocate a one-realisation core-state container
        allocate(core_state_slice.measures_, source=input_core_state.measures_)
        allocate(core_state_slice.max_degrees_, source=input_core_state.max_degrees_)

        !# MPI rank
        comm = MPI_COMM_WORLD
        call MPI_Comm_rank(comm, rank, ierr)

        !# Each process computes its attributed models
        do i_idx = 1, SIZE(algo_attributed_models)
            global_i_real = algo_attributed_models(i_idx)

            !# Slice the core state for the local model
            do j = 1, SIZE(input_core_state.measures_)
                deallocate(core_state_slice.measures_(j).measure_data)
                allocate( &
                    core_state_slice.measures_(j).measure_data( &
                        1, &
                        1, &
                        SIZE(input_core_state.measures_(j).measure_data, 3) &
                    ) &
                )
                core_state_slice.measures_(j).measure_data = &
                    input_core_state.measures_(j).measure_data(i_idx:i_idx,:,:)
            end do

            call self.forecast_step_AR3( &
                core_state_slice,              &
                algo_nb_realisations,          &
                algo_config,                   &
                algo_pcaU_operator,            &
                algo_avg_prior,                &
                algo_cov_prior,                &
                Z_AR3(i_idx,:,:),              &
                global_seed,                   &
                global_i_real,                 &
                i_t,                           &
                next_core_state                &
            )

            do j = 1, SIZE(input_core_state.measures_)
                forecast_at_t.measures_(j).measure_data(i_idx:i_idx,1:1,:) = &
                    next_core_state.measures_(j).measure_data
            end do

            call coef_print(next_core_state, 1, to_print)
            write (10,'(a,i2,a,a)') &
                'Process-'//'rank:', rank, '  ', TRIM(to_print)
            write (*,'(a,i2,a,a)') &
                'Process-'//'rank:', rank, '  ', TRIM(to_print)
        end do

        !# Update the three-state AR-3 history
        call self.update_Z_AR3( &
            forecast_at_t.measures_(5).measure_data(:,1,:), &
            Z_AR3                                           &
        )

    end subroutine parallel_forecast_step_AR3
!==========================================================================================================================
    
end module
