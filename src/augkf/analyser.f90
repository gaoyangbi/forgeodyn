module analyser
    use mpi
    use common
    use computer
    use observations
    use blas95
    use lapack95
    use f95_precision
    use, intrinsic :: ieee_arithmetic
    use config
    use corestate
    use pca
    implicit none
    
    
    type, extends(GenericComputer),public :: AugkfAnalyserAR1
        character(len=10) :: type_(2)
        class(Observation), allocatable :: measure_observations_SV(:)
        class(Observation), allocatable :: measure_observations_MF(:)
        real(kind=8) :: current_misfits(5)
        character(len=10) :: keys(5) 
        logical :: do_backward_analysis
        logical, allocatable :: ana_sv(:), ana_mf(:)
        class(measure_observations_mat), allocatable :: sv_X(:), sv_H(:), sv_RXX(:)
        class(measure_observations_mat), allocatable :: mf_X(:), mf_H(:), mf_RXX(:)
    contains
        procedure :: init_AugkfAnalyserAR, invalid_misfits
        procedure :: extract_observations, check_if_analysis_data
        procedure :: sv_analysis, mf_analysis, analysis_step, analyse_B, analyse_Z
        procedure :: remove_small_correlations, setup_Hz
    end type AugkfAnalyserAR1
    
    type, extends(AugkfAnalyserAR1),public :: AugkfAnalyserAR3
        
    contains
        procedure :: compute_full_state_AR3
        procedure :: build_H_operator_same_AR3
        procedure :: build_H_operator_different_AR3
        generic :: build_H_operator_AR3 => build_H_operator_same_AR3, build_H_operator_different_AR3
        procedure :: build_obs_operator_AR3
        procedure :: build_err_operator_AR3
        procedure :: analyse_B_AR3
        procedure :: build_full_Hz_AR3
        procedure :: analyse_Z_AR3
        procedure :: dZ_to_dU_dER_AR3
        procedure :: init_forecast_AR3
        procedure :: analysis_step_AR3
    end type AugkfAnalyserAR3
    
contains    
!==========================================================================================================================
    subroutine compute_full_state_AR3(self, X_state, algo_cfg, nb_realisations, attributed_models, X_full)
    !*****************************************************************************************************************
    !"""
    !Compute the least-squares fit of X, dX/dt and d2X/dt2
    !at the analysis time from a forecast time window.
    !
    !:param X_state: forecast states around the analysis time
    !:type X_state: np.array (Nreal x Ntimes x Ncoef)
    !:param nb_realisations: total number of realisations
    !:type nb_realisations: int
    !:param attributed_models: global model indices attributed to this MPI process
    !:type attributed_models: np.array
    !:return: fitted states [X, dX/dt, d2X/dt2]
    !:rtype: np.array (Nreal x 3*Ncoef)
    !"""
    !*****************************************************************************************************************
        class(AugkfAnalyserAR3), intent(in) :: self
        real(kind=8), intent(in) :: X_state(:,:,:)
        class(ComputationConfig), intent(in) :: algo_cfg
        integer, intent(in) :: nb_realisations
        integer, intent(in) :: attributed_models(:)
        real(kind=8), allocatable, intent(out) :: X_full(:,:)

        integer :: ratio, N_times, N_coef
        integer :: i, i_idx, global_idx, ierr
        real(kind=8), allocatable :: dt(:), G(:,:), GTG(:,:)
        real(kind=8), allocatable :: GTG_inv(:,:), fit_operator(:,:), S(:,:)

        ratio = algo_cfg.dt_a_f_ratio
        N_times = SIZE(X_state, 2)
        N_coef = SIZE(X_state, 3)

        if (N_times .ne. 2*ratio + 1) then
            write (10,'(A, i4, A, i4, A)') &
                'AR3 full-state fit requires ', 2*ratio + 1, &
                ' forecast times. Got ', N_times, ' instead.'
            write (*,'(A, i4, A, i4, A)') &
                'AR3 full-state fit requires ', 2*ratio + 1, &
                ' forecast times. Got ', N_times, ' instead.'
            stop
        end if

        allocate(dt(N_times))

        if (self.do_backward_analysis) then
            do i = 1, N_times
                dt(i) = REAL(i - N_times, kind=8) * algo_cfg.dt_f
            end do
        else
            do i = 1, N_times
                dt(i) = REAL(i - ratio - 1, kind=8) * algo_cfg.dt_f
            end do
        end if

        !# Build the second-order Taylor operator
        allocate(G(N_times, 3), source=0.0d0)
        G(:,1) = 1.0d0
        G(:,2) = dt
        G(:,3) = dt**2 / 2.0d0

        !# Least-squares operator: inv(G.T G) G.T
        allocate(GTG, source=MATMUL(TRANSPOSE(G), G))
        call max_inv(GTG, GTG_inv)
        allocate(fit_operator, source=MATMUL(GTG_inv, TRANSPOSE(G)))

        allocate(X_full(nb_realisations, 3*N_coef), source=0.0d0)
        allocate(S(3, N_coef), source=0.0d0)

        do i_idx = 1, SIZE(attributed_models)
            global_idx = attributed_models(i_idx) + 1
            S = MATMUL(fit_operator, X_state(global_idx,:,:))

            !# Preserve NumPy S[:3].flatten() ordering:
            !# [X coefficients, dX/dt coefficients, d2X/dt2 coefficients]
            X_full(global_idx,1:N_coef) = S(1,:)
            X_full(global_idx,N_coef+1:2*N_coef) = S(2,:)
            X_full(global_idx,2*N_coef+1:3*N_coef) = S(3,:)
        end do

        call MPI_ALLREDUCE( &
            MPI_IN_PLACE,         &
            X_full,               &
            SIZE(X_full),         &
            MPI_DOUBLE_PRECISION, &
            MPI_SUM,              &
            MPI_COMM_WORLD,       &
            ierr                  &
        )

    end subroutine compute_full_state_AR3
!==========================================================================================================================

!==========================================================================================================================
    subroutine build_H_operator_same_AR3(self, X, H_operator)
    !*****************************************************************************************************************
    !"""
    !Build the AR3 observation operator when the same matrix X is used
    !at the three analysis times (Python row_diff=False).
    !"""
    !*****************************************************************************************************************
        class(AugkfAnalyserAR3), intent(in) :: self
        real(kind=8), intent(in) :: X(:,:)
        real(kind=8), allocatable, intent(out) :: H_operator(:,:)

        integer :: N_rows, N_cols
        real(kind=8) :: dt_a

        N_rows = SIZE(X, 1)
        N_cols = SIZE(X, 2)
        dt_a = self.cfg.dt_a

        allocate(H_operator(3*N_rows, 3*N_cols), source=0.0d0)

        if (.not. self.do_backward_analysis) then
            H_operator(1:N_rows,1:N_cols) = X
            H_operator(1:N_rows,N_cols+1:2*N_cols) = -dt_a * X
            H_operator(1:N_rows,2*N_cols+1:3*N_cols) = dt_a**2 / 2.0d0 * X

            H_operator(N_rows+1:2*N_rows,1:N_cols) = X

            H_operator(2*N_rows+1:3*N_rows,1:N_cols) = X
            H_operator(2*N_rows+1:3*N_rows,N_cols+1:2*N_cols) = dt_a * X
            H_operator(2*N_rows+1:3*N_rows,2*N_cols+1:3*N_cols) = dt_a**2 / 2.0d0 * X
        else
            H_operator(1:N_rows,1:N_cols) = X
            H_operator(1:N_rows,N_cols+1:2*N_cols) = -2.0d0 * dt_a * X
            H_operator(1:N_rows,2*N_cols+1:3*N_cols) = 2.0d0 * dt_a**2 * X

            H_operator(N_rows+1:2*N_rows,1:N_cols) = X
            H_operator(N_rows+1:2*N_rows,N_cols+1:2*N_cols) = -dt_a * X
            H_operator(N_rows+1:2*N_rows,2*N_cols+1:3*N_cols) = dt_a**2 / 2.0d0 * X

            H_operator(2*N_rows+1:3*N_rows,1:N_cols) = X
        end if

    end subroutine build_H_operator_same_AR3
!==========================================================================================================================

!==========================================================================================================================
    subroutine build_H_operator_different_AR3(self, X, H_operator)
    !*****************************************************************************************************************
    !"""
    !Build the AR3 observation operator when a different matrix X is used
    !at each of the three analysis times (Python row_diff=True).
    !"""
    !*****************************************************************************************************************
        class(AugkfAnalyserAR3), intent(in) :: self
        real(kind=8), intent(in) :: X(:,:,:)
        real(kind=8), allocatable, intent(out) :: H_operator(:,:)

        integer :: N_rows, N_cols
        real(kind=8) :: dt_a

        if (SIZE(X, 1) .ne. 3) then
            write (10,'(A, i4, A)') &
                'AR3 H operator requires 3 input matrices. Got ', SIZE(X, 1), ' instead.'
            write (*,'(A, i4, A)') &
                'AR3 H operator requires 3 input matrices. Got ', SIZE(X, 1), ' instead.'
            stop
        end if

        N_rows = SIZE(X, 2)
        N_cols = SIZE(X, 3)
        dt_a = self.cfg.dt_a

        allocate(H_operator(3*N_rows, 3*N_cols), source=0.0d0)

        if (.not. self.do_backward_analysis) then
            H_operator(1:N_rows,1:N_cols) = X(1,:,:)
            H_operator(1:N_rows,N_cols+1:2*N_cols) = -dt_a * X(1,:,:)
            H_operator(1:N_rows,2*N_cols+1:3*N_cols) = dt_a**2 / 2.0d0 * X(1,:,:)

            H_operator(N_rows+1:2*N_rows,1:N_cols) = X(2,:,:)

            H_operator(2*N_rows+1:3*N_rows,1:N_cols) = X(3,:,:)
            H_operator(2*N_rows+1:3*N_rows,N_cols+1:2*N_cols) = dt_a * X(3,:,:)
            H_operator(2*N_rows+1:3*N_rows,2*N_cols+1:3*N_cols) = dt_a**2 / 2.0d0 * X(3,:,:)
        else
            H_operator(1:N_rows,1:N_cols) = X(1,:,:)
            H_operator(1:N_rows,N_cols+1:2*N_cols) = -2.0d0 * dt_a * X(1,:,:)
            H_operator(1:N_rows,2*N_cols+1:3*N_cols) = 2.0d0 * dt_a**2 * X(1,:,:)

            H_operator(N_rows+1:2*N_rows,1:N_cols) = X(2,:,:)
            H_operator(N_rows+1:2*N_rows,N_cols+1:2*N_cols) = -dt_a * X(2,:,:)
            H_operator(N_rows+1:2*N_rows,2*N_cols+1:3*N_cols) = dt_a**2 / 2.0d0 * X(2,:,:)

            H_operator(2*N_rows+1:3*N_rows,1:N_cols) = X(3,:,:)
        end if

    end subroutine build_H_operator_different_AR3
!==========================================================================================================================

!==========================================================================================================================
    subroutine build_obs_operator_AR3(self, obs, obs_operator)
    !*****************************************************************************************************************
    !"""
    !Concatenate the observation data at the three AR3 analysis times.
    !Equivalent to NumPy concatenate((obs[0], obs[1], obs[2]), axis=1).
    !"""
    !*****************************************************************************************************************
        class(AugkfAnalyserAR3), intent(in) :: self
        class(measure_observations_mat), intent(in) :: obs(:)
        real(kind=8), allocatable, intent(out) :: obs_operator(:,:)

        integer :: N_rows, N_cols_1, N_cols_2, N_cols_3

        if (SIZE(obs) .ne. 3) then
            write (10,'(A, i4, A)') &
                'AR3 observation operator requires 3 input matrices. Got ', SIZE(obs), ' instead.'
            write (*,'(A, i4, A)') &
                'AR3 observation operator requires 3 input matrices. Got ', SIZE(obs), ' instead.'
            stop
        end if

        N_rows = SIZE(obs(1).mat, 1)
        if (SIZE(obs(2).mat, 1) .ne. N_rows .or. SIZE(obs(3).mat, 1) .ne. N_rows) then
            write (10,'(A)') 'AR3 observation matrices must have the same number of rows.'
            write (*,'(A)') 'AR3 observation matrices must have the same number of rows.'
            stop
        end if

        N_cols_1 = SIZE(obs(1).mat, 2)
        N_cols_2 = SIZE(obs(2).mat, 2)
        N_cols_3 = SIZE(obs(3).mat, 2)

        allocate(obs_operator(N_rows, N_cols_1+N_cols_2+N_cols_3))
        obs_operator(:,1:N_cols_1) = obs(1).mat
        obs_operator(:,N_cols_1+1:N_cols_1+N_cols_2) = obs(2).mat
        obs_operator(:,N_cols_1+N_cols_2+1:N_cols_1+N_cols_2+N_cols_3) = obs(3).mat

    end subroutine build_obs_operator_AR3
!==========================================================================================================================

!==========================================================================================================================
    subroutine build_err_operator_AR3(self, error, err_operator)
    !*****************************************************************************************************************
    !"""
    !Build the block-diagonal observation-error matrix for the three
    !AR3 analysis times. Equivalent to scipy.linalg.block_diag.
    !"""
    !*****************************************************************************************************************
        class(AugkfAnalyserAR3), intent(in) :: self
        class(measure_observations_mat), intent(in) :: error(:)
        real(kind=8), allocatable, intent(out) :: err_operator(:,:)

        integer :: N_1, N_2, N_3

        if (SIZE(error) .ne. 3) then
            write (10,'(A, i4, A)') &
                'AR3 error operator requires 3 input matrices. Got ', SIZE(error), ' instead.'
            write (*,'(A, i4, A)') &
                'AR3 error operator requires 3 input matrices. Got ', SIZE(error), ' instead.'
            stop
        end if

        if (SIZE(error(1).mat, 1) .ne. SIZE(error(1).mat, 2) .or. &
            SIZE(error(2).mat, 1) .ne. SIZE(error(2).mat, 2) .or. &
            SIZE(error(3).mat, 1) .ne. SIZE(error(3).mat, 2)) then
            write (10,'(A)') 'AR3 observation-error matrices must be square.'
            write (*,'(A)') 'AR3 observation-error matrices must be square.'
            stop
        end if

        N_1 = SIZE(error(1).mat, 1)
        N_2 = SIZE(error(2).mat, 1)
        N_3 = SIZE(error(3).mat, 1)

        allocate(err_operator(N_1+N_2+N_3, N_1+N_2+N_3), source=0.0d0)
        err_operator(1:N_1,1:N_1) = error(1).mat
        err_operator(N_1+1:N_1+N_2,N_1+1:N_1+N_2) = error(2).mat
        err_operator(N_1+N_2+1:N_1+N_2+N_3,N_1+N_2+1:N_1+N_2+N_3) = error(3).mat

    end subroutine build_err_operator_AR3
!==========================================================================================================================

!==========================================================================================================================
    subroutine analyse_B_AR3(self, B_full, algo_cfg, nb_realisations, attributed_models, &
                             b_minus, b, b_plus)
    !*****************************************************************************************************************
    !"""
    !Analyse the complete AR3 magnetic-field state [B, dB/dt, d2B/dt2]
    !using observations at three analysis times, then reconstruct
    !B at the previous, current and following analysis times.
    !"""
    !*****************************************************************************************************************
        class(AugkfAnalyserAR3), intent(inout) :: self
        real(kind=8), intent(in) :: B_full(:,:)
        class(ComputationConfig), intent(in) :: algo_cfg
        integer, intent(in) :: nb_realisations
        integer, intent(in) :: attributed_models(:)
        real(kind=8), allocatable, intent(out) :: b_minus(:,:), b(:,:), b_plus(:,:)

        real(kind=8), allocatable :: mf_H_3D(:,:,:)
        real(kind=8), allocatable :: Hb(:,:), Yb(:,:), Rbb(:,:)
        real(kind=8), allocatable :: Pbb(:,:), Kbb(:,:), analysed_B(:,:)
        real(kind=8), allocatable :: P_eig_val(:), P_eig_vec(:,:)
        real(kind=8), allocatable :: HX_b(:,:), inv_Rbb(:,:)
        real(kind=8), allocatable :: identity_B(:,:), time_operator(:,:), B_a(:,:)
        integer :: Nb, N_H_rows, N_H_cols
        integer :: i, i_real, ierr, info

        Nb = algo_cfg.Nb()

        if (SIZE(B_full, 1) .ne. nb_realisations .or. SIZE(B_full, 2) .ne. 3*Nb) then
            write (10,'(A)') 'AR3 B_full must have dimensions (nb_realisations, 3*Nb).'
            write (*,'(A)') 'AR3 B_full must have dimensions (nb_realisations, 3*Nb).'
            stop
        end if

        if (SIZE(self.mf_H) .ne. 3) then
            write (10,'(A)') 'AR3 B analysis requires MF observation operators at 3 times.'
            write (*,'(A)') 'AR3 B analysis requires MF observation operators at 3 times.'
            stop
        end if

        N_H_rows = SIZE(self.mf_H(1).mat, 1)
        N_H_cols = SIZE(self.mf_H(1).mat, 2)
        if (SIZE(self.mf_H(2).mat, 1) .ne. N_H_rows .or. &
            SIZE(self.mf_H(3).mat, 1) .ne. N_H_rows .or. &
            SIZE(self.mf_H(2).mat, 2) .ne. N_H_cols .or. &
            SIZE(self.mf_H(3).mat, 2) .ne. N_H_cols) then
            write (10,'(A)') 'The 3 AR3 MF observation operators must have matching dimensions.'
            write (*,'(A)') 'The 3 AR3 MF observation operators must have matching dimensions.'
            stop
        end if

        allocate(mf_H_3D(3, N_H_rows, N_H_cols))
        mf_H_3D(1,:,:) = self.mf_H(1).mat
        mf_H_3D(2,:,:) = self.mf_H(2).mat
        mf_H_3D(3,:,:) = self.mf_H(3).mat

        !# Observation operator, observation data and block-diagonal error covariance
        call self.build_H_operator_AR3(mf_H_3D, Hb)
        call self.build_obs_operator_AR3(self.mf_X, Yb)
        call self.build_err_operator_AR3(self.mf_Rxx, Rbb)

        if (SIZE(Hb, 2) .ne. SIZE(B_full, 2) .or. &
            SIZE(Yb, 1) .ne. nb_realisations .or. &
            SIZE(Yb, 2) .ne. SIZE(Hb, 1) .or. &
            SIZE(Rbb, 1) .ne. SIZE(Hb, 1)) then
            write (10,'(A)') 'Inconsistent matrix dimensions in AR3 B analysis.'
            write (*,'(A)') 'Inconsistent matrix dimensions in AR3 B analysis.'
            stop
        end if

        !# Forecast covariance of the complete state [B, dB/dt, d2B/dt2]
        call self.remove_small_correlations(B_full, 1.0d-10, algo_cfg, Pbb)

        write (*,*) 'Getting best linear unbiased estimate of AR3 B...'
        write (10,*) 'Getting best linear unbiased estimate of AR3 B...'
        allocate(analysed_B(nb_realisations, 3*Nb), source=0.0d0)

        if (TRIM(algo_cfg.kalman_norm) == 'l2') then
            call compute_Kalman_gain_matrix(Pbb, Hb, Rbb, .True., Kbb)

            do i = 1, SIZE(attributed_models)
                i_real = attributed_models(i) + 1
                call get_BLUE(B_full(i_real,:), Yb(i_real,:), Pbb, Hb, Rbb, &
                              Kbb, .True., analysed_B(i_real,:))
            end do

        else if (TRIM(algo_cfg.kalman_norm) == 'huber') then
            allocate(P_eig_val(SIZE(Pbb, 1)))
            allocate(P_eig_vec, source=Pbb)
            call syevd(P_eig_vec, P_eig_val, 'V', 'U', info)
            if (info .ne. 0) then
                write (10,'(A, i6)') 'Eigenvalue decomposition failed in AR3 B analysis. INFO = ', info
                write (*,'(A, i6)') 'Eigenvalue decomposition failed in AR3 B analysis. INFO = ', info
                stop
            end if
            P_eig_val = MAX(P_eig_val, 1.0d-10)

            do i = 1, SIZE(attributed_models)
                i_real = attributed_models(i) + 1
                call compute_Kalman_huber(B_full(i_real,:), Yb(i_real,:), &
                                           P_eig_val, P_eig_vec, Hb, Rbb, 50, &
                                           analysed_B(i_real,:))
            end do

        else
            write (10,'(A,A)') 'Invalid kalman_norm in AR3 B analysis: ', TRIM(algo_cfg.kalman_norm)
            write (*,'(A,A)') 'Invalid kalman_norm in AR3 B analysis: ', TRIM(algo_cfg.kalman_norm)
            stop
        end if

        !# Each MPI process filled only its attributed global-realisation rows.
        call MPI_ALLREDUCE( &
            MPI_IN_PLACE,         &
            analysed_B,           &
            SIZE(analysed_B),     &
            MPI_DOUBLE_PRECISION, &
            MPI_SUM,              &
            MPI_COMM_WORLD,       &
            ierr                  &
        )

        !# MF misfit in the combined three-time observation space
        allocate(HX_b, source=TRANSPOSE(MATMUL(Hb, TRANSPOSE(analysed_B))))
        call max_inv(Rbb, inv_Rbb)
        call compute_misfit(Yb, HX_b, inv_Rbb, self.current_misfits(1))

        !# Reconstruct [B(t-), B(t), B(t+)] from [B, dB/dt, d2B/dt2].
        allocate(identity_B(Nb, Nb), source=0.0d0)
        do i = 1, Nb
            identity_B(i,i) = 1.0d0
        end do
        call self.build_H_operator_AR3(identity_B, time_operator)
        allocate(B_a, source=MATMUL(analysed_B, TRANSPOSE(time_operator)))

        allocate(b_minus(nb_realisations, Nb), source=B_a(:,1:Nb))
        allocate(b(nb_realisations, Nb), source=B_a(:,Nb+1:2*Nb))
        allocate(b_plus(nb_realisations, Nb), source=B_a(:,2*Nb+1:3*Nb))

    end subroutine analyse_B_AR3
!==========================================================================================================================

!==========================================================================================================================
    subroutine build_full_Hz_AR3(self, input_core_state, b_minus, b, b_plus, sv_H, &
                                 algo_cfg, pcaU_operator, H, Ab_out)
    !*****************************************************************************************************************
    !"""
    !Build the complete Z observation matrices at the three AR3 analysis
    !times using B(t-), B(t) and B(t+) respectively.
    !"""
    !*****************************************************************************************************************
        class(AugkfAnalyserAR3), intent(in) :: self
        class(CoreState_type), intent(inout) :: input_core_state
        real(kind=8), intent(in) :: b_minus(:), b(:), b_plus(:)
        class(measure_observations_mat), intent(in) :: sv_H(:)
        class(ComputationConfig), intent(in) :: algo_cfg
        class(NormedPCAOperator), intent(in) :: pcaU_operator
        real(kind=8), allocatable, intent(out) :: H(:,:,:)
        real(kind=8), allocatable, intent(out) :: Ab_out(:,:,:)

        type(input_core_state_type) :: CoreState_temp
        real(kind=8), allocatable :: B_at_time(:)
        real(kind=8), allocatable :: Ab(:,:), S_u(:,:), Hz(:,:), complete_H(:,:)
        integer :: i, N_obs, N_state

        if (SIZE(sv_H) .ne. 3) then
            write (10,'(A)') 'AR3 full Hz construction requires SV operators at 3 times.'
            write (*,'(A)') 'AR3 full Hz construction requires SV operators at 3 times.'
            stop
        end if

        if (SIZE(b_minus) .ne. algo_cfg.Nb() .or. &
            SIZE(b) .ne. algo_cfg.Nb() .or. &
            SIZE(b_plus) .ne. algo_cfg.Nb()) then
            write (10,'(A)') 'Invalid magnetic-field dimensions in AR3 full Hz construction.'
            write (*,'(A)') 'Invalid magnetic-field dimensions in AR3 full Hz construction.'
            stop
        end if

        N_obs = SIZE(sv_H(1).mat, 1)
        if (SIZE(sv_H(2).mat, 1) .ne. N_obs .or. &
            SIZE(sv_H(3).mat, 1) .ne. N_obs .or. &
            SIZE(sv_H(1).mat, 2) .ne. algo_cfg.Nsv() .or. &
            SIZE(sv_H(2).mat, 2) .ne. algo_cfg.Nsv() .or. &
            SIZE(sv_H(3).mat, 2) .ne. algo_cfg.Nsv()) then
            write (10,'(A)') 'Inconsistent SV operator dimensions in AR3 full Hz construction.'
            write (*,'(A)') 'Inconsistent SV operator dimensions in AR3 full Hz construction.'
            stop
        end if

        if (algo_cfg.pca == 1) then
            N_state = algo_cfg.N_pca_u + algo_cfg.Nsv()
            call pcaU_operator.S_u(S_u)
            if (SIZE(S_u, 1) .ne. algo_cfg.Nu2() .or. &
                SIZE(S_u, 2) .ne. algo_cfg.N_pca_u) then
                write (10,'(A)') 'Inconsistent PCA operator dimensions in AR3 full Hz construction.'
                write (*,'(A)') 'Inconsistent PCA operator dimensions in AR3 full Hz construction.'
                stop
            end if
        else
            N_state = algo_cfg.Nu2() + algo_cfg.Nsv()
        end if

        allocate(H(3, N_obs, N_state), source=0.0d0)
        allocate(Ab_out(3, algo_cfg.Nsv(), algo_cfg.Nu2()), source=0.0d0)
        allocate(B_at_time(algo_cfg.Nb()))

        CoreState_temp.Lsv = input_core_state.cs_Lsv()
        CoreState_temp.Lu = input_core_state.cs_Lu()
        CoreState_temp.Lb = input_core_state.cs_Lb()
        CoreState_temp.Nsv = input_core_state.cs_Nsv()
        CoreState_temp.Nu2 = input_core_state.cs_Nu2()
        CoreState_temp.Nb = input_core_state.cs_Nb()
        allocate(CoreState_temp.B(algo_cfg.Nb()))

        do i = 1, 3
            select case (i)
            case (1)
                B_at_time = b_minus
            case (2)
                B_at_time = b
            case (3)
                B_at_time = b_plus
            end select

            CoreState_temp.B = B_at_time
            call self.compute_Ab(CoreState_temp, Ab)
            Ab_out(i,:,:) = Ab

            if (algo_cfg.pca == 1) then
                call self.setup_Hz(MATMUL(Ab, S_u), algo_cfg.N_pca_u, algo_cfg, Hz)
            else
                call self.setup_Hz(Ab, algo_cfg.Nu2(), algo_cfg, Hz)
            end if

            allocate(complete_H, source=MATMUL(sv_H(i).mat, Hz))
            H(i,:,:) = complete_H
            deallocate(complete_H)
        end do

    end subroutine build_full_Hz_AR3
!==========================================================================================================================

!==========================================================================================================================
    subroutine analyse_Z_AR3(self, Z_full, input_core_state, algo_cfg, nb_realisations, &
                             attributed_models, pcaU_operator, algo_avg_prior,          &
                             b_minus, b, b_plus, ana_Z, ana_dZ, ana_d2Z,                &
                             analysed_U, analysed_ER, analysed_SV)
    !*****************************************************************************************************************
    !"""
    !Analyse the complete AR3 state [Z, dZ/dt, d2Z/dt2] using SV
    !observations at three analysis times.
    !"""
    !*****************************************************************************************************************
        class(AugkfAnalyserAR3), intent(inout) :: self
        real(kind=8), intent(in) :: Z_full(:,:)
        class(CoreState_type), intent(inout) :: input_core_state
        class(ComputationConfig), intent(in) :: algo_cfg
        integer, intent(in) :: nb_realisations
        integer, intent(in) :: attributed_models(:)
        class(NormedPCAOperator), intent(in) :: pcaU_operator
        class(set_prior_type), intent(in) :: algo_avg_prior
        real(kind=8), intent(in) :: b_minus(:,:), b(:,:), b_plus(:,:)
        real(kind=8), allocatable, intent(out) :: ana_Z(:,:), ana_dZ(:,:), ana_d2Z(:,:)
        real(kind=8), allocatable, intent(out) :: analysed_U(:,:), analysed_ER(:,:), analysed_SV(:,:)

        real(kind=8), allocatable :: Rzz(:,:), PZZ(:,:)
        real(kind=8), allocatable :: H_3D(:,:,:), Ab_3D(:,:,:), Hz(:,:)
        real(kind=8), allocatable :: PzzHT(:,:), HPzzHT(:,:)
        real(kind=8), allocatable :: analysed_Z_full(:,:), Yz(:,:), HX_z(:,:), inv_Rzz(:,:)
        real(kind=8), allocatable :: avg_U(:), avg_ER(:), mean_SV(:)
        integer :: Nz, Nsv, Nu2, N_obs, N_total_obs
        integer :: i, i_idx, i_real, col_start, col_end, ierr

        Nz = algo_cfg.Nz()
        Nsv = algo_cfg.Nsv()
        Nu2 = algo_cfg.Nu2()

        if (SIZE(Z_full, 1) .ne. nb_realisations .or. SIZE(Z_full, 2) .ne. 3*Nz) then
            write (10,'(A)') 'AR3 Z_full must have dimensions (nb_realisations, 3*Nz).'
            write (*,'(A)') 'AR3 Z_full must have dimensions (nb_realisations, 3*Nz).'
            stop
        end if

        if (SIZE(b_minus, 1) .ne. nb_realisations .or. &
            SIZE(b, 1) .ne. nb_realisations .or. &
            SIZE(b_plus, 1) .ne. nb_realisations .or. &
            SIZE(b_minus, 2) .ne. algo_cfg.Nb() .or. &
            SIZE(b, 2) .ne. algo_cfg.Nb() .or. &
            SIZE(b_plus, 2) .ne. algo_cfg.Nb()) then
            write (10,'(A)') 'Invalid magnetic-field dimensions in AR3 Z analysis.'
            write (*,'(A)') 'Invalid magnetic-field dimensions in AR3 Z analysis.'
            stop
        end if

        if (SIZE(self.sv_X) .ne. 3 .or. SIZE(self.sv_H) .ne. 3 .or. &
            SIZE(self.sv_Rxx) .ne. 3) then
            write (10,'(A)') 'AR3 Z analysis requires SV observations at 3 times.'
            write (*,'(A)') 'AR3 Z analysis requires SV observations at 3 times.'
            stop
        end if

        N_obs = SIZE(self.sv_X(1).mat, 2)
        if (SIZE(self.sv_X(2).mat, 2) .ne. N_obs .or. &
            SIZE(self.sv_X(3).mat, 2) .ne. N_obs) then
            write (10,'(A)') 'The 3 AR3 SV observation vectors must have matching dimensions.'
            write (*,'(A)') 'The 3 AR3 SV observation vectors must have matching dimensions.'
            stop
        end if
        N_total_obs = 3*N_obs

        call self.build_err_operator_AR3(self.sv_Rxx, Rzz)
        if (SIZE(Rzz, 1) .ne. N_total_obs) then
            write (10,'(A)') 'SV observation and error dimensions do not match in AR3 Z analysis.'
            write (*,'(A)') 'SV observation and error dimensions do not match in AR3 Z analysis.'
            stop
        end if

        !# Forecast covariance of the complete state [Z, dZ/dt, d2Z/dt2]
        call self.remove_small_correlations(Z_full, 1.0d-10, algo_cfg, PZZ)

        allocate(analysed_Z_full(nb_realisations, 3*Nz), source=0.0d0)
        allocate(analysed_U(nb_realisations, Nu2), source=0.0d0)
        allocate(analysed_ER(nb_realisations, Nsv), source=0.0d0)
        allocate(analysed_SV(nb_realisations, Nsv), source=0.0d0)
        allocate(Yz(nb_realisations, N_total_obs), source=0.0d0)
        allocate(HX_z(nb_realisations, N_total_obs), source=0.0d0)

        allocate(avg_U, source=RESHAPE(algo_avg_prior.U, [SIZE(algo_avg_prior.U)]))
        allocate(avg_ER, source=RESHAPE(algo_avg_prior.ER, [SIZE(algo_avg_prior.ER)]))
        if (SIZE(avg_U) .ne. Nu2 .or. SIZE(avg_ER) .ne. Nsv) then
            write (10,'(A)') 'Invalid prior-mean dimensions in AR3 Z analysis.'
            write (*,'(A)') 'Invalid prior-mean dimensions in AR3 Z analysis.'
            stop
        end if
        allocate(mean_SV(N_obs))

        do i_idx = 1, SIZE(attributed_models)
            i_real = attributed_models(i_idx) + 1

            call self.build_full_Hz_AR3( &
                input_core_state,         &
                b_minus(i_real,:),        &
                b(i_real,:),              &
                b_plus(i_real,:),         &
                self.sv_H,                &
                algo_cfg,                 &
                pcaU_operator,            &
                H_3D,                     &
                Ab_3D                     &
            )
            call self.build_H_operator_AR3(H_3D, Hz)

            if (SIZE(Hz, 1) .ne. N_total_obs .or. SIZE(Hz, 2) .ne. 3*Nz) then
                write (10,'(A)') 'Invalid complete observation-operator dimensions in AR3 Z analysis.'
                write (*,'(A)') 'Invalid complete observation-operator dimensions in AR3 Z analysis.'
                stop
            end if

            !# Z is centred on zero. Remove the physical prior means U0 and ER0
            !# independently at each of the three observation times.
            do i = 1, 3
                col_start = (i-1)*N_obs + 1
                col_end = i*N_obs
                mean_SV = MATMUL(self.sv_H(i).mat, MATMUL(Ab_3D(i,:,:), avg_U) + avg_ER)
                Yz(i_real,col_start:col_end) = self.sv_X(i).mat(i_real,:) - mean_SV
            end do

            allocate(PzzHT, source=MATMUL(PZZ, TRANSPOSE(Hz)))
            allocate(HPzzHT, source=MATMUL(Hz, PzzHT))

            !# As in Python AR1/AR3, Z analysis is intentionally fixed to Huber.
            call compute_Kalman_huber_parameter_basis( &
                Z_full(i_real,:),                       &
                Yz(i_real,:),                           &
                HPzzHT,                                 &
                PzzHT,                                  &
                Hz,                                     &
                Rzz,                                    &
                'huber',                                &
                50,                                     &
                analysed_Z_full(i_real,:)                &
            )

            call Z_to_U_ER1( &
                algo_cfg,                               &
                algo_avg_prior,                         &
                pcaU_operator,                          &
                analysed_Z_full(i_real,1:Nz),           &
                analysed_U(i_real,:),                   &
                analysed_ER(i_real,:)                   &
            )

            if (self.do_backward_analysis) then
                analysed_SV(i_real,:) = MATMUL(Ab_3D(3,:,:), analysed_U(i_real,:)) &
                                      + analysed_ER(i_real,:)
            else
                analysed_SV(i_real,:) = MATMUL(Ab_3D(2,:,:), analysed_U(i_real,:)) &
                                      + analysed_ER(i_real,:)
            end if

            !# Hz depends on B and therefore on the realization. Keep the
            !# corresponding predicted observations for a correct global misfit.
            HX_z(i_real,:) = MATMUL(Hz, analysed_Z_full(i_real,:))

            deallocate(PzzHT, HPzzHT)
        end do

        call MPI_ALLREDUCE(MPI_IN_PLACE, analysed_Z_full, SIZE(analysed_Z_full), &
                           MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, analysed_U, SIZE(analysed_U), &
                           MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, analysed_ER, SIZE(analysed_ER), &
                           MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, analysed_SV, SIZE(analysed_SV), &
                           MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, Yz, SIZE(Yz), &
                           MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, HX_z, SIZE(HX_z), &
                           MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

        call max_inv(Rzz, inv_Rzz)
        call compute_misfit(Yz, HX_z, inv_Rzz, self.current_misfits(2))

        allocate(ana_Z(nb_realisations, Nz), source=analysed_Z_full(:,1:Nz))
        allocate(ana_dZ(nb_realisations, Nz), source=analysed_Z_full(:,Nz+1:2*Nz))
        allocate(ana_d2Z(nb_realisations, Nz), source=analysed_Z_full(:,2*Nz+1:3*Nz))

    end subroutine analyse_Z_AR3
!==========================================================================================================================

!==========================================================================================================================
    subroutine dZ_to_dU_dER_AR3(self, dZ, algo_cfg, pcaU_operator, dU, dER)
    !*****************************************************************************************************************
    !"""
    !Convert a first or second derivative of Z into the corresponding
    !derivatives of U and ER. No prior mean is added to derivatives.
    !"""
    !*****************************************************************************************************************
        class(AugkfAnalyserAR3), intent(in) :: self
        real(kind=8), intent(in) :: dZ(:,:)
        class(ComputationConfig), intent(in) :: algo_cfg
        class(NormedPCAOperator), intent(in) :: pcaU_operator
        real(kind=8), allocatable, intent(out) :: dU(:,:), dER(:,:)

        integer :: N_reals, Nz

        N_reals = SIZE(dZ, 1)
        Nz = algo_cfg.Nz()

        if (SIZE(dZ, 2) .ne. Nz) then
            write (10,'(A)') 'Invalid dZ dimensions in AR3 derivative conversion.'
            write (*,'(A)') 'Invalid dZ dimensions in AR3 derivative conversion.'
            stop
        end if

        if (algo_cfg.pca == 1) then
            call pcaU_operator.inverse_transform_deriv(dU, dZ(:,1:algo_cfg.N_pca_u))
            allocate(dER(N_reals, algo_cfg.Nsv()), &
                     source=dZ(:,algo_cfg.N_pca_u+1:Nz))
        else
            allocate(dU(N_reals, algo_cfg.Nu2()), &
                     source=dZ(:,1:algo_cfg.Nu2()))
            allocate(dER(N_reals, algo_cfg.Nsv()), &
                     source=dZ(:,algo_cfg.Nu2()+1:Nz))
        end if

    end subroutine dZ_to_dU_dER_AR3
!==========================================================================================================================

!==========================================================================================================================
    subroutine init_forecast_AR3(self, ana_Z, ana_dZ, ana_d2Z, algo_cfg, &
                                 attributed_models, Z_AR3)
    !*****************************************************************************************************************
    !"""
    !Initialise the three local AR3 forecast-history states from analysed
    !Z, dZ/dt and d2Z/dt2 at the analysis time.
    !"""
    !*****************************************************************************************************************
        class(AugkfAnalyserAR3), intent(in) :: self
        real(kind=8), intent(in) :: ana_Z(:,:), ana_dZ(:,:), ana_d2Z(:,:)
        class(ComputationConfig), intent(in) :: algo_cfg
        integer, intent(in) :: attributed_models(:)
        real(kind=8), allocatable, intent(out) :: Z_AR3(:,:,:)

        integer :: i_idx, i_real, N_reals, Nz
        real(kind=8) :: dt_f

        N_reals = SIZE(ana_Z, 1)
        Nz = algo_cfg.Nz()
        dt_f = algo_cfg.dt_f

        if (SIZE(ana_Z, 2) .ne. Nz .or. &
            SIZE(ana_dZ, 1) .ne. N_reals .or. SIZE(ana_dZ, 2) .ne. Nz .or. &
            SIZE(ana_d2Z, 1) .ne. N_reals .or. SIZE(ana_d2Z, 2) .ne. Nz) then
            write (10,'(A)') 'Invalid analysed Z derivative dimensions in AR3 forecast initialisation.'
            write (*,'(A)') 'Invalid analysed Z derivative dimensions in AR3 forecast initialisation.'
            stop
        end if

        allocate(Z_AR3(SIZE(attributed_models), 3, Nz), source=0.0d0)

        do i_idx = 1, SIZE(attributed_models)
            i_real = attributed_models(i_idx) + 1
            if (i_real < 1 .or. i_real > N_reals) then
                write (10,'(A)') 'Invalid global realization index in AR3 forecast initialisation.'
                write (*,'(A)') 'Invalid global realization index in AR3 forecast initialisation.'
                stop
            end if

            Z_AR3(i_idx,3,:) = ana_Z(i_real,:)
            Z_AR3(i_idx,2,:) = ana_Z(i_real,:) - dt_f*ana_dZ(i_real,:) &
                             + dt_f**2/2.0d0*ana_d2Z(i_real,:)
            Z_AR3(i_idx,1,:) = ana_Z(i_real,:) - 2.0d0*dt_f*ana_dZ(i_real,:) &
                             + (2.0d0*dt_f)**2/2.0d0*ana_d2Z(i_real,:)
        end do

    end subroutine init_forecast_AR3
!==========================================================================================================================

!==========================================================================================================================
    subroutine analysis_step_AR3(self, input_core_state, Z_AR3, do_backward, algo_cfg, &
                                 nb_realisations, attributed_models, pcaU_operator,     &
                                 algo_avg_prior, ana_core_state_slice)
    !*****************************************************************************************************************
    !"""
    !Perform one AR3 analysis using a complete forecast window and return
    !the analysed state rows attributed to this MPI process.
    !"""
    !*****************************************************************************************************************
        class(AugkfAnalyserAR3), intent(inout) :: self
        class(CoreState_type), intent(inout) :: input_core_state
        real(kind=8), allocatable, intent(inout) :: Z_AR3(:,:,:)
        logical, intent(in) :: do_backward
        class(ComputationConfig), intent(in) :: algo_cfg
        integer, intent(in) :: nb_realisations
        integer, intent(in) :: attributed_models(:)
        class(NormedPCAOperator), intent(in) :: pcaU_operator
        class(set_prior_type), intent(in) :: algo_avg_prior
        class(CoreState_type), intent(inout) :: ana_core_state_slice

        real(kind=8), allocatable :: B_full(:,:), Z_full(:,:)
        real(kind=8), allocatable :: b_minus(:,:), b(:,:), b_plus(:,:)
        real(kind=8), allocatable :: ana_Z(:,:), ana_dZ(:,:), ana_d2Z(:,:)
        real(kind=8), allocatable :: analysed_U(:,:), analysed_ER(:,:), analysed_SV(:,:)
        real(kind=8), allocatable :: analysed_dU(:,:), analysed_dER(:,:)
        real(kind=8), allocatable :: analysed_d2U(:,:), analysed_d2ER(:,:)
        integer :: ratio, N_times, analysis_time_idx
        integer :: i, i_idx, i_real

        self.do_backward_analysis = do_backward
        ratio = algo_cfg.dt_a_f_ratio
        N_times = 2*ratio + 1

        if (SIZE(input_core_state.measures_) < 9) then
            write (10,'(A)') 'AR3 analysis requires the 9 AR3 core-state measures.'
            write (*,'(A)') 'AR3 analysis requires the 9 AR3 core-state measures.'
            stop
        end if

        if (SIZE(input_core_state.measures_(1).measure_data, 1) .ne. nb_realisations .or. &
            SIZE(input_core_state.measures_(1).measure_data, 2) .ne. N_times .or. &
            SIZE(input_core_state.measures_(5).measure_data, 1) .ne. nb_realisations .or. &
            SIZE(input_core_state.measures_(5).measure_data, 2) .ne. N_times) then
            write (10,'(A)') 'Invalid global forecast-window dimensions in AR3 analysis.'
            write (*,'(A)') 'Invalid global forecast-window dimensions in AR3 analysis.'
            stop
        end if

        if (do_backward) then
            analysis_time_idx = N_times
        else
            analysis_time_idx = ratio + 1
        end if

        !# Start from the forecast state at the actual analysis time.
        do i = 1, SIZE(ana_core_state_slice.measures_)
            do i_idx = 1, SIZE(attributed_models)
                i_real = attributed_models(i_idx) + 1
                ana_core_state_slice.measures_(i).measure_data(i_idx,1,:) = &
                    input_core_state.measures_(i).measure_data(i_real,analysis_time_idx,:)
            end do
        end do

        if (.not. self.sv_analysis() .and. .not. self.mf_analysis()) return

        if (self.mf_analysis()) then
            call self.compute_full_state_AR3( &
                input_core_state.measures_(1).measure_data, &
                algo_cfg,                                  &
                nb_realisations,                           &
                attributed_models,                         &
                B_full                                     &
            )
            call self.analyse_B_AR3( &
                B_full,                                    &
                algo_cfg,                                  &
                nb_realisations,                           &
                attributed_models,                         &
                b_minus, b, b_plus                         &
            )

            do i_idx = 1, SIZE(attributed_models)
                i_real = attributed_models(i_idx) + 1
                if (do_backward) then
                    ana_core_state_slice.measures_(1).measure_data(i_idx,1,:) = b_plus(i_real,:)
                else
                    ana_core_state_slice.measures_(1).measure_data(i_idx,1,:) = b(i_real,:)
                end if
            end do
        else
            allocate(b_minus(nb_realisations, algo_cfg.Nb()), &
                     source=input_core_state.measures_(1).measure_data(:,1,:))
            allocate(b(nb_realisations, algo_cfg.Nb()), &
                     source=input_core_state.measures_(1).measure_data(:,ratio+1,:))
            allocate(b_plus(nb_realisations, algo_cfg.Nb()), &
                     source=input_core_state.measures_(1).measure_data(:,N_times,:))
        end if

        if (self.sv_analysis()) then
            call self.compute_full_state_AR3( &
                input_core_state.measures_(5).measure_data, &
                algo_cfg,                                  &
                nb_realisations,                           &
                attributed_models,                         &
                Z_full                                     &
            )
            call self.analyse_Z_AR3( &
                Z_full,                                    &
                input_core_state,                           &
                algo_cfg,                                  &
                nb_realisations,                           &
                attributed_models,                         &
                pcaU_operator,                             &
                algo_avg_prior,                            &
                b_minus, b, b_plus,                        &
                ana_Z, ana_dZ, ana_d2Z,                    &
                analysed_U, analysed_ER, analysed_SV       &
            )

            call self.dZ_to_dU_dER_AR3( &
                ana_dZ, algo_cfg, pcaU_operator, analysed_dU, analysed_dER)
            call self.dZ_to_dU_dER_AR3( &
                ana_d2Z, algo_cfg, pcaU_operator, analysed_d2U, analysed_d2ER)

            do i_idx = 1, SIZE(attributed_models)
                i_real = attributed_models(i_idx) + 1
                ana_core_state_slice.measures_(2).measure_data(i_idx,1,:) = analysed_U(i_real,:)
                ana_core_state_slice.measures_(3).measure_data(i_idx,1,:) = analysed_SV(i_real,:)
                ana_core_state_slice.measures_(4).measure_data(i_idx,1,:) = analysed_ER(i_real,:)
                ana_core_state_slice.measures_(5).measure_data(i_idx,1,:) = ana_Z(i_real,:)
                ana_core_state_slice.measures_(6).measure_data(i_idx,1,:) = analysed_dU(i_real,:)
                ana_core_state_slice.measures_(7).measure_data(i_idx,1,:) = analysed_d2U(i_real,:)
                ana_core_state_slice.measures_(8).measure_data(i_idx,1,:) = analysed_dER(i_real,:)
                ana_core_state_slice.measures_(9).measure_data(i_idx,1,:) = analysed_d2ER(i_real,:)
            end do

            call self.init_forecast_AR3( &
                ana_Z, ana_dZ, ana_d2Z, algo_cfg, attributed_models, Z_AR3)
        end if

    end subroutine analysis_step_AR3
!==========================================================================================================================

!==========================================================================================================================
    subroutine init_AugkfAnalyserAR(self, config, legendre_polys, nb_realisations, seed)
    !*****************************************************************************************************************
    !"""
    !Class that handles the analyses of the Augmented State Kalman Filter algorithm with DIFF treated as a contribution to ER.
    !"""
    !*****************************************************************************************************************
        class(AugkfAnalyserAR1), intent(inout) :: self
        class(ComputationConfig), intent(in) :: config
        class(legendre_polys_type), intent(in) :: legendre_polys
        integer :: nb_realisations, seed, nb_obs_mf, nb_obs_sv
        character(len=10) :: keys(5) = [ &
                                    "MF", &
                                    "SV", &
                                    "U", &
                                    "ER", &
                                    "Z" ]
        
        call self.init_GenericComputer(config, legendre_polys, nb_realisations, seed)
        
        !# Date-based dicts for observations, stored in measure_observations dict
        call self.extract_observations(nb_obs_mf, nb_obs_sv)
        
        if (nb_obs_mf == 0) then
            write(10, '(A)') "No observation was extracted for MF! Analyses on b will be completely skipped"
            write(*, '(A)') "No observation was extracted for MF! Analyses on b will be completely skipped"
            stop
        else if (nb_obs_sv == 0) then
            write(10, '(A)') "No observation was extracted for SV! Analyses on u, e and d/dt b will be completely skipped"
            write(*, '(A)') "No observation was extracted for SV! Analyses on u, e and d/dt b will be completely skipped"
            stop
        else
            write(10, '(A, i3, A, i3, A)') "Finished extracting the observations: MF (",nb_obs_mf, ") and SV (", nb_obs_sv, ")!"
            write(*, '(A, i3, A, i3, A)') "Finished extracting the observations: MF (",nb_obs_mf, ") and SV (", nb_obs_sv, ")!"
        end if
        
        !# Dict to store the misfits of the current analysis
        self.keys = keys
        self.current_misfits = 0.0d0
        call self.invalid_misfits(keys)
        
        !# Performing or not a backward analysis        
        self.do_backward_analysis = .False.

    end subroutine init_AugkfAnalyserAR
!==========================================================================================================================
    
!==========================================================================================================================
    subroutine invalid_misfits(self, keys)
        class(AugkfAnalyserAR1), intent(inout) :: self
        character(len=*), intent(in) :: keys(:)
        real(kind=8) :: quiet_nan
        integer :: i, idx
        character(len=10) :: possible_keys(5) =[ &
                                    "MF", &
                                    "SV", &
                                    "U", &
                                    "ER", &
                                    "Z" ]
        quiet_nan = ieee_value(quiet_nan, ieee_quiet_nan)
        !self.current_misfits = 0.0d0
        do i = 1, SIZE(keys)
            idx = findloc(possible_keys, TRIM(keys(i)), dim=1)
            self.current_misfits(idx) = quiet_nan
        end do
    end subroutine
!==========================================================================================================================
    
!==========================================================================================================================
    subroutine extract_observations(self, nb_obs_mf, nb_obs_sv)
    !*****************************************************************************************************************
    !"""
    !Extracts the observations for all obs types in the config. Updates the internal dictionaries observations_mf and observations_sv.
    !
    !:return: The numbers of dates for MF and SV for eventual checking.
    !:rtype: int, int
    !"""
    !*****************************************************************************************************************
        class(AugkfAnalyserAR1), intent(inout) :: self
        integer, intent(out) :: nb_obs_mf, nb_obs_sv
        integer :: i
        
        !# building function can be either, build_go_vo_observations, build_covobs_observations,
        !# build_covobs_observations or build_covobs_hdf5_observations
        write(10, '(A,A,A)') "Reading ", trim(self.cfg.obs_type), " data as observations..."
        write(*, '(A,A,A)') "Reading ", trim(self.cfg.obs_type), " data as observations..."
        self.type_(1) = 'SV'
        self.type_(2) = 'MF'
        
        call build_chaos_hdf5_observations(self.cfg, self.algo_nb_realisations, 'SV', self.algo_seed, self.measure_observations_SV)
        call build_chaos_hdf5_observations(self.cfg, self.algo_nb_realisations, 'MF', self.algo_seed, self.measure_observations_MF)
        
        nb_obs_mf = SIZE(self.measure_observations_MF)
        nb_obs_sv = SIZE(self.measure_observations_SV)
    end subroutine extract_observations
!==========================================================================================================================
    
!==========================================================================================================================    
    function sv_analysis(self) result (log)
    !"""
    !Shortcut to check if sv_analysis is performed
    !"""
        class(AugkfAnalyserAR1), intent(in) :: self
        logical :: log
        !# 1 means analysis on sv, 0 means no analysis on sv
        log = all(self.ana_sv)
    end function
!==========================================================================================================================
    
!==========================================================================================================================    
    function mf_analysis(self) result (log)
    !"""
    !Shortcut to check if mf_analysis is performed
    !"""
        class(AugkfAnalyserAR1), intent(in) :: self
        logical :: log
        
        log = all(self.ana_mf)        
    end function
!==========================================================================================================================
    
!==========================================================================================================================
    subroutine check_if_analysis_data(self, algo_config, i_analysis, do_backward)
    !*****************************************************************************************************************
    !"""
    !check if there is mf and/or sv observation at next analysis time (times if AR3)
    !and prepare observation data for either AR1 or AR3 analysis
    !
    !:param i_analysis: analysis time iteration
    !:type i_analysis: int
    !:return: update self.ana_sv and self.ana_mf
    !"""
    !*****************************************************************************************************************
        class(AugkfAnalyserAR1), intent(inout) :: self
        class(ComputationConfig), intent(in) :: algo_config
        integer, intent(in) :: i_analysis
        logical :: do_backward
        integer :: i, Nt, t
        character(len=20) :: AR_type
        integer, allocatable :: times(:)
        
        
        AR_type = algo_config.AR_type
        
        !# set number of times (Nt) involved in analysis depending on AR_type
        if (trim(AR_type) == "AR3") then
            Nt = 3
        else
            Nt = 1
        end if
        
        !# init boolean vectors 
        if (ALLOCATED(self.ana_sv)) deallocate(self.ana_sv)
        if (ALLOCATED(self.ana_mf)) deallocate(self.ana_mf)
        allocate(self.ana_sv(Nt), source=.false.)
        allocate(self.ana_mf(Nt), source=.false.)
        
        if (i_analysis >= algo_config.nb_analyses()) then
            return
        endif 
        
        !# set times vector depending on AR_type
        if (trim(AR_type) == "AR3") then
            allocate(times(3))
            if (do_backward) then
                times = [i_analysis - 1, i_analysis, i_analysis + 1]
            else
                times = [i_analysis, i_analysis + 1, i_analysis + 2]
            end if
        else
            allocate(times(1))
            times (1)= i_analysis+1
        end if
                
        !# loop over Nt
        do i = 1, Nt
            t = times(i)
            !# Observation indices in times are zero-based, while Fortran arrays
            !# are one-based.  The observation container is allocated for every
            !# analysis date, but dates without data have unallocated X/H/Rxx
            !# members.  Check both the index and the members, as the Python
            !# implementation does by testing whether the date key exists.
            if ((t < SIZE(self.measure_observations_SV, 1)) .and. (t>=0)) then
                if (ALLOCATED(self.measure_observations_SV(t+1).X) .and. &
                    ALLOCATED(self.measure_observations_SV(t+1).H) .and. &
                    ALLOCATED(self.measure_observations_SV(t+1).Rxx)) then
                    self.ana_sv(i) = .true.
                end if
            end if
            if ((t < SIZE(self.measure_observations_MF, 1)) .and. (t>=0)) then
                if (ALLOCATED(self.measure_observations_MF(t+1).X) .and. &
                    ALLOCATED(self.measure_observations_MF(t+1).H) .and. &
                    ALLOCATED(self.measure_observations_MF(t+1).Rxx)) then
                    self.ana_mf(i) = .true.
                end if
            end if
        end do
        
        ! #setup self.measure_observations
        if (trim(AR_type) == "AR3") then
            if (self.sv_analysis()) then
                if (ALLOCATED(self.sv_X)) deallocate(self.sv_X)
                if (ALLOCATED(self.sv_H)) deallocate(self.sv_H)
                if (ALLOCATED(self.sv_Rxx)) deallocate(self.sv_Rxx)
                allocate(self.sv_X(3), self.sv_H(3), self.sv_Rxx(3))
                allocate(self.sv_X(1).mat, source=self.measure_observations_SV(times(1)+1).X)
                allocate(self.sv_X(2).mat, source=self.measure_observations_SV(times(2)+1).X)
                allocate(self.sv_X(3).mat, source=self.measure_observations_SV(times(3)+1).X)
                
                allocate(self.sv_H(1).mat, source=self.measure_observations_SV(times(1)+1).H)
                allocate(self.sv_H(2).mat, source=self.measure_observations_SV(times(2)+1).H)
                allocate(self.sv_H(3).mat, source=self.measure_observations_SV(times(3)+1).H)
                
                allocate(self.sv_Rxx(1).mat, source=self.measure_observations_SV(times(1)+1).Rxx)
                allocate(self.sv_Rxx(2).mat, source=self.measure_observations_SV(times(2)+1).Rxx)
                allocate(self.sv_Rxx(3).mat, source=self.measure_observations_SV(times(3)+1).Rxx)
            end if
            if (self.mf_analysis()) then
                if (ALLOCATED(self.mf_X)) deallocate(self.mf_X)
                if (ALLOCATED(self.mf_H)) deallocate(self.mf_H)
                if (ALLOCATED(self.mf_Rxx)) deallocate(self.mf_Rxx)
                allocate(self.mf_X(3), self.mf_H(3), self.mf_Rxx(3))
                allocate(self.mf_X(1).mat, source=self.measure_observations_MF(times(1)+1).X)
                allocate(self.mf_X(2).mat, source=self.measure_observations_MF(times(2)+1).X)
                allocate(self.mf_X(3).mat, source=self.measure_observations_MF(times(3)+1).X)
                
                allocate(self.mf_H(1).mat, source=self.measure_observations_MF(times(1)+1).H)
                allocate(self.mf_H(2).mat, source=self.measure_observations_MF(times(2)+1).H)
                allocate(self.mf_H(3).mat, source=self.measure_observations_MF(times(3)+1).H)
                
                allocate(self.mf_Rxx(1).mat, source=self.measure_observations_MF(times(1)+1).Rxx)
                allocate(self.mf_Rxx(2).mat, source=self.measure_observations_MF(times(2)+1).Rxx)
                allocate(self.mf_Rxx(3).mat, source=self.measure_observations_MF(times(3)+1).Rxx)
            end if
        else
            if (self.sv_analysis()) then
                if (ALLOCATED(self.sv_X)) deallocate(self.sv_X)
                if (ALLOCATED(self.sv_H)) deallocate(self.sv_H)
                if (ALLOCATED(self.sv_Rxx)) deallocate(self.sv_Rxx)
                allocate(self.sv_X(1), self.sv_H(1), self.sv_Rxx(1))
                allocate(self.sv_X(1).mat, source=self.measure_observations_SV(times(1)+1).X)                
                allocate(self.sv_H(1).mat, source=self.measure_observations_SV(times(1)+1).H)                
                allocate(self.sv_Rxx(1).mat, source=self.measure_observations_SV(times(1)+1).Rxx)
            end if
            if (self.mf_analysis()) then
                if (ALLOCATED(self.mf_X)) deallocate(self.mf_X)
                if (ALLOCATED(self.mf_H)) deallocate(self.mf_H)
                if (ALLOCATED(self.mf_Rxx)) deallocate(self.mf_Rxx)
                allocate(self.mf_X(1), self.mf_H(1), self.mf_Rxx(1))
                allocate(self.mf_X(1).mat, source=self.measure_observations_MF(times(1)+1).X)                
                allocate(self.mf_H(1).mat, source=self.measure_observations_MF(times(1)+1).H)              
                allocate(self.mf_Rxx(1).mat, source=self.measure_observations_MF(times(1)+1).Rxx)
            end if
        end if
        
        !# if no mf analysis
        if (.not. self.mf_analysis()) then
            write(10,'(A)') 'Skipping MF analysis'
            write(*,'(A)') 'Skipping MF analysis'
            call self.invalid_misfits(['MF'])
        end if
        
        !# if no sv analysis
        if (.not. self.sv_analysis()) then
            write(10,'(A)') 'Skipping SV analysis'
            write(*,'(A)') 'Skipping SV analysis'
            call self.invalid_misfits(['SV'])
        end if
    end subroutine check_if_analysis_data
!==========================================================================================================================

!========================================================================================================================== 
    subroutine analysis_step(self, input_core_state, algo_cfg, nb_realisations, attributed_models, pcaU_operator, algo_avg_prior, ana_core_state_slice)
    !*****************************************************************************************************************
    !""" Does the analysis at time t on the B and Z=[UE] part of the input_core_state.
    !Updates SV = A(B)U - ER in consequence.
    !
    !:param input_core_state: Core state at time t
    !:type input_core_state: corestates.CoreState (dim: nb_realisations x Ncorestate)
    !:return: the analysed core state
    !:rtype: corestates.CoreState (dim: nb_realisations x Ncorestate)
    !"""
    !*****************************************************************************************************************
        class(AugkfAnalyserAR1), intent(inout) :: self
        class(CoreState_type), intent(in) :: input_core_state
        class(ComputationConfig), intent(in) :: algo_cfg
        class(set_prior_type), intent(in) :: algo_avg_prior
        class(CoreState_type), intent(inout) :: ana_core_state_slice
        integer, intent(in) :: attributed_models(:)
        integer, intent(in) :: nb_realisations
        class(NormedPCAOperator), intent(in) :: pcaU_operator
        class(CoreState_type), allocatable :: ana_core_state
        integer :: i, n_rea, n_t, n_coef, nprocs, rank, ierr, comm, local_idx, global_idx
        
        !# copy core state
        allocate(ana_core_state, source=input_core_state)        
        !call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
        !call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
        
        !# if no analysis
        if (.not. self.sv_analysis() .and. .not. self.mf_analysis()) then 
            do i = 1, size(ana_core_state.measures_, 1)                
                n_rea = SIZE(attributed_models)
                !n_t = SIZE(ana_core_state.measures_(i).measure_data, 2)
                !n_coef = SIZE(ana_core_state.measures_(i).measure_data, 3)
                !deallocate(ana_core_state.measures_(i).measure_data)
                !allocate(ana_core_state.measures_(i).measure_data(n_rea, n_t, n_coef))                
                do local_idx = 1, n_rea

                    global_idx = attributed_models(local_idx) + 1

                    ana_core_state_slice.measures_(i).measure_data(local_idx,:,:) = &
                         input_core_state.measures_(i).measure_data(global_idx,:,:)
                end do
            end do
            
            print *, "No analysis performed at this step, returning the input core state as analysed core state!"
            return
        end if
        
        if (self.mf_analysis()) then
            !# perform MF analysis
            call self.analyse_B(ana_core_state, algo_cfg, nb_realisations, attributed_models)
        end if
        
        if (self.sv_analysis()) then
            !# perform SV analysis
            call self.analyse_Z(ana_core_state, algo_cfg, nb_realisations, attributed_models, pcaU_operator, algo_avg_prior)
        end if
        
        ! # ana_core_state[self.algo.attributed_models]        
        do i = 1, size(ana_core_state.measures_, 1)                
            n_rea = SIZE(attributed_models)
            !n_t = SIZE(ana_core_state.measures_(i).measure_data, 2)
            !n_coef = SIZE(ana_core_state.measures_(i).measure_data, 3)
            !deallocate(ana_core_state.measures_(i).measure_data)
            !allocate(ana_core_state.measures_(i).measure_data(n_rea, n_t, n_coef))                
            do local_idx = 1, n_rea

                global_idx = attributed_models(local_idx) + 1

                ana_core_state_slice.measures_(i).measure_data(local_idx,:,:) = &
                        ana_core_state.measures_(i).measure_data(global_idx,:,:)

            end do
        end do
        
    end subroutine
!==========================================================================================================================
    
!========================================================================================================================== 
    subroutine analyse_B(self, inout_core_state, algo_cfg, nb_realisations, attributed_models)
    !*****************************************************************************************************************
    !"""
    !Returns the analysed data for B by a BLUE given the observations.
    !
    !:param input_core_state: NumPy array containing the coefficient data of B
    !:type input_core_state: np.array (dim: nb_realisations x Nb)
    !:param mf_X: Observation data to use for the BLUE
    !:type mf_X: Observation
    !:param mf_H: Observation matrix to use for the BLUE
    !:type mf_H: Observation
    !:param mf_Rxx: Observation error to use for the BLUE
    !:type mf_Rxx: Observation
    !:return: NumPy array containing the analysed coefficient data of B
    !:rtype: np.array (dim: nb_realisations x Nb)
    !"""
    !*****************************************************************************************************************
        class(AugkfAnalyserAR1), intent(inout) :: self
        class(CoreState_type), intent(inout) :: inout_core_state
        class(ComputationConfig), intent(in) :: algo_cfg
        integer, intent(in) :: nb_realisations
        integer, intent(in) :: attributed_models(:)
        REAL(kind=8), allocatable :: analysed_B(:, :)
        class(measure_observations_mat), allocatable :: mf_X(:), Hb(:), Rbb(:)
        REAL(KIND=8), allocatable :: Pbb_forecast(:,:), Kbb(:,:)
        REAL(KIND=8), allocatable :: P_eig_val(:), P_eig_vec(:,:), HX_b(:,:), inv_Rbb(:,:)
        integer :: i_real, i, info, rank, ierr
        
        allocate(mf_X, source=self.mf_X)
        !# obs operator
        allocate(Hb, source=self.mf_H)
        !# obs error
        allocate(Rbb, source=self.mf_Rxx)
        !# compute Pbb from B state
        call self.remove_small_correlations(inout_core_state.measures_(1).measure_data(:,1,:), 1.0d-10, algo_cfg, Pbb_forecast)
        !# Updates the B part of the core_state by the result of the Kalman filter for each model
        write(*,*) "Getting best linear unbiased estimate of B..."
        write(10,*) "Getting best linear unbiased estimate of B..."
        allocate(analysed_B(nb_realisations, algo_cfg.Nb()))
        analysed_B = 0.0d0
        
        if (TRIM(algo_cfg.kalman_norm) == 'l2') then  !  # for non least square norm, iteration are needed
            call compute_Kalman_gain_matrix(Pbb_forecast, Hb(1).mat, Rbb(1).mat, .True., Kbb)
            
            do i = 1, SIZE(attributed_models)
                i_real = attributed_models(i) + 1
                !analysed_B(i_real, :) = input_core_state.measures_(1).measure_data(i_real, 1, :) + matmul(Kbb, (mf_X(1).mat(:, i) - matmul(Hb(1).mat, input_core_state.measures_(1).measure_data(i_real, 1, :))))
                call get_BLUE(inout_core_state.measures_(1).measure_data(i_real,1,:), &
                                mf_X(1).mat(i_real,:), &
                                Pbb_forecast, Hb(1).mat, &
                                Rbb(1).mat, &
                                Kbb, &
                                .True., &
                                analysed_B(i_real,:))
            end do
            
        else if (TRIM(algo_cfg.kalman_norm) == 'huber') then
            !# compute inverse of P_bb before loop on reals using its symmetry
            allocate(P_eig_val(SIZE(Pbb_forecast, 1)))
            allocate(P_eig_vec, source=Pbb_forecast)
            call syevd(P_eig_vec, P_eig_val, 'V', 'U', info)
            P_eig_val = max(P_eig_val, 1.0d-10)
        
            do i = 1, SIZE(attributed_models)
                i_real = attributed_models(i) + 1
                call compute_Kalman_huber(inout_core_state.measures_(1).measure_data(i_real,1,:), &
                                    mf_X(1).mat(i_real,:), &
                                    P_eig_val, &
                                    P_eig_vec, &
                                    Hb(1).mat, &
                                    Rbb(1).mat, &
                                    50, &
                                    analysed_B(i_real,:))
            end do
            
        else
            write(*,*) "Invalid value of param kalman_norm, should be equal to huber or l2."
            write(10,*) "Invalid value of param kalman_norm, should be equal to huber or l2."
        end if            
        
        call MPI_ALLREDUCE( MPI_IN_PLACE, &
                            analysed_B, &
                            SIZE(analysed_B), &
                            MPI_DOUBLE_PRECISION, &
                            MPI_SUM, &
                            MPI_COMM_WORLD, &
                            ierr)
        
        ! # Compute the misfits for B (Y - HX)
        allocate(HX_b, source=TRANSPOSE(MATMUL(Hb(1).mat, TRANSPOSE(analysed_B))))
        call max_inv(Rbb(1).mat, inv_Rbb)
        call compute_misfit(mf_X(1).mat, HX_b, inv_Rbb, self.current_misfits(1))
        inout_core_state.measures_(1).measure_data(:,1,:) = analysed_B
        !print *, 'imput', inout_core_state.measures_(1).measure_data(1,1,2)
    end subroutine
!==========================================================================================================================
    
!========================================================================================================================== 
    subroutine setup_Hz(self, Ab, Nu, algo_cfg, Hz)
    !*****************************************************************************************************************
    !"""
    !Compute the matrix [Ab | I_e] where Ab is the contains the Gaunt elasser integrals,
    !while I_e is the identity matrix of size Nsv
    !    
    !:param Nu: dimension of the flow
    !:type Nu: int
    !:return: matrix Nsv x (Nu + Ne)
    !:rtype: numpy array
    !"""
    !*****************************************************************************************************************
        class(AugkfAnalyserAR1), intent(in) :: self
        real(kind=8), intent(in) :: Ab(:,:)
        class(ComputationConfig), intent(in) :: algo_cfg
        integer, intent(in) :: Nu
        REAL(KIND=8), allocatable, intent(out) :: Hz(:, :)
        REAL(KIND=8), allocatable :: I_e(:,:)
        integer :: Nsv, Ne, i
        
        if (SIZE(Ab, 2) .ne. Nu) then
            write(*,*) "mismatch in dimension of Ab and flow"
            write(10,*) "mismatch in dimension of Ab and flow"
            stop
        end if
        
        Nsv = algo_cfg.Nsv()
        allocate(Hz(Nsv, Nu + Nsv))
        Hz = 0.0d0
        
        !# Set the observation operator to A(B) for U
        Hz(1:Nsv, 1:Nu) = Ab
        !# and identity for E
        allocate(I_e(Nsv, Nsv))
        I_e = 0.0d0
        do i = 1, Nsv
            I_e(i, i) = algo_cfg.compute_e
        end do
        
        Hz(1:Nsv, Nu+1:Nu+Nsv) = I_e
        
    end subroutine
!==========================================================================================================================
    
!========================================================================================================================== 
    subroutine remove_small_correlations(self, input_core_state_matrix, eps, algo_cfg, result_)
    !*****************************************************************************************************************
    !"""
    !Apply the graphical lasso to the correlation matrix. The correlation matrix is computed from
    !the covariance matrix, either Pzz or Pbb in practice, with C[i, j] = P[i, j] / (P[i, i] P[j, j]).
    !Warning: In some cases, some variance elements can be zero, for instance if the initialisation
    !parameter, core_state_init, is set to constant. Then the correlation matrix cannot be computed
    !and the Glasso is not applied.
    !
    !If the glasso parameter, self.cfg.remove_spurious, is set to 0 (np.inf), then the resp. diagonal (empirical)
    !covariance matrix is returned.
    !Otherwise the glasso is applied on the correlation matrix.
    !
    !:param input_core_state: Corestate which can either be Z or B,
    !                            at a given time for all realizations (dim: nb_realisations x Ncorestate, matrix)
    !:param eps: threshold that determines if a value should be considered as null. During the division to
    !            get the correlation matrix, null values are replaced by eps
    !:type eps: float
    !"""
    !*****************************************************************************************************************
        class(AugkfAnalyserAR1), intent(in) :: self
        real(kind=8), intent(in) :: input_core_state_matrix(:,:)
        real(kind=8), intent(in) :: eps
        class(ComputationConfig), intent(in) :: algo_cfg
        real(kind=8), allocatable, intent(out) :: result_(:, :)
        real(kind=8), allocatable :: diag_A(:)
        real(kind=8), allocatable :: P_forecast(:,:), D_(:,:), C_forecast(:,:)
        integer :: n, maxIt, msg, warm, info
        real(kind=8) :: thr
        real(kind=8), allocatable :: L(:,:), X(:,:), C_lasso(:,:), Wd(:), WXj(:), D_2(:,:)
        integer :: i
        
        ! # computation of the empirical Pzz_forecast
        call cov(input_core_state_matrix, P_forecast)
        
        if (algo_cfg.remove_spurious < eps) then
            allocate(result_, source=P_forecast)
            return
        end if
        
        if (.not. ieee_is_finite(algo_cfg.remove_spurious)) then
            allocate(result_, source=P_forecast)
            result_ = 0.0d0
            
            do i = 1, SIZE(result_, 1)
                result_(i, i) = P_forecast(i, i)
            end do
            
            return
        end if
        
        allocate(diag_A(SIZE(P_forecast, 1)))
        do i = 1, SIZE(P_forecast, 1)
            diag_A(i) = P_forecast(i, i)
        end do
        
        if (any(abs(diag_A) < eps)) then
            !# avoid division by zeros (exactly zeros, small numbers are left) by regularization
            !# with many zeros, might give a hard time to the graphical lasso algo, as result may not converge.
            write(10,*) 'Some coefficients in the diagonal sample covariance matrix are very close to zero'
            write(*,*) 'Some coefficients in the diagonal sample covariance matrix are very close to zero'
            !$omp parallel do default(shared) private(i)
            do i = 1, SIZE(P_forecast, 1)
                if (abs(diag_A(i)) < eps) then                    
                    P_forecast(i, i) = eps
                end if
            end do
            !$omp end parallel do
        end if
        
        
        !# if some values in the diagonal of the covariance matrix are zero, it will still give a 1 in the diagonal of the correlation matrix
        !# compute correlation matrix
        call diag_sq_inv(P_forecast, D_)
        allocate(C_forecast(SIZE(P_forecast, 1), SIZE(P_forecast, 2)))
        C_forecast = matmul(matmul(D_, P_forecast), D_)
        
        ! # Compute the lasso approximation
        n = SIZE(C_forecast, 1)
        maxIt = 100
        msg = 0
        warm = 0
        thr = 1.0d-5
        
        allocate(L(n,n), X(n,n), C_lasso(n,n), Wd(n), WXj(n))
        
        L = algo_cfg.remove_spurious
    
        do i = 1, n
            L(i, i) = 0.0d0
        end do
        call glassofast(n,C_forecast,L,thr,maxIt,msg,warm,X,C_lasso,Wd,WXj,info)
        
        !# compute the P_lasso from P_forecast
        call diag_sq(P_forecast, D_2)
        result_ = matmul(matmul(D_2, C_lasso), D_2)
        
    end subroutine
!==========================================================================================================================
    
!========================================================================================================================== 
    subroutine analyse_Z(self, inout_core_state, algo_cfg, nb_realisations, attributed_models, pcaU_operator, algo_avg_prior)
    !*****************************************************************************************************************
    !"""
    !Returns the analysed data for the augmented state Z = [U ER] and SV by a BLUE given the observations.
    !
    !:param input_core_state: 2D CoreState containing the coefficient data
    !:type input_core_state: CoreState
    !:param sv_X: Observation data to use for the BLUE
    !:type sv_X: Observation
    !:param sv_H: Observation matrix to use for the BLUE
    !:type sv_H: Observation
    !:param sv_Rxx: Observation error to use for the BLUE
    !:type sv_Rxx: Observation
    !:return: 2D analysed Z U ER SV
    !:rtype: 2D arrays (Nreal x Ncoef)
    !
    !"""
    !*****************************************************************************************************************
        class(AugkfAnalyserAR1), intent(inout) :: self
        class(CoreState_type), intent(inout) :: inout_core_state
        class(ComputationConfig), intent(in) :: algo_cfg
        class(set_prior_type), intent(in) :: algo_avg_prior
        integer, intent(in) :: nb_realisations
        integer, intent(in) :: attributed_models(:)
        REAL(kind=8), allocatable :: analysed_Z(:, :), analysed_SV(:, :), analysed_ER(:, :), analysed_U(:, :)
        REAL(kind=8), allocatable :: Pzz_forecast(:, :)
        REAL(kind=8), allocatable :: sv_X(:,:), sv_H(:,:), sv_Rxx(:,:), sv_X_real(:,:)
        REAL(kind=8), allocatable :: Ab(:,:), S_u(:,:), setup_Hz_mat(:,:)
        REAL(kind=8), allocatable :: complete_H(:,:), PzzHT(:,:), HPzzHT(:,:), HX_z(:,:), inv_Rzz(:,:)
        integer :: i_real, i_idx, info, rank, ierr
        type(input_core_state_type) :: CoreState_temp
        class(NormedPCAOperator), intent(in) :: pcaU_operator
        
        
        allocate(sv_X, source=self.sv_X(1).mat)
        allocate(sv_H, source=self.sv_H(1).mat)
        allocate(sv_Rxx, source=self.sv_RXX(1).mat)
        
        !# compute necessary matrices for Kalman filter
        call self.remove_small_correlations(inout_core_state.measures_(5).measure_data(:,1,:), 1.0d-10, algo_cfg, Pzz_forecast)
        
        !initialize analysed_Z,analysed_SV,analysed_ER,analysed_U
        allocate(analysed_Z(nb_realisations, algo_cfg.Nz()))
        analysed_Z = 0.0d0
        allocate(analysed_SV(nb_realisations, algo_cfg.Nsv()))
        analysed_SV = 0.0d0
        allocate(analysed_ER(nb_realisations, algo_cfg.Nsv()))
        analysed_ER = 0.0d0
        allocate(analysed_U(nb_realisations, algo_cfg.Nu2()))
        analysed_U = 0.0d0
        allocate(sv_X_real, source=sv_X)
        sv_X_real = 0.0d0
        allocate(HX_z(nb_realisations, SIZE(sv_X, 2)), source=0.0d0)
        
        
        do i_idx = 1, SIZE(attributed_models)
            i_real = attributed_models(i_idx) + 1
            
            !# Compute A(b)
            CoreState_temp.Lsv = inout_core_state.cs_Lsv()
            CoreState_temp.Lu = inout_core_state.cs_Lu()
            CoreState_temp.Lb = inout_core_state.cs_Lb()
            CoreState_temp.Nsv = inout_core_state.cs_Nsv()
            CoreState_temp.Nu2 = inout_core_state.cs_Nu2()
            CoreState_temp.Nb = inout_core_state.cs_Nb()
            if (.not. ALLOCATED(CoreState_temp.B)) allocate(CoreState_temp.B, source=inout_core_state.measures_(1).measure_data(i_real, 1, :))
            CoreState_temp.B = inout_core_state.measures_(1).measure_data(i_real, 1, :)
            call self.compute_Ab(CoreState_temp, Ab)
            
            
            !# The complete H operator (No x Ncoefs of U) is:
            !# Hsv (No x (Ncoefs of core + Ncoefs SV)) * Hz ((Ncoefs of core + Ncoefs SV) x Ncoefs of U)
            if (algo_cfg.pca == 1) then
                !# if PCA
                call pcaU_operator.S_u(S_u)
                call self.setup_Hz(MATMUL(Ab, S_u), algo_cfg.N_pca_u, algo_cfg, setup_Hz_mat)
                if (.not. ALLOCATED(complete_H)) allocate(complete_H, source=MATMUL(sv_H, setup_Hz_mat))
                complete_H = MATMUL(sv_H, setup_Hz_mat)
            else                
                !# if no PCA
                call self.setup_Hz(Ab, algo_cfg.Nu2(), algo_cfg, setup_Hz_mat)
                if (.not. ALLOCATED(complete_H)) allocate(complete_H, source=MATMUL(sv_H, setup_Hz_mat))
                complete_H = MATMUL(sv_H, setup_Hz_mat)
            end if
            
            if (.not. ALLOCATED(PzzHT)) allocate(PzzHT, source=MATMUL(Pzz_forecast, TRANSPOSE(complete_H)))
            PzzHT = MATMUL(Pzz_forecast, TRANSPOSE(complete_H))
            if (.not. ALLOCATED(HPzzHT)) allocate(HPzzHT, source=MATMUL(complete_H, PzzHT))
            HPzzHT = MATMUL(complete_H, PzzHT)
            
            !# Z is centered on 0 so we must remove the mean from the observation data
            !# Y = Ab (U+U0) + (ER + ER0) => Y - Ab U0 - ER0 = Ab U + ER
            sv_X_real(i_real, :) = sv_X(i_real, :) - (MATMUL(sv_H, MATMUL(Ab, RESHAPE(algo_avg_prior.U, [SIZE(algo_avg_prior.U)]))) + MATMUL(sv_H, RESHAPE(algo_avg_prior.ER, [SIZE(algo_avg_prior.ER)])))
            call compute_Kalman_huber_parameter_basis(inout_core_state.measures_(5).measure_data(i_real,1,:), &
                                                        sv_X_real(i_real, :), &
                                                        HPzzHT, PzzHT, &
                                                        complete_H, &
                                                        sv_Rxx, &
                                                        'huber', 50, &
                                                        analysed_Z(i_real,:))
            !# complete_H depends on the magnetic field of this realization.
            !# Save its own predicted observations for the global SV misfit.
            HX_z(i_real,:) = MATMUL(complete_H, analysed_Z(i_real,:))
            call Z_to_U_ER1(algo_cfg, algo_avg_prior, pcaU_operator, analysed_Z(i_real,:), analysed_U(i_real,:), analysed_ER(i_real,:))
            analysed_SV(i_real,:) = MATMUL(Ab, analysed_U(i_real,:)) + analysed_ER(i_real,:)            
        end do
        
        call MPI_ALLREDUCE( MPI_IN_PLACE, &
                            analysed_Z, &
                            SIZE(analysed_Z), &
                            MPI_DOUBLE_PRECISION, &
                            MPI_SUM, &
                            MPI_COMM_WORLD, &
                            ierr)
        call MPI_ALLREDUCE( MPI_IN_PLACE, &
                            analysed_U, &
                            SIZE(analysed_U), &
                            MPI_DOUBLE_PRECISION, &
                            MPI_SUM, &
                            MPI_COMM_WORLD, &
                            ierr)
        call MPI_ALLREDUCE( MPI_IN_PLACE, &
                            analysed_SV, &
                            SIZE(analysed_SV), &
                            MPI_DOUBLE_PRECISION, &
                            MPI_SUM, &
                            MPI_COMM_WORLD, &
                            ierr)
        call MPI_ALLREDUCE( MPI_IN_PLACE, &
                            analysed_ER, &
                            SIZE(analysed_ER), &
                            MPI_DOUBLE_PRECISION, &
                            MPI_SUM, &
                            MPI_COMM_WORLD, &
                            ierr)
        call MPI_ALLREDUCE( MPI_IN_PLACE, &
                            sv_X_real, &
                            SIZE(sv_X_real), &
                            MPI_DOUBLE_PRECISION, &
                            MPI_SUM, &
                            MPI_COMM_WORLD, &
                            ierr)
        call MPI_ALLREDUCE( MPI_IN_PLACE, &
                            HX_z, &
                            SIZE(HX_z), &
                            MPI_DOUBLE_PRECISION, &
                            MPI_SUM, &
                            MPI_COMM_WORLD, &
                            ierr)
        
        ! # Compute the misfits for SV (Y - HX)
        call max_inv(sv_Rxx, inv_Rzz)
        call compute_misfit(sv_X_real, HX_z, inv_Rzz, self.current_misfits(2))
        
        inout_core_state.measures_(2).measure_data(:,1,:) = analysed_U
        inout_core_state.measures_(3).measure_data(:,1,:) = analysed_SV
        inout_core_state.measures_(4).measure_data(:,1,:) = analysed_ER
        inout_core_state.measures_(5).measure_data(:,1,:) = analysed_Z        
    end subroutine
!==========================================================================================================================
end module
