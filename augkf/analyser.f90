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
        procedure :: sv_analysis, mf_analysis, analysis_step, analyse_B
        procedure :: remove_small_correlations
    end type AugkfAnalyserAR1
    
    type, extends(AugkfAnalyserAR1),public :: AugkfAnalyserAR3
        
    contains
        !procedure :: init_AugkfForecasterAR3
        !procedure :: forecast_step
    end type AugkfAnalyserAR3
    
contains    
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
            !# if observation found in Fortran  idx:1-7; in Python idx:0-6 
            if ((t < SIZE(self.measure_observations_SV, 1)) .and. (t>=0)) then
                !#update self.ana_sv
                self.ana_sv(i) = .true.
            end if
            if ((t < SIZE(self.measure_observations_MF, 1)) .and. (t>=0)) then
                !#update self.ana_mf
                self.ana_mf(i) = .true.
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
    subroutine analysis_step(self, input_core_state, algo_cfg, nb_realisations, attributed_models)
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
        class(AugkfAnalyserAR1), intent(in) :: self
        class(CoreState_type), intent(in) :: input_core_state
        class(ComputationConfig), intent(in) :: algo_cfg
        integer, intent(in) :: attributed_models(:)
        integer, intent(in) :: nb_realisations
        class(CoreState_type), allocatable :: ana_core_state
        integer :: i, n_rea, n_t, n_coef, nprocs, rank, ierr, comm, local_idx, global_idx
        
        !# copy core state
        allocate(ana_core_state, source=input_core_state)        
        call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
        
        !# if no analysis
        if (.not. self.sv_analysis() .and. .not. self.mf_analysis()) then 
            do i = 1, size(ana_core_state.measures_, 1)
                deallocate(ana_core_state.measures_(i).measure_data)
                n_rea = SIZE(attributed_models)
                n_t = SIZE(ana_core_state.measures_(i).measure_data, 2)
                n_coef = SIZE(ana_core_state.measures_(i).measure_data, 3)
                allocate(ana_core_state.measures_(i).measure_data(n_rea, n_t, n_coef))
                
                do local_idx = 1, n_rea

                    global_idx = attributed_models(local_idx) + 1

                    ana_core_state.measures_(i).measure_data(local_idx,:,:) = &
                         input_core_state.measures_(i).measure_data(global_idx,:,:)

                end do
            end do
            
            print *, "No analysis performed at this step, returning the input core state as analysed core state!"
        end if
        
        if (self.mf_analysis()) then
            !# perform MF analysis
            call self.analyse_B(ana_core_state, algo_cfg, nb_realisations, attributed_models)
        end if
        
        
        !print *, "rank", rank, ana_core_state.measures_(3).measure_data(1,1,1), input_core_state.measures_(3).measure_data(1,1,1)
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
        class(AugkfAnalyserAR1), intent(in) :: self
        class(CoreState_type), intent(inout) :: inout_core_state
        class(ComputationConfig), intent(in) :: algo_cfg
        integer, intent(in) :: nb_realisations
        integer, intent(in) :: attributed_models(:)
        REAL(kind=8), allocatable :: analysed_B(:, :)
        class(measure_observations_mat), allocatable :: mf_X(:), Hb(:), Rbb(:)
        REAL(KIND=8), allocatable :: Pbb_forecast(:,:), Kbb(:,:)
        REAL(KIND=8), allocatable :: P_eig_val(:), P_eig_vec(:,:)
        integer :: i_real, i, info
        
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
        
        !print *, P_eig_vec(2,:)
        
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
        do i = 1, SIZE(result_, 1)
            diag_A(i) = P_forecast(i, i)
        end do
        
        if (any(abs(diag_A) < eps)) then
            !# avoid division by zeros (exactly zeros, small numbers are left) by regularization
            !# with many zeros, might give a hard time to the graphical lasso algo, as result may not converge.
            write(10,*) 'Some coefficients in the diagonal sample covariance matrix are very close to zero'
            write(*,*) 'Some coefficients in the diagonal sample covariance matrix are very close to zero'
            !$omp parallel do default(shared) private(i)
            do i = 1, SIZE(result_, 1)
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
end module