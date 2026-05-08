module analyser
    use mpi
    use common
    use computer
    use observations
    use, intrinsic :: ieee_arithmetic
    implicit none
    
    
    integer, private :: comm
    integer, private :: rank
    
    type, extends(GenericComputer),public :: AugkfAnalyserAR1
        character(len=10) :: type_(2)
        class(Observation), allocatable :: measure_observations_SV(:)
        class(Observation), allocatable :: measure_observations_MF(:)
        real(kind=8) :: current_misfits(5)
        character(len=10) :: keys(5) 
        logical :: do_backward_analysis
    contains
        procedure :: init_AugkfAnalyserAR, invalid_misfits
        procedure :: extract_observations
    end type AugkfAnalyserAR1
    
    type, extends(AugkfAnalyserAR1),public :: AugkfAnalyserAR3
        
    contains
        !procedure :: init_AugkfForecasterAR3
        !procedure :: forecast_step
    end type AugkfAnalyserAR3
    
contains
!==========================================================================================================================
    subroutine init_MPI(comm, rank)
        integer, intent(out) :: comm, rank
        integer :: nb_proc, ierr
        comm = MPI_COMM_WORLD
        call MPI_Comm_size(comm, nb_proc, ierr)
        call MPI_Comm_rank(comm, rank, ierr)    
    end subroutine
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

    
end module