module analyser
    use mpi
    use common
    use computer
    use observations
    use, intrinsic :: ieee_arithmetic
    use config
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
        procedure :: sv_analysis, mf_analysis
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

    
end module