module analyser
    use mpi
    use common
    use computer
    implicit none
    
    
    integer, private :: comm
    integer, private :: rank
    
    type, extends(GenericComputer),public :: AugkfAnalyserAR1
    
    contains
        procedure :: init_AugkfAnalyserAR
        procedure :: extract_observations
    end type AugkfAnalyserAR1
    
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
    subroutine init_AugkfAnalyserAR(self, config, legendre_polys)
    !*****************************************************************************************************************
    !"""
    !Class that handles the analyses of the Augmented State Kalman Filter algorithm with DIFF treated as a contribution to ER.
    !"""
    !*****************************************************************************************************************
        class(AugkfAnalyserAR1), intent(inout) :: self
        class(ComputationConfig), intent(in) :: config
        class(legendre_polys_type), intent(in) :: legendre_polys
        
        call self.init_GenericComputer(config, legendre_polys)
        
        !# Bool to deactivate checks on AR processes
        
    end subroutine init_AugkfAnalyserAR
!==========================================================================================================================
    
!==========================================================================================================================
    subroutine extract_observations(self)
    !*****************************************************************************************************************
    !"""
    !Extracts the observations for all obs types in the config. Updates the internal dictionaries observations_mf and observations_sv.
    !
    !:return: The numbers of dates for MF and SV for eventual checking.
    !:rtype: int, int
    !"""
    !*****************************************************************************************************************
        class(AugkfAnalyserAR1), intent(in) :: self
        
    end subroutine extract_observations
!==========================================================================================================================

    
end module