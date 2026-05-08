module forecaster
    use mpi
    use computer
    implicit none
    
    
    integer, private :: comm
    integer, private :: rank
    
    type, extends(GenericComputer),public :: AugkfForecasterAR1
        logical :: Cholesky_AR_check
    contains
        procedure :: init_AugkfForecasterAR
        !procedure :: forecast_step
    end type AugkfForecasterAR1
    
    type, extends(AugkfForecasterAR1),public :: AugkfForecasterAR3
        
    contains
        !procedure :: init_AugkfForecasterAR3
        !procedure :: forecast_step
    end type AugkfForecasterAR3
    
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
    !subroutine forecast_step(self, input_core_state, Z_AR, seed, i_real, i_t)
    !!*****************************************************************************************************************
    !!"""
    !!Forecasts the input_core_state using AR processes for Z, computation of SV and Euler scheme for B.
    !!
    !!:param input_core_state: core_state of a single realisation at a single date 
    !!:type input_core_state: corestates.CoreState
    !!:param Z_AR: forecast state AR
    !!:type Z_AR: np.array 1D if AR1 (Ncoef) or 2D if AR3 (3 x Ncoef)
    !!:param seed: random seed
    !!:type seed: int
    !!:param i_real: model realisation index
    !!:type i_real: int
    !!:param i_t: time index
    !!:type i_t: int
    !!:return: CoreState containing the result from the forecast
    !!:rtype: corestates.CoreState
    !!"""
    !!*****************************************************************************************************************
    !    class(AugkfForecasterAR1), intent(in) :: self
    !    class(input_core_state_type), intent(in) :: input_core_state
    !    call self.init_GenericComputer(algo)
    !    
    !    !# Bool to deactivate checks on AR processes
    !    self.Cholesky_AR_check = .false.
    !end subroutine forecast_step
!==========================================================================================================================

    
end module