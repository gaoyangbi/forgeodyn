module generic_algo
    use utilities
    use config
    implicit none
    
    type, abstract :: GenericAlgo
        class(ComputationConfig), allocatable :: config
        integer :: nb_realisations
        
    contains
        procedure :: init_GenericAlgo
        procedure(init_corestates_iface), deferred :: init_corestates
        procedure(analysis_step_iface),  deferred :: analysis_step_algo
        procedure(forecast_step_iface),  deferred :: forecast_step_algo
        procedure :: get_current_misfits
    end type 
    
    abstract interface
        subroutine init_corestates_iface(self, random_state, computed_states, forecast_states, analysed_states, misfits, Z_AR3)
            use corestate
            import :: GenericAlgo
            class(GenericAlgo), intent(inout) :: self
            real(kind=8), intent(in) :: random_state
            type(CoreState_type), intent(out) :: computed_states, forecast_states, analysed_states, misfits
            REAL(kind=8), allocatable, intent(out) :: Z_AR3(:,:,:)
        end subroutine

        subroutine analysis_step_iface(self, input_core_state, analysis_time)
            import :: GenericAlgo
            class(GenericAlgo), intent(inout) :: self
            real(8), intent(in) :: input_core_state
            real(8), intent(in) :: analysis_time
        end subroutine

        subroutine forecast_step_iface(self, input_core_state, random_state)
            import :: GenericAlgo
            class(GenericAlgo), intent(inout) :: self
            real(8), intent(in) :: input_core_state
            integer, intent(in), optional :: random_state
        end subroutine
    end interface
    
    
contains
!==========================================================================================================================    
    subroutine init_GenericAlgo(self, config, nb_realisations)
    !*****************************************************************************************************************
    !"""
    !Base class that defines the interface of Algo for it to be used in run.py
    !"""
    !*****************************************************************************************************************
        class(GenericAlgo), intent(inout) :: self
        class(ComputationConfig), intent(in) :: config
        integer, intent(in) :: nb_realisations
        
        allocate(self.config)
        self.config = config
        
        if (nb_realisations <= 0) then
            write(10,'(A)') "Number of realisation should be a positive int ! Got {} instead."
            stop
        end if
        self.nb_realisations = nb_realisations
    end subroutine    
!==========================================================================================================================
    
!==========================================================================================================================
    function get_current_misfits(this, measure) result(res)
        class(GenericAlgo), intent(in) :: this
        integer, intent(in) :: measure
        real(8) :: res

        res = 0.0d0
    end function get_current_misfits
!==========================================================================================================================
    
    
end module
    
    