module corestate
    use common
    use utilities
    use reads
    implicit none
    
    type :: CoreState_
        !*****************************************************************************************************************
        !"""
        !CoreState object. It stores the measure datas in a dict of ndarrays.
        !"""
        !*****************************************************************************************************************
        
        integer :: b
    contains
        procedure :: init_CoreState
    end type
    
    
contains
!==========================================================================================================================  
    subroutine init_CoreState(self)
    !*****************************************************************************************************************
    !"""
    !Initiates the CoreState. Initial measures can be given in the init_measures arg:
    !    - Either with only the data (max_degree will be inferred from the data)
    !        Ex: CoreState({SV: np.zeros(224), ...})
    !    - Or by giving also the max_degree as a second member of a 2-tuple/list (the first being the data):
    !        Ex: CoreState({SV: [np.zeros(224), 14], ...})
    !
    !:param init_measures: The measures to add in a dict with members can be np.ndarray or 2-tuple/lists.
    !:type init_measures: dict or None
    !"""
    !*****************************************************************************************************************
        class(CoreState_), intent(inout) :: self        
        
        
    end subroutine    
!==========================================================================================================================
    
!==========================================================================================================================  
    !subroutine addMeasure(self, meas_id, meas_data)
    !!*****************************************************************************************************************
    !!"""
    !!Adds a measure to the CoreState.
    !!
    !!:param meas_id: name of the measure. Used as key of dict for internal storing.
    !!:type meas_id: str
    !!:param meas_data: data of the measure.
    !!:type meas_data: np.ndarray or list
    !!:param meas_max_degree: Max degree of the measure. If None (default), it will be inferred from the last dimension of the data.
    !!:type meas_max_degree: int or None
    !!"""
    !!*****************************************************************************************************************
    !    class(CoreState_), intent(inout) :: self        
    !    
    !    
    !end subroutine    
!==========================================================================================================================
    
end module
    