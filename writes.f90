module writes
    use utilities   
    use HDF5
    use iso_fortran_env
    implicit none
    
    type :: Hdf5GroupHandler
        type(key_measures), allocatable :: measures_(:)
        type(key_max_degrees), allocatable :: max_degrees_(:)
    contains
        !procedure :: init_Hdf5GroupHandler
        !procedure :: update_all
        !procedure :: close_
    end type
    
contains

    
    
    
end module