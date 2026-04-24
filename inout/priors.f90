    !*****************************************************************************************************************
    !""" Contains the reading methods for prior types """
    !*****************************************************************************************************************
module priors
    use utilities
    use hdf5
    implicit none
    
    
contains

!========================================================================================================================== 
    subroutine read_0_path(data_directory, prior_type, dt_sampling, measures, prior_path_data)
    !*****************************************************************************************************************
    !Reads prior data under the 0 path format.
    !
    !:param data_directory: directory where the prior data is stored
    !:type data_directory: str
    !:param prior_type: type of the prior data
    !:type prior_type: str
    !:param dt_sampling: smoothing window duration in second
    !:type dt_sampling: float
    !:param measures: sequence of string saying which measures should be extracted. Allowed
    !          values are 'times', 'B', 'U' and 'ER', order do not matter.
    !:return: MF, U, ER, times, dt_samp, tag
    !:rtype: lists
    !*****************************************************************************************************************
        character(len=*), intent(in) :: data_directory, prior_type
        real(kind=8), intent(in) :: dt_sampling
        character(len=*), intent(in) :: measures(:)
        type(prior_data), intent(inout) :: prior_path_data(:)
        character(len=200) :: path
        integer :: i, j
        logical :: valid
        character(len=20) :: allowed_measures(4) = ['times', 'B', 'U', 'ER']
        ! examine the measure---------------------------------------------------
        valid = .false.
        do i = 1, size(measures)
            if (any(measures(i) == allowed_measures)) valid = .true.
            
            if (.not. valid) then
                print *, "Invalid measure given in measures: ", trim(measures(i))
                print *, "Allowed values are: ", "'times', 'B', 'U', 'ER'"
                stop
            end if
            valid = .false.
        end do
        !-----------------------------------------------------------------------
        
        path = trim(data_directory) // '/' // 'Real.hdf5'
        call init_prior_data(prior_path_data(1), path, prior_type, dt_sampling)           
    end subroutine
!==========================================================================================================================   
    
!==========================================================================================================================        
    subroutine read_50_and_71_path(data_directory, prior_type, dt_sampling, measures, prior_path_data)
    !*****************************************************************************************************************
    !Reads prior data under the 50% or 71% path format.
    !
    !:param data_directory: directory where the prior data is stored
    !:type data_directory: str
    !:param prior_type: type of the prior data
    !:type prior_type: str
    !:param dt_sampling: smoothing window duration in second
    !:type dt_sampling: float
    !:param measures: sequence of string saying which measures should be extracted. Allowed
    !          values are 'times', 'B', 'U' and 'ER', order do not matter.
    !:return: MF, U, ER, times, dt_samp, tag
    !:rtype: lists
    !*****************************************************************************************************************
        character(len=*), intent(in) :: data_directory, prior_type
        real(kind=8), intent(in) :: dt_sampling
        character(len=*), intent(in) :: measures(:)
        type(prior_data), intent(inout) :: prior_path_data(:)
        character(len=200) :: path
        integer :: i, j
        logical :: valid
        character(len=20) :: allowed_measures(4) = ['times', 'B', 'U', 'ER']

        ! examine the measure---------------------------------------------------
        valid = .false.
        do i = 1, size(measures)
            if (any(measures(i) == allowed_measures)) valid = .true.
            
            if (.not. valid) then
                print *, "Invalid measure given in measures: ", trim(measures(i))
                print *, "Allowed values are: ", "'times', 'B', 'U', 'ER'"
                stop
            end if
            valid = .false.
        end do
        !-----------------------------------------------------------------------
        path = trim(data_directory) // '/' // 'Real.hdf5'
        call init_prior_data(prior_path_data(1), path, prior_type, dt_sampling)           
    end subroutine
!==========================================================================================================================   
    
!==========================================================================================================================    
    subroutine read_100_path(data_directory, prior_type, dt_sampling, measures, prior_path_data)
    !*****************************************************************************************************************
    !Reads prior data under the 100% path format.
    !
    !:param data_directory: directory where the prior data is stored
    !:type data_directory: str
    !:param prior_type: type of the prior data
    !:type prior_type: str
    !:param dt_sampling: smoothing window duration in second
    !:type dt_sampling: float
    !:param measures: sequence of string saying which measures should be extracted. Allowed
    !          values are 'times', 'B', 'U' and 'ER', order do not matter.
    !:return: MF, U, ER, times, dt_samp, tag
    !:rtype: lists
    !*****************************************************************************************************************
        character(len=*), intent(in) :: data_directory, prior_type
        real(kind=8), intent(in) :: dt_sampling
        character(len=*), intent(in) :: measures(:)
        type(prior_data), intent(inout) :: prior_path_data(:)
        character(len=200) :: path
        integer :: i, j
        logical :: valid
        character(len=20) :: allowed_measures(4) = ['times', 'B', 'U', 'ER']
        character(len=50) :: R(9) = ["Real_A1.hdf5","Real_A2.hdf5", "Real_A3.hdf5", "Real_A4.hdf5", &
         "Real_B1.hdf5", "Real_B2.hdf5", "Real_B3_clean_1.hdf5", "Real_B3_clean_2.hdf5","Real_B4.hdf5"]
        
        ! examine the measure---------------------------------------------------
        valid = .false.
        do i = 1, size(measures)
            if (any(measures(i) == allowed_measures)) valid = .true.
            
            if (.not. valid) then
                print *, "Invalid measure given in measures: ", trim(measures(i))
                print *, "Allowed values are: ", "'times', 'B', 'U', 'ER'"
                stop
            end if
            valid = .false.
        end do
        !-----------------------------------------------------------------------
        do i = 2, 10
            path = trim(data_directory) // '/' // TRIM(R(i-1))
            call init_prior_data(prior_path_data(i), path, prior_type, dt_sampling)
        end do
        path = trim(data_directory) // '/' // '../71path' // '/' //'Real.hdf5'
        call init_prior_data(prior_path_data(1), path, '71path', dt_sampling)
                   
    end subroutine
!==========================================================================================================================     
    
    subroutine read_0_path_test(data_directory, prior_type, dt_sampling, measures, MF, U, ER, times, dt_samp, tag)
    !*****************************************************************************************************************
    !Reads prior data under the 0 path format.
    !
    !:param data_directory: directory where the prior data is stored
    !:type data_directory: str
    !:param prior_type: type of the prior data
    !:type prior_type: str
    !:param dt_sampling: smoothing window duration in second
    !:type dt_sampling: float
    !:param measures: sequence of string saying which measures should be extracted. Allowed
    !          values are 'times', 'B', 'U' and 'ER', order do not matter.
    !:return: MF, U, ER, times, dt_samp, tag
    !:rtype: lists
    !*****************************************************************************************************************
        character(len=*), intent(in) :: data_directory, prior_type
        real(kind=8), intent(in) :: dt_sampling
        character(len=*), intent(in) :: measures(:)
        real(kind=8), intent(inout) :: MF(:,:), U(:,:), ER(:,:), times(:), dt_samp(:)
        character(len=20), intent(inout) :: tag
        logical :: valid
        integer :: i, j, row, column
        character(len=20) :: allowed_measures(4) = ['times', 'B', 'U', 'ER']
        character(len=30) :: m, measure_fn
        character(len=200) :: path
        ! examine the measure---------------------------------------------------
        valid = .false.
        do i = 1, size(measures)
            if (any(measures(i) == allowed_measures)) valid = .true.
            
            if (.not. valid) then
                print *, "Invalid measure given in measures: ", trim(measures(i))
                print *, "Allowed values are: ", "'times', 'B', 'U', 'ER'"
                stop
            end if
            valid = .false.
        end do
        !-----------------------------------------------------------------------

        ! Read the prior data---------------------------------------------------
        do i = 1, SIZE(measures)            
            measure_fn = 'Real' // trim(measures(i)) // '.dat'
            if (measures(i) == 'U') then
                path = trim(data_directory) // '/' // measure_fn
                open(20, file = path)
                call read_line_row(path, row, column)        
                read(20, *) 
                do j = 1, row
                    read(20, *) U(j, :)
                end do                
                close(20)
                U = U * (-1.0d0)
            else if (measures(i) == 'B') then
                path = trim(data_directory) // '/' // measure_fn
                open(20, file = path)
                call read_line_row(path, row, column)   
                read(20, *)
                do j = 1, row
                    read(20, *) MF(j, :)
                end do                
                close(20)
            else if (measures(i) == 'ER') then
                path = trim(data_directory) // '/' // measure_fn
                open(20, file = path)
                call read_line_row(path, row, column)          
                read(20, *)
                do j = 1, row
                    read(20, *) ER(j, :)
                end do                
                close(20)                
            end if            
        end do        
        tag = TRIM(prior_type)
    end subroutine    
!==========================================================================================================================   
    
!==========================================================================================================================        
    subroutine read_50_and_71_path_test(data_directory, prior_type, dt_sampling, measures, MF, U, ER, times, dt_samp, tag)
    !*****************************************************************************************************************
    !Reads prior data under the 50% or 71% path format.
    !
    !:param data_directory: directory where the prior data is stored
    !:type data_directory: str
    !:param prior_type: type of the prior data
    !:type prior_type: str
    !:param dt_sampling: smoothing window duration in second
    !:type dt_sampling: float
    !:param measures: sequence of string saying which measures should be extracted. Allowed
    !          values are 'times', 'B', 'U' and 'ER', order do not matter.
    !:return: MF, U, ER, times, dt_samp, tag
    !:rtype: lists
    !*****************************************************************************************************************
    
        character(len=*), intent(in) :: data_directory, prior_type
        real(kind=8), intent(in) :: dt_sampling
        character(len=*), intent(in) :: measures(:)
        real(kind=8), intent(inout) :: MF(:,:), U(:,:), ER(:,:), times(:), dt_samp(:)
        character(len=20), intent(inout) :: tag
        logical :: valid
        integer :: i, j
        character(len=20) :: allowed_measures(4) = ['times', 'B', 'U', 'ER']
        character(len=30) :: m, measure_fn
        character(len=200) :: path
        INTEGER(HID_T) :: file, space, dset
        INTEGER(HSIZE_T), allocatable :: dims(:)
        integer :: rank, hdferr
        real(kind=8), allocatable :: data_2(:,:), data_1(:)
        ! examine the measure---------------------------------------------------
        valid = .false.
        do i = 1, size(measures)
            if (any(measures(i) == allowed_measures)) valid = .true.
            
            if (.not. valid) then
                print *, "Invalid measure given in measures: ", trim(measures(i))
                print *, "Allowed values are: ", "'times', 'B', 'U', 'ER'"
                stop
            end if
            valid = .false.
        end do
        !-----------------------------------------------------------------------
        
        ! Read the prior data---------------------------------------------------
        path = trim(data_directory) // '/' // 'Real.hdf5'
        ! Initialize FORTRAN interface.
        CALL h5open_f(hdferr)
        
        ! times==========================================================================
        ! Open file
        CALL h5fopen_f(path, H5F_ACC_RDONLY_F, file, hdferr)      
        
        ! Open the dataset
        CALL h5dopen_f(file, '/times', dset, hdferr)
        
        ! Catch the dataspace of the dataset
        CALL h5dget_space_f(dset, space, hdferr)
        
        ! Catch the rank of the dataset
        CALL h5sget_simple_extent_ndims_f(space, rank, hdferr)
        
        ! Catch the dimensions of the dataset
        ALLOCATE(dims(rank))
        CALL h5sget_simple_extent_dims_f(space, dims, dims, hdferr)
        ! Read the data
        ALLOCATE(data_1(dims(1)))
        CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, data_1, dims, hdferr)
        times = data_1
        deallocate(data_1, dims)
        CALL h5dclose_f(dset, hdferr)
        CALL h5fclose_f(file, hdferr)
        !===============================================================================
        
        ! MF============================================================================
        ! Open file
        CALL h5fopen_f(path, H5F_ACC_RDONLY_F, file, hdferr)  
        
        if (trim(prior_type) == '71path') then
            ! Open the dataset
            CALL h5dopen_f(file, '/MF', dset, hdferr)
        
            ! Catch the dataspace of the dataset
            CALL h5dget_space_f(dset, space, hdferr)
        
            ! Catch the rank of the dataset
            CALL h5sget_simple_extent_ndims_f(space, rank, hdferr)
        
            ! Catch the dimensions of the dataset
            ALLOCATE(dims(rank))
            CALL h5sget_simple_extent_dims_f(space, dims, dims, hdferr)
            ! Read the data
            ALLOCATE(data_2(dims(1), dims(2)))
            CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, data_2, dims, hdferr)
            
            MF = TRANSPOSE(data_2) * -1.0d0
            deallocate(data_2, dims)
            CALL h5dclose_f(dset, hdferr)
            CALL h5fclose_f(file, hdferr)
        else
            ! Open the dataset
            CALL h5dopen_f(file, '/MF', dset, hdferr)
        
            ! Catch the dataspace of the dataset
            CALL h5dget_space_f(dset, space, hdferr)
        
            ! Catch the rank of the dataset
            CALL h5sget_simple_extent_ndims_f(space, rank, hdferr)
        
            ! Catch the dimensions of the dataset
            ALLOCATE(dims(rank))
            CALL h5sget_simple_extent_dims_f(space, dims, dims, hdferr)
            ! Read the data
            ALLOCATE(data_2(dims(1), dims(2)))
            CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, data_2, dims, hdferr)
            
            MF = TRANSPOSE(data_2)
            deallocate(data_2, dims)
            CALL h5dclose_f(dset, hdferr)
            CALL h5fclose_f(file, hdferr)
        end if
        !===============================================================================
        
        ! U ============================================================================
        ! Open file
        CALL h5fopen_f(path, H5F_ACC_RDONLY_F, file, hdferr)      
        
        ! Open the dataset
        CALL h5dopen_f(file, '/U', dset, hdferr)
        
        ! Catch the dataspace of the dataset
        CALL h5dget_space_f(dset, space, hdferr)
        
        ! Catch the rank of the dataset
        CALL h5sget_simple_extent_ndims_f(space, rank, hdferr)
        
        ! Catch the dimensions of the dataset
        ALLOCATE(dims(rank))
        CALL h5sget_simple_extent_dims_f(space, dims, dims, hdferr)
        ! Read the data
        ALLOCATE(data_2(dims(1), dims(2)))
        CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, data_2, dims, hdferr)
        U = TRANSPOSE(data_2)
        deallocate(data_2, dims)
        CALL h5dclose_f(dset, hdferr)
        CALL h5fclose_f(file, hdferr)
        !===============================================================================
        
        ! ER ===========================================================================
        ! Open file
        CALL h5fopen_f(path, H5F_ACC_RDONLY_F, file, hdferr)      
        
        ! Open the dataset
        CALL h5dopen_f(file, '/ER', dset, hdferr)
        
        ! Catch the dataspace of the dataset
        CALL h5dget_space_f(dset, space, hdferr)
        
        ! Catch the rank of the dataset
        CALL h5sget_simple_extent_ndims_f(space, rank, hdferr)
        
        ! Catch the dimensions of the dataset
        ALLOCATE(dims(rank))
        CALL h5sget_simple_extent_dims_f(space, dims, dims, hdferr)
        ! Read the data
        ALLOCATE(data_2(dims(1), dims(2)))
        CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, data_2, dims, hdferr)
        ER = TRANSPOSE(data_2)
        deallocate(data_2, dims)
        CALL h5dclose_f(dset, hdferr)
        CALL h5fclose_f(file, hdferr)
        CALL h5close_f(hdferr)
        !===============================================================================
        
        tag = TRIM(prior_type)
        dt_samp = dt_sampling        
    end subroutine
!========================================================================================================================== 
    
!==========================================================================================================================        
    subroutine read_100_path_test(data_directory, prior_type, dt_sampling, measures, MF, U, ER, times, dt_samp, tag)
    !*****************************************************************************************************************
    !Reads prior data under the 100% path format.
    !
    !:param data_directory: directory where the prior data is stored
    !:type data_directory: str
    !:param prior_type: type of the prior data
    !:type prior_type: str
    !:param dt_sampling: smoothing window duration in second
    !:type dt_sampling: float
    !:param measures: sequence of string saying which measures should be extracted. Allowed
    !          values are 'times', 'B', 'U' and 'ER', order do not matter.
    !:return: MF, U, ER, times, dt_samp, tag
    !:rtype: lists
    !*****************************************************************************************************************  
    
        character(len=*), intent(in) :: data_directory, prior_type
        real(kind=8), intent(in) :: dt_sampling
        character(len=*), intent(in) :: measures(:)
        real(kind=8), intent(inout) :: MF(:,:), U(:,:), ER(:,:), times(:), dt_samp(:)
        character(len=20), intent(inout) :: tag
        logical :: valid
        integer :: i, j
        character(len=20) :: allowed_measures(4) = ['times', 'B', 'U', 'ER']
        character(len=30) :: m, measure_fn
        character(len=200) :: path
        INTEGER(HID_T) :: file, space, dset
        INTEGER(HSIZE_T), allocatable :: dims(:, :)
        integer :: rank, hdferr
        real(kind=8), allocatable :: data_2(:,:), data_1(:)
        character(len=50) :: R(9) = ["Real_A1.hdf5","Real_A2.hdf5", "Real_A3.hdf5", "Real_A4.hdf5", &
         "Real_B1.hdf5", "Real_B2.hdf5", "Real_B3_clean_1.hdf5", "Real_B3_clean_2.hdf5","Real_B4.hdf5"]
        ! examine the measure---------------------------------------------------
        valid = .false.
        do i = 1, size(measures)
            if (any(measures(i) == allowed_measures)) valid = .true.
            
            if (.not. valid) then
                print *, "Invalid measure given in measures: ", trim(measures(i))
                print *, "Allowed values are: ", "'times', 'B', 'U', 'ER'"
                stop
            end if
            valid = .false.
        end do
        !-----------------------------------------------------------------------
        
        
        ! loop over time series-------------------------------------------------
        ! times==========================================================================  

        ! ===============================================================================
    end subroutine
!========================================================================================================================== 
    
end module
    