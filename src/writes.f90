module writes
    use utilities   
    use HDF5
    use iso_fortran_env
    use corestate
    implicit none
    
    type :: Hdf5GroupHandler
        character(len=200) :: mother_file
        real(kind=8), allocatable :: dates(:)
    contains
        procedure :: init_Hdf5GroupHandler
        !procedure :: update_all
        !procedure :: close_
    end type
    
contains
!==========================================================================================================================
    subroutine init_Hdf5GroupHandler(self, dates, nb_realisations, core_state, hdf5_file, group_name, format, exclude)
        class(Hdf5GroupHandler), intent(inout) :: self
        real(kind=8), intent(in) :: dates(:)
        integer, intent(in) :: nb_realisations
        class(CoreState_type), intent(in) :: core_state
        character(len=*), intent(in) :: hdf5_file
        character(len=*), intent(in) :: group_name
        character(len=*), intent(in) :: format
        character(len=*), intent(in) :: exclude(:)
        
    
    end subroutine
!==========================================================================================================================
    
!==========================================================================================================================
    subroutine write_hdf5(core_state, dates, hdf5_file, group_name, format_, exclude)
    ! Note:
    ! The precision control parameter `format_` is reserved for specifying
    ! floating-point precision in HDF5 output. For simplicity, this functionality
    ! is not implemented currently. All floating-point data are written using
    ! H5T_NATIVE_DOUBLE (double precision).
        real(kind=8), intent(in) :: dates(:)
        class(CoreState_type), intent(in) :: core_state
        character(len=*), intent(in) :: hdf5_file
        character(len=*), intent(in) :: group_name
        character(len=*), intent(in) :: format_
        character(len=*), intent(in), optional :: exclude(:)
        
        integer(HID_T) :: file_id
        integer(HID_T) :: group_id
        integer :: ierr, i, j
        logical :: is_excluded
        
        
        call h5fopen_f( &
            trim(hdf5_file), &
            H5F_ACC_RDWR_F, &
            file_id, &
            ierr)
        
        call h5gcreate_f( &
            file_id, &
            trim(group_name), &
            group_id, &
            ierr)
        
        ! Write every measure stored in the CoreState. This also supports
        ! CoreStates such as misfits, which contain only MF and SV.
        do i = 1, SIZE(core_state.measures_)
            is_excluded = .false.
            if (PRESENT(exclude)) then
                do j = 1, SIZE(exclude)
                    if (TRIM(core_state.measures_(i).key) == TRIM(exclude(j))) then
                        is_excluded = .true.
                        exit
                    end if
                end do
            end if

            if (.not. is_excluded) then
                call write_matrix( &
                    group_id, &
                    TRIM(core_state.measures_(i).key), &
                    core_state.measures_(i).measure_data &
                )
            end if
        end do
        
        ! write times dataset
        call write_1d(group_id,"times",dates)
        
        call h5gclose_f(group_id,ierr)
        call h5fclose_f(file_id,ierr)
    end subroutine
!==========================================================================================================================    

!==========================================================================================================================    
    subroutine write_matrix(group_id, name, data)
        use hdf5
        implicit none

        integer(HID_T), intent(in) :: group_id
        character(len=*), intent(in) :: name
        real(kind=8), intent(in) :: data(:,:,:)
        real(kind=8), allocatable :: data_tmp(:,:,:)
        integer(HID_T) :: dataset_id
        integer(HID_T) :: dataspace_id
        integer(HSIZE_T) :: dims(3)
        integer :: ierr
        integer :: n1, n2, n3, i, j, k

        n1=size(data,1)
        n2=size(data,2)
        n3=size(data,3)
        allocate(data_tmp(n3,n2,n1))
        dims = shape(data_tmp)
        
        !$omp parallel do collapse(3) default(shared) private(i,j,k)
        do i=1,n1
            do j=1,n2
                do k=1,n3
                    data_tmp(k,j,i) = data(i,j,k)
                end do
            end do
        end do
        !$omp end parallel do
        
        call h5screate_simple_f( &
            3, &
            dims, &
            dataspace_id, &
            ierr)

        call h5dcreate_f( &
            group_id, &
            trim(name), &
            H5T_NATIVE_DOUBLE, &
            dataspace_id, &
            dataset_id, &
            ierr)

        call h5dwrite_f( &
            dataset_id, &
            H5T_NATIVE_DOUBLE, &
            data_tmp, &
            dims, &
            ierr)

        call h5dclose_f(dataset_id,ierr)
        call h5sclose_f(dataspace_id,ierr)
    end subroutine
!==========================================================================================================================
    
!==========================================================================================================================
    subroutine write_1d(group_id, name, data)
        implicit none

        integer(HID_T), intent(in) :: group_id
        character(len=*), intent(in) :: name
        real(kind=8), intent(in) :: data(:)
        integer(HID_T) :: dataspace_id
        integer(HID_T) :: dataset_id
        integer(HSIZE_T) :: dims(1)
        integer :: ierr

        ! dimension
        dims(1)=size(data)

        ! create dataspace
        call h5screate_simple_f( &
            1, &
            dims, &
            dataspace_id, &
            ierr)
        
        ! create dataset
        call h5dcreate_f( &
            group_id, &
            trim(name), &
            H5T_NATIVE_DOUBLE, &
            dataspace_id, &
            dataset_id, &
            ierr)

        ! write data
        call h5dwrite_f( &
            dataset_id, &
            H5T_NATIVE_DOUBLE, &
            data, &
            dims, &
            ierr)

        ! close
        call h5dclose_f(dataset_id,ierr)
        call h5sclose_f(dataspace_id,ierr)
    end subroutine write_1d
!==========================================================================================================================
end module
