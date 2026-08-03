!""" Contains the functions building observations """
module observations
    use reads
    use utilities
    use common
    use constants
    use HDF5
    use config
    implicit none
    
    type :: Observation
    !*****************************************************************************************************************
    !"""
    !Class storing the observation data, operator and errors in members:
    !    * Observation.data
    !    * Observation.H_core
    !    * Observation.errors
    !"""
    !*****************************************************************************************************************
        real(kind=8), allocatable :: X(:,:)
        integer :: nb_obs
        real(kind=8), allocatable :: H(:,:)
        real(kind=8), allocatable :: Rxx(:,:)
    contains
        procedure :: init_Observation
    end type

contains
!==========================================================================================================================    
    subroutine init_Observation(self, X, H, Rxx)
    !*****************************************************************************************************************
    !"""
    !The observation is expected to be used as Y = HX with Y the observation data, X the observed data under spectral form and H the observation operator.
    !
    !
    !:param obs_data: 2D numpy array (nb_real x nb_obs) storing the observation data Y
    !:type obs_data: numpy.ndarray
    !:param H: 2D numpy array (nb_obs x nb_coeffs) storing the observation operator
    !:type H: numpy.ndarray
    !:param R: 2D numpy array (nb_obs x nb_obs) storing the observation errors
    !:type R: numpy.ndarray
    !"""
    !*****************************************************************************************************************
        class(Observation), intent(inout) :: self
        real(kind=8), intent(in) :: X(:,:), H(:,:), Rxx(:,:)
        
        allocate(self.X, source=X)
        self.nb_obs = size(X, 2)
        !# Check that H has first dim nb_obs
        if (SIZE(H, 1) .ne. self.nb_obs) then
            write(10,'(A)') "Observation operator's first dimension must be equal to number of observations."
            write(*,'(A)') "Observation operator's first dimension must be equal to number of observations."
        end if
        allocate(self.H, source=H)
        
        !# Check that R has shape (nb_obs x nb_obs)
        allocate(self.Rxx, source=Rxx)        
    end subroutine    
!==========================================================================================================================
    
!==========================================================================================================================    
    function find_obs_analysis_match(obs_date, t_analysis, dt_f) result(idx)
    !"""
    !return index if obs date found in analysis times +-dt_forecast/2 
    ! -1 means no match found
    !"""
        real(kind=8), intent(in) :: obs_date, dt_f
        real(kind=8), intent(in) :: t_analysis(:)
        integer :: idx
        real(kind=8) :: ana_date
        
        idx = -1
        do idx = 1, size(t_analysis)
            ana_date = t_analysis(idx)
            if (ABS(obs_date - ana_date) <= dt_f/2.0d0) then
                return
            end if
        end do
        idx = -1
    end function   
!==========================================================================================================================
    
!==========================================================================================================================    
    subroutine build_go_vo_observations(cfg, nb_realisations, measure_type, seed)
    !*****************************************************************************************************************
    !"""
    !Builds the observations (including direct observation operators H) from direct sources
    !(Ground observatories (GO), virtual observatories of CHAMP (VO_CHAMP) and SWARM (VO_SWARM).
    !The dates where no analysis takes place are discarded.
    !
    !:param cfg: configuration of the computation
    !:type cfg: inout.config.ComputationConfig
    !:param nb_realisations: Number of realisations of the computation
    !:type nb_realisations: int
    !:param measure_type: Type of the measure to extract ('MF' or 'SV')
    !:type measure_type: str
    !:return: dictionary of Observation objects with dates (np.datetime64) as keys.
    !:rtype: dict[np.datetime64, Observation]
    !"""
    !*****************************************************************************************************************
        class(ComputationConfig), intent(in) :: cfg
        integer, intent(in) :: nb_realisations
        character(len=*), intent(in) :: measure_type
        integer, intent(in) :: seed
        character(len=200) :: datadir
        character(len=20) :: possible_ids(7) = ['GROUND', 'CHAMP', 'SWARM', 'OERSTED', 'CRYOSAT', 'GRACE', 'COMPOSITE'], &
                         cdf_name_shortcuts(7) = ['GO', 'CH', 'SW', 'OR', 'CR', 'GR', 'CO']
        
        datadir = cfg.obs_dir
        
        !# if adding a new satellite, the first two letters of the cdf file should
        !# match the first two letters after '_' in the variable obs_types_to_use.
        !# If not possible, update code lines below that find cdf filename.
        
        !# Adding a correspondence dict
    end subroutine    
!==========================================================================================================================
    
!==========================================================================================================================    
    subroutine build_chaos_hdf5_observations(cfg, nb_realisations, measure_type, seed, observations_)
    !*****************************************************************************************************************
    !"""
    !Builds the observations from a hdf5 file storing CHAOS' spectral coefficients.
    !The dates where no analysis takes place are discarded.
    !
    !:param cfg: configuration of the computation
    !:type cfg: inout.config.ComputationConfig
    !:param nb_realisations: Number of realisations of the computation
    !:type nb_realisations: int
    !:param measure_type: Type of the measure to extract ('MF' or 'SV')
    !:type measure_type: str
    !:return: dictionary of Observation objects with dates (np.datetime64) as keys.
    !:rtype: dict[np.datetime64, Observation]
    !"""
    !*****************************************************************************************************************
        class(ComputationConfig), intent(in) :: cfg
        class(Observation), allocatable, intent(out) :: observations_(:)
        integer, intent(in) :: nb_realisations
        character(len=*), intent(in) :: measure_type
        integer, intent(in) :: seed
        character(len=200) :: datadir, dataset_name
        integer :: max_degree, nb_coefs, nb_coefs_data, max_degree_data, i, j
        logical status
        INTEGER(HID_T) :: file, space, dset
        INTEGER(HSIZE_T), allocatable :: dims(:)
        integer :: rank, hdferr, i_d, match_idx, idx, current_No
        real(kind=8), allocatable :: dates(:), chaos_data(:,:), variance_data(:,:), matrix_temp(:,:)
        real(kind=8), allocatable :: dates_kalmag(:), real_data(:,:,:), matrix_temp2(:,:,:)
        real(kind=8) :: date
        real(kind=8), allocatable :: obs_data(:,:), mean_vec(:), H(:,:), R(:,:)
        
                
        datadir = cfg.obs_dir
        if (TRIM(measure_type) == 'SV') then
            max_degree = cfg.Lsv
            nb_coefs = cfg.Nsv()
            dataset_name = 'dgnm'
        else
            max_degree = cfg.Lb
            nb_coefs = cfg.Nb()
            dataset_name = 'gnm'
        end if
        
        !# Find the hdf5 files (sort to have reproducible order of realisations_files)
        !Ensure the path directory exists; create it if necessary.
        inquire(file=trim(datadir), exist=status) 
        ! only intel fortran have directory option
    
        if (.not. status) then
            write(10,'(A)') "No hdf5 files were found.Check that it is the valid directory path."
            print *, "No hdf5 files were found.Check that it is the valid directory path."
            stop
        end if
        
        write(10,'(A,A,A,A)') "Observations ", TRIM(measure_type), " are read from ", TRIM(datadir)
        write(*,'(A,A,A,A)') "Observations ", TRIM(measure_type), " are read from ", TRIM(datadir)
        
        ! Initialize FORTRAN interface.
        CALL h5open_f(hdferr)
        
        ! Open file
        CALL h5fopen_f(trim(datadir), H5F_ACC_RDONLY_F, file, hdferr)
        
        ! times==========================================================================
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
        ALLOCATE(dates(dims(1)))
        CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, dates, dims, hdferr)
        deallocate(dims)
        CALL h5dclose_f(dset, hdferr)
        !===============================================================================
        
        ! dataset_name==================================================================
        ! Open the dataset
        CALL h5dopen_f(file, '/'//TRIM(dataset_name), dset, hdferr)
        
        ! Catch the dataspace of the dataset
        CALL h5dget_space_f(dset, space, hdferr)
        
        ! Catch the rank of the dataset
        CALL h5sget_simple_extent_ndims_f(space, rank, hdferr)
        
        ! Catch the dimensions of the dataset
        ALLOCATE(dims(rank))
        CALL h5sget_simple_extent_dims_f(space, dims, dims, hdferr)
        
        ! Read the data   
        ALLOCATE(matrix_temp(dims(1), dims(2)), chaos_data(dims(2), dims(1)))
        CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, matrix_temp, dims, hdferr)
        chaos_data = TRANSPOSE(matrix_temp)
        deallocate(dims, matrix_temp)
        CALL h5dclose_f(dset, hdferr)
        !===============================================================================
        
        ! var_dataset_name =============================================================       
        ! Open the dataset
        CALL h5dopen_f(file, '/var_'//TRIM(dataset_name), dset, hdferr)
        
        ! Catch the dataspace of the dataset
        CALL h5dget_space_f(dset, space, hdferr)
        
        ! Catch the rank of the dataset
        CALL h5sget_simple_extent_ndims_f(space, rank, hdferr)
        
        ! Catch the dimensions of the dataset
        ALLOCATE(dims(rank))
        CALL h5sget_simple_extent_dims_f(space, dims, dims, hdferr)
        
        ! Read the data   
        ALLOCATE(matrix_temp(dims(1), dims(2)), variance_data(dims(2), dims(1)))
        CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, matrix_temp, dims, hdferr)
        variance_data = TRANSPOSE(matrix_temp)
        deallocate(dims, matrix_temp)
        CALL h5dclose_f(dset, hdferr)
        !===============================================================================
      
        ! Close file
        CALL h5fclose_f(file, hdferr)
        CALL h5close_f(hdferr)
        
        !---------------------------------------------------------------------------
        nb_coefs_data = SIZE(chaos_data, 2)
        ! # N = L(L+2) => L = sqrt(N+1) - 1
        max_degree_data = INT(SQRT(REAL(nb_coefs_data + 1, kind=8)) - 1)
        
        if (max_degree_data < max_degree) then
            write(10,'(A,I3,A,A,A,I3,A)') "There are not enough coefficients in the CHAOS file to handle the max degree", max_degree, " asked for ", trim(measure_type), ".  Please retry with a lower degree (max: ", max_degree_data ,")."
            write(*,'(A,I3,A,A,A,I3,A)') "There are not enough coefficients in the CHAOS file to handle the max degree", max_degree, " asked for ", trim(measure_type), ".  Please retry with a lower degree (max: ", max_degree_data ,")."
        endif
        
        ! # find kalmag file        
        datadir = './data/observations/KALMAG/KALMAG.hdf5'
        !Ensure the path directory exists; create it if necessary.
        inquire(file=trim(datadir), exist=status) 
        ! only intel fortran have directory option
    
        if (.not. status) then
            write(10,'(A)') "KALMAG not found! Check the path!"
            print *, "KALMAG not found! Check the path!"
            stop
        end if
        
        ! Initialize FORTRAN interface.
        CALL h5open_f(hdferr)
        
        ! Open file
        CALL h5fopen_f(trim(datadir), H5F_ACC_RDONLY_F, file, hdferr)
        
        ! times==========================================================================
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
        ALLOCATE(dates_kalmag(dims(1)))
        CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, dates_kalmag, dims, hdferr)
        deallocate(dims)
        CALL h5dclose_f(dset, hdferr)
        !===============================================================================
        
        ! measure_type =================================================================       
        ! Open the dataset
        CALL h5dopen_f(file, '/'//TRIM(measure_type), dset, hdferr)
        
        ! Catch the dataspace of the dataset
        CALL h5dget_space_f(dset, space, hdferr)
        
        ! Catch the rank of the dataset
        CALL h5sget_simple_extent_ndims_f(space, rank, hdferr)
        
        ! Catch the dimensions of the dataset
        ALLOCATE(dims(rank))
        CALL h5sget_simple_extent_dims_f(space, dims, dims, hdferr)
        
        ! Read the data   
        ALLOCATE(matrix_temp2(dims(1), dims(2), dims(3)), real_data(dims(3), dims(2), dims(1)))
        CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, matrix_temp2, dims, hdferr)
        do i = 1, dims(1)
            real_data(:,:,i) = TRANSPOSE(matrix_temp2(i,:,:))
        end do
        deallocate(dims, matrix_temp2)
        CALL h5dclose_f(dset, hdferr)
        !===============================================================================
      
        ! Close file
        CALL h5fclose_f(file, hdferr)  
        CALL h5close_f(hdferr)
        
        ! # Format the data as a dict of Observations with dates as keys
        allocate(obs_data(nb_realisations, nb_coefs), source=0.0d0)
        allocate(mean_vec(nb_coefs), source=0.0d0)
        allocate(H(nb_coefs, nb_coefs), source=0.0d0)
        allocate(R(nb_coefs, nb_coefs), source=0.0d0)
        allocate(observations_(size(cfg.t_analyses_full)))
        do i_d = 1, size(dates)
            date = dates(i_d)
            !# Shift the date by one month (different convention)
            match_idx = find_obs_analysis_match(date, cfg.t_analyses_full, cfg.dt_f)
            if (match_idx == -1) then
                cycle
            end if
            !# Get asked nb_realisations and nb_coefs
            do j = 1, nb_realisations
                obs_data(j,:) = chaos_data(i_d, 1:nb_coefs)
            end do
            
            !#Now we build full obs_data = obs_data + N(0,sigma_mods_error) + deviation_from_mean_reals_kalmag
            idx = FINDLOC(dates_kalmag, date, dim=1)
            !# Add covobs deviation from mean to obs_data
            mean_vec = 0.0d0
            do i = 1, size(real_data,1)
                mean_vec = mean_vec + real_data(i, idx, 1:nb_coefs)
            end do
            mean_vec = mean_vec / real(size(real_data,1), kind=8)
            
            do i = 1, nb_realisations
                obs_data(i, :) = obs_data(i, :) + real_data(i, idx, 1:nb_coefs) - mean_vec
            end do
            
            !# Get the number of observed coefficients (here equal to nb_coefs)
            current_No = SIZE(obs_data, 2)
            !# H is identity (spectral data)
            call dlaset('A', nb_coefs, nb_coefs, 0.0d0, 1.0d0, H, nb_coefs)
            !# R is a diagonal matrix with variances as diagonal
            do i = 1, nb_coefs
                R(i, i) = variance_data(i_d, i)
            end do
            call init_Observation(observations_(match_idx), obs_data, H, R)
        end do
    end subroutine    
!==========================================================================================================================
end module
    