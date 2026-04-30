module utilities
!==================================
!Utility functions
!In python, index begin at 0
!In fortran, index begin at 1
!==================================
    use HDF5
    use blas95
    use lapack95
    use f95_precision
    implicit none
    public
    
    real(kind=8), parameter :: pi = 4.0d0 * atan(1.0d0)
    
    type :: obs_type
        integer :: seed
        character(len=20) :: num_obs
        character(len=20) :: measure_name
    end type
    
    type :: prior_data
        integer :: dim_times, dim_MF_1, dim_MF_2, dim_U_1, dim_U_2, dim_ER_1, dim_ER_2
        real(kind=8), allocatable :: times(:), MF(:, :), U(:, :), ER(:, :)
        character(len=20) :: tag
        real(kind=8) :: dt_samp
    contains
        procedure :: init_prior_data
    end type
       
    type :: config_data
        integer :: Lu = -1, Lb = -1, Lsv = -1, Ly = -1, dt_a_f_ratio = -1, N_dta_start_forecast = -1, compute_e = -1
        integer :: N_dta_end_forecast = -1, Nth_legendre = -1, N_pca_u = -1
        integer :: out_analysis = -1, out_forecast = -1, out_computed = -1, out_misfits = -1, combined_U_ER_forecast = -1
        integer :: export_every_analysis = -1, last_analysis_backward = -1
        real(kind=8) :: TauU = -1.0d0, TauE = -1.0d0, TauG = -1.0d0
        real(kind=8) :: t_start_analysis = -1.0d0, t_end_analysis = -1.0d0, dt_f = -1.0d0, dt_smoothing = -1.0d0, dt_sampling = -1.0d0
        real(kind=8) :: init_date = -1.0d0, remove_spurious = -1.0d0, remove_spurious_shear_u = -1.0d0, remove_spurious_shear_err = -1.0d0, discard_high_lat = -1.0d0
        character(len=200) :: prior_dir = 'none', prior_dir_shear = 'none', prior_type = 'none', prior_type_shear = 'none', init_file = 'none', SW_err = 'none', out_format = 'none'
        character(len=200) :: obs_dir = 'none', obs_type = 'none', AR_type = 'none', pca_norm = 'none', core_state_init = 'none', kalman_norm = 'none'
    end type
    
    type :: cov_prior_type
        real(kind=8), allocatable :: Z_Z(:,:), B_B(:,:), U_U(:,:), ER_ER(:,:)
        real(kind=8), allocatable :: dZ_dZ(:,:), d2Z_d2Z(:,:)
        real(kind=8), allocatable :: A(:,:), B(:,:), C(:,:), Chol(:,:)
    end type
    
    type :: set_prior_type ! augkf_algo.f90 line:280 all U, MF, times, ER are stored in this type
        real(kind=8), allocatable :: U(:,:), ER(:,:), MF(:,:)
        real(kind=8), allocatable :: times(:)
    end type
    
    type :: container_type ! augkf_algo.f90 line:470 all Z.T, dZ.T * dt_prior, dt_prior, Nt are stored in this type
        real(kind=8), allocatable :: z_T(:,:), dz_T(:,:)
        real(kind=8) :: dt_prior
        integer :: Nt
    end type
    
    type :: WWT ! common.f90 in subroutine compute_AR_coefs_avg  
        real(kind=8), allocatable :: matrix(:,:)
    end type
    
    type :: legendre_polys_type
        real(kind=8), allocatable :: thetas(:), weights(:)
        real(kind=8), allocatable :: MF(:,:,:), U(:,:,:), SV(:,:,:)
    end type
    
    type :: input_core_state_type
        !real(kind=8), allocatable :: MF(:,:), U(:,:), ER(:,:)
        real(kind=8), allocatable :: B(:)
        integer :: Nu2, Nsv, Nb
        integer :: Lsv, Lu, Lb
    end type
    
contains

!========================================================================================================================== 
    function sum_blackman(Nt) result(blackman_w)
        integer :: Nt,i
        real(kind=8) :: blackman_w
        real(kind=8) :: w(Nt)
        
        do i = 1, Nt
            w(i) = 0.42d0 - 0.5d0 * cos(2.0d0*pi*(i-1)/(Nt-1)) &
                         + 0.08d0 * cos(4.0d0*pi*(i-1)/(Nt-1))
        end do
        blackman_w = SUM(w)
    end function
!==========================================================================================================================
    
!==========================================================================================================================
    subroutine init_container(container, zT, dzT, dt_prior, Nt)
        class(container_type), intent(out) :: container
        real(kind=8), intent(in) :: zT(:,:), dzT(:,:)
        real(kind=8), intent(in) :: dt_prior
        integer, intent(in) :: Nt
        
        allocate(container.z_T, source=zT)
        allocate(container.dz_T, source=dzT)
        container.dt_prior = dt_prior
        container.Nt = Nt
    end subroutine
!========================================================================================================================== 
    
!========================================================================================================================== 
    subroutine np_concatenate0(matrix1, matrix2, matrix3)
    !***************************************************************************************************************"""
    !np_concatenate
    !matrix3 = [matrix1
    !           matrix2]
    !*****************************************************************************************************************"""
        real(kind=8), intent(in) :: matrix1(:,:), matrix2(:,:)
        real(kind=8), allocatable, intent(out) :: matrix3(:,:)
        integer :: i
        
        allocate(matrix3(SIZE(matrix1, 1)+SIZE(matrix2, 1), SIZE(matrix1, 2)))
        matrix3(1:SIZE(matrix1, 1), :) = matrix1
        matrix3(SIZE(matrix1, 1)+1:SIZE(matrix1, 1)+SIZE(matrix2, 1), :) = matrix2
    end subroutine
!==========================================================================================================================
    
!========================================================================================================================== 
    subroutine block_diag(matrix1, matrix2, matrix3)
    !***************************************************************************************************************"""
    !sc.linalg.block_diag
    !matrix3 = [matrix1 0
    !           0       matrix2]
    !*****************************************************************************************************************"""
        real(kind=8), intent(in) :: matrix1(:,:), matrix2(:,:)
        real(kind=8), allocatable, intent(out) :: matrix3(:,:)
        integer :: i
        
        allocate(matrix3(SIZE(matrix1, 1)+SIZE(matrix2, 1), SIZE(matrix1, 2)+SIZE(matrix2, 2)))
        matrix3 = 0.0d0
        matrix3(1:SIZE(matrix1, 1), 1:SIZE(matrix1, 2)) = matrix1
        matrix3(SIZE(matrix1, 1)+1:SIZE(matrix1, 1)+SIZE(matrix2, 1), SIZE(matrix1, 2)+1:SIZE(matrix1, 2)+SIZE(matrix2, 2)) = matrix2
    end subroutine
!==========================================================================================================================

!==========================================================================================================================    
    subroutine max_inv(x, inv_x)
    ! matrix_inverse
        real(kind=8), intent(in) :: x(:, :)
        real(kind=8), allocatable, intent(out) :: inv_x(:, :)
        integer, allocatable :: ipiv(:)
        REAL(kind=8), allocatable :: matrix(:, :)
        
        allocate(ipiv(SIZE(x, 1)))
        allocate(inv_x, source=x)
        allocate(matrix, source=x)
        call getrf(matrix,ipiv)
        call getri(matrix,ipiv)
        inv_x = matrix    
    end subroutine
!==========================================================================================================================
    
!========================================================================================================================== 
    subroutine blackman(Nt, result_) 
        integer, intent(in) :: Nt
        integer :: i
        real(kind=8), allocatable, intent(out) :: result_(:)
        
        allocate(result_(Nt))
        if (Nt == 1) then
            result_(1) = 1.0d0
        else
            do concurrent (i = 1: Nt)
                result_(i) = 0.42d0 - 0.5d0 * cos(2.0d0*pi*(i-1)/(Nt-1)) &
                             + 0.08d0 * cos(4.0d0*pi*(i-1)/(Nt-1))
            end do            
        end if
    end subroutine
!==========================================================================================================================

!========================================================================================================================== 
    subroutine init_prior_data(this, path, prior_type, dt_sampling)
        class(prior_data), intent(inout) :: this
        character(len=*), intent(in) :: path, prior_type
        real(kind=8), intent(in) :: dt_sampling
        INTEGER(HID_T) :: file, space, dset
        INTEGER(HSIZE_T), allocatable :: dims(:)
        integer :: rank, hdferr
        real(kind=8), allocatable :: data_2(:,:), data_1(:)
        
        ! Initialize FORTRAN interface.
        CALL h5open_f(hdferr)
        
        ! times==========================================================================
        ! Open file
        CALL h5fopen_f(trim(path), H5F_ACC_RDONLY_F, file, hdferr)
        
        !! Start SWMR mode (mpi mode)
        !call h5fstart_swmr_read_f(file, hdferr)
        
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
        this.dim_times = dims(1)
        
        if (trim(prior_type) == '0path') then
            ALLOCATE(data_1(dims(1)), this.times(dims(1)))
            CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, data_1, dims, hdferr)
            this.times = data_1
        else if (trim(prior_type) == '50path') then
            ALLOCATE(data_1(dims(1)), this.times(dims(1)))
            CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, data_1, dims, hdferr)
            this.times = data_1
        else if (trim(prior_type) == '71path') then
            ALLOCATE(data_1(dims(1)), this.times(dims(1)))
            CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, data_1, dims, hdferr)
            this.times = data_1
        else if (trim(prior_type) == '100path') then
            ALLOCATE(data_1(dims(1)), this.times(dims(1) - 50))
            CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, data_1, dims, hdferr)
            this.times = data_1(51:)
        end if

        deallocate(data_1, dims)
        CALL h5dclose_f(dset, hdferr)
        CALL h5fclose_f(file, hdferr)
        !===============================================================================
        
        ! MF============================================================================
        ! Open file
        CALL h5fopen_f(trim(path), H5F_ACC_RDONLY_F, file, hdferr)  
        
        !! Start SWMR mode (mpi mode)
        !call h5fstart_swmr_read_f(file, hdferr)
        
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
        this.dim_MF_1 = dims(1)
        this.dim_MF_2 = dims(2)
        
        if (trim(prior_type) == '0path') then
            ALLOCATE(data_2(dims(1), dims(2)), this.MF(dims(2), dims(1)))
            CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, data_2, dims, hdferr)
            this.MF = TRANSPOSE(data_2)
        else if (trim(prior_type) == '50path') then
            ALLOCATE(data_2(dims(1), dims(2)), this.MF(dims(2), dims(1)))
            CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, data_2, dims, hdferr)
            this.MF = TRANSPOSE(data_2)
        else if (trim(prior_type) == '71path') then
            ALLOCATE(data_2(dims(1), dims(2)), this.MF(dims(2), dims(1)))
            CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, data_2, dims, hdferr)
            this.MF = TRANSPOSE(data_2) * -1.0d0
        else if (trim(prior_type) == '100path') then
            ALLOCATE(data_2(dims(1), dims(2)), this.MF(dims(2) - 50, dims(1)))
            CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, data_2, dims, hdferr)
            this.MF = TRANSPOSE(data_2(:,51:))
        end if
        
        deallocate(data_2, dims)
        CALL h5dclose_f(dset, hdferr)
        CALL h5fclose_f(file, hdferr)
        !===============================================================================
        
        ! U ============================================================================
        ! Open file
        CALL h5fopen_f(trim(path), H5F_ACC_RDONLY_F, file, hdferr)
        
        !! Start SWMR mode (mpi mode)
        !call h5fstart_swmr_read_f(file, hdferr)
        
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
        this.dim_U_1 = dims(1)
        this.dim_U_2 = dims(2)
        
        if (trim(prior_type) == '0path') then
            ALLOCATE(data_2(dims(1), dims(2)), this.U(dims(2), dims(1)))
            CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, data_2, dims, hdferr)        
            this.U = TRANSPOSE(data_2)
        else if (trim(prior_type) == '50path') then
            ALLOCATE(data_2(dims(1), dims(2)), this.U(dims(2), dims(1)))
            CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, data_2, dims, hdferr)        
            this.U = TRANSPOSE(data_2)
        else if (trim(prior_type) == '71path') then
            ALLOCATE(data_2(dims(1), dims(2)), this.U(dims(2), dims(1)))
            CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, data_2, dims, hdferr)        
            this.U = TRANSPOSE(data_2)
        else if (trim(prior_type) == '100path') then
            ALLOCATE(data_2(dims(1), dims(2)), this.U(dims(2) - 50, dims(1)))
            CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, data_2, dims, hdferr)        
            this.U = TRANSPOSE(data_2(:,51:))
        end if        
        
        deallocate(data_2, dims)
        CALL h5dclose_f(dset, hdferr)
        CALL h5fclose_f(file, hdferr)
        !===============================================================================
        
        ! ER ===========================================================================
        ! Open file
        CALL h5fopen_f(trim(path), H5F_ACC_RDONLY_F, file, hdferr)
        
        !! Start SWMR mode (mpi mode)
        !call h5fstart_swmr_read_f(file, hdferr)
        
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
        this.dim_ER_1 = dims(1)
        this.dim_ER_2 = dims(2)
        
        if (trim(prior_type) == '0path') then
            ALLOCATE(data_2(dims(1), dims(2)), this.ER(dims(2), dims(1)))
            CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, data_2, dims, hdferr)        
            this.ER = TRANSPOSE(data_2)
        else if (trim(prior_type) == '50path') then
            ALLOCATE(data_2(dims(1), dims(2)), this.ER(dims(2), dims(1)))
            CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, data_2, dims, hdferr)        
            this.ER = TRANSPOSE(data_2)
        else if (trim(prior_type) == '71path') then
            ALLOCATE(data_2(dims(1), dims(2)), this.ER(dims(2), dims(1)))
            CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, data_2, dims, hdferr)        
            this.ER = TRANSPOSE(data_2)
        else if (trim(prior_type) == '100path') then
            ALLOCATE(data_2(dims(1), dims(2)), this.ER(dims(2) - 50, dims(1)))
            CALL h5dread_f(dset, H5T_NATIVE_DOUBLE, data_2, dims, hdferr)        
            this.ER = TRANSPOSE(data_2(:,51:))
        end if           
        
        deallocate(data_2, dims)
        CALL h5dclose_f(dset, hdferr)
        CALL h5fclose_f(file, hdferr)
        CALL h5close_f(hdferr)
        !===============================================================================     
        
        if (trim(prior_type) == '0path') then
            ! tag===========================================================================
            this.tag = TRIM(prior_type)
            this.dt_samp = 0
        else if (trim(prior_type) == '50path' .OR. (trim(prior_type) == '71path')) then
            this.tag = TRIM(prior_type)
            this.dt_samp = dt_sampling
        else if (trim(prior_type) == '100path') then
            this.tag = TRIM(prior_type)
            this.dt_samp = this.times(2) - this.times(1)  
        end if
    end subroutine
!========================================================================================================================== 
    
!==========================================================================================================================  
    subroutine decimal_to_date(decimal_date, month_shift, year, month)
    !***************************************************************************************************************"""
    !Converts a decimal date in a datetime64 object.
    !Assumes decimal_date = YEAR + (MONTH + SHIFT)/12.
    !
    !:param decimal_date: date to convert
    !:type decimal_date: double
    !:param month_shift: shift of the month. Integer dates have Dec as month with 0, whereas they have Jan with 1.
    !:type month_shift: 0 or 1
    !:return: date converted
    !:rtype: int year month
    !*****************************************************************************************************************"""
        real(kind=8), intent(in) :: decimal_date
        integer, intent(in), optional :: month_shift
        integer, intent(out) :: year, month        
        integer :: month_shift_value
        ! 设置默认值
        if (present(month_shift)) then
            month_shift_value = month_shift
        else
            month_shift_value = 0  
        end if
    
        year = INT(decimal_date)
        month = NINT((decimal_date - year) * 12.0) + month_shift_value
    
        ! For Dec (year + 12/12), the previous formula returns 0 (if shift is 0) so need to treat the case apart
        if (month == 0) then
            month = 12
            year = year - 1
        end if    
    end subroutine
!==========================================================================================================================    
    

!==========================================================================================================================
    subroutine date_to_decimal(year, month, nb_decimals, decimal_date)
    !***************************************************************************************************************"""
    !Converts a datetime64 object in a decimal date.
    !Returns YEAR + MONTH/12. rounded to the number of decimals given
    !
    !:param date: date to convert
    !:type date: int  year, month 
    !:param nb_decimals: number of decimals for round-off (default: None)
    !:type nb_decimals: int or None(-1)
    !:return: date in decimal format
    !:rtype: double decimal_date
    !*****************************************************************************************************************"""
        integer, intent(in) :: year, month
        integer, intent(in), optional :: nb_decimals
        real(kind=8), intent(out) :: decimal_date
        integer :: nb_decimals_value
        
        ! 设置默认值
        if (present(nb_decimals)) then
            nb_decimals_value = nb_decimals
        else
            nb_decimals_value = -1 
        end if
        
        if (nb_decimals_value == -1) then
            decimal_date = year + month/12.0d0
        else
            decimal_date = (NINT((year + month/12.0d0) * (10.0d0 ** nb_decimals_value))) / (10.0d0 ** nb_decimals_value)
        end if
    end subroutine
!==========================================================================================================================

!==========================================================================================================================
    subroutine date_array_to_decimal(year_array, month_array, nb_decimals, decimal_date_array)
    !***************************************************************************************************************"""
    !Uses date_to_decimal to convert a 1D array of datetime64 in a 1D array of floating numbers.

    !:param year_array, month_array: array to convert
    !:type year_array, month_array: 1D array of int
    !:param nb_decimals: Number of decimals for round-off (default: None)
    !:type nb_decimals: int or None(-1)
    !:return: Same array with dates as floating numbers
    !:rtype: 1D array double
    !*****************************************************************************************************************"""
        implicit none
        integer, intent(in) :: year_array(:), month_array(:)
        integer, intent(in), optional :: nb_decimals
        real(kind=8), intent(out) :: decimal_date_array(:)
        
        integer :: i, nb_decimals_value
        
        if (present(nb_decimals)) then
            nb_decimals_value = nb_decimals
        else
            nb_decimals_value = -1 
        end if
        
        do i = 1, size(year_array)
            call date_to_decimal(year_array(i), month_array(i), nb_decimals_value, decimal_date_array(i))
        end do
                
    end subroutine
!==========================================================================================================================
    
!==========================================================================================================================
    subroutine eval_Plm(z, thetas, legendre_poly_values_at_thetas, parity, result_out)
    !***************************************************************************************************************"""
    !Evaluate the Legendre polynomial at l,m at an arbitrary z using interpolation.

    !:param z: Value where the Legendre polynomial must be evaluated
    !:type z: double
    !:param thetas: List of the angles used for the computation of Legendre polynomials
    !:type thetas: array double
    !:param legendre_poly_values_at_thetas: List of the Legendre polynomials value for all cos(thetas) (at fixed l,m) (same dimension as thetas)
    !:type legendre_poly_values_at_thetas: array double
    !:param parity: Parity of the Legendre polynomial (0 for even, 1 for odd)
    !:type legendre_poly_values_at_thetas: array double
    !:return: Interpolated value of the Legendre polynomial at z
    !:rtype: result_out double
    !待定
    !*****************************************************************************************************************"""
        real(kind=8), intent(in) :: z
        real(kind=8), intent(in) :: thetas(:), legendre_poly_values_at_thetas(:)
        integer, intent(in) :: parity
        real(kind=8), intent(out) :: result_out
        
        if (z == 1.0d0) then
            result_out = 1.0d0
        else if (z == -1.0d0) then
            result_out = -1.0d0 ** parity
        else
            result_out = 0.0d0
        end if
    end subroutine
    
!==========================================================================================================================

!==========================================================================================================================    
    subroutine get_degree_order(k, degree, order, name)
    !*****************************************************************************************************************
    !Gets the degree, order and coefficient name ("g" or "h") of a coefficient.
    
    !(the same as python , index begin at 0)  eg: k = 5, return 2, 1, "h"
    !k = 5 => 2, 1, "h"   means:
    !coef_index(6)  is h_21 in fortran index 
    
    !:param k: index of the coefficient
    !:type k: int
    !:return: degree, order and coef name of the coefficient
    !:rtype: int, int, str
    !*****************************************************************************************************************
        integer, intent(in) :: k
        integer, intent(out) :: degree, order
        character(len=1), intent(out) :: name        
        real(kind=8) :: sqrt_k
        integer :: twice_order
        
        if (k .lt. 0) then
            print *, "k must be .ge. 0"
            call exit()
        end if
        
        sqrt_k = (k+1) ** 0.5d0
        degree = INT(sqrt_k)
        
        if (degree == sqrt_k) then
            order = 0
            name = "g"
        else 
            twice_order = k - degree ** 2 + 2
            if (mod(twice_order, 2) == 0) then
                name = "g"
            else
                name = "h"
                twice_order = twice_order - 1
            end if
            order = twice_order / 2                        
        end if
        !# We need now to find m verifying:
        !#    for g : k = l**2 + 2m - 2 => 2m = k - l**2 + 2
        !# OR for h : k = l**2 + 2m - 1 => 2m = k - l**2 + 1

    end subroutine
!==========================================================================================================================

!==========================================================================================================================
    subroutine north_geomagnetic_pole_angle(g_10, g_11, h_11, theta_0, phi_0)
    !*****************************************************************************************************************
    !Given the first three elements of the spherical harmonic decomposition,
    !return the angle theta and phi of the geomagnetic pole.
    !Details in:
    !Hulot, G., et al. "The present and future geomagnetic field." (2015): 33-78.
    !input: double g_10, g_11, h_11
    !output: double theta0, phi0 (rad)
    !*****************************************************************************************************************
        real(kind=8), intent(in) :: g_10, g_11, h_11
        real(kind=8), intent(out) :: theta_0, phi_0
        real(kind=8) :: m_0
        
        m_0 = SQRT(g_10 ** 2.0d0 + g_11 ** 2.0d0 + h_11 ** 2.0d0)
        theta_0 = pi - ACOS(g_10 / m_0)
        phi_0   = -pi + ATAN2(h_11, g_11)
    
    end subroutine
!==========================================================================================================================

!==========================================================================================================================
    subroutine geomagnetic_latitude(theta, phi, theta_0, phi_0, geo_latitude)
    !*****************************************************************************************************************
    !Computes the geomagnetic (or dipole) latitude.
    !Angle must be given in radians.
    !theta_0: latitude of the north geomagnetic pole
    !phi_0: longitude of the north geomagnetic pole
    !type: double theta, phi, theta_0, phi_0, geo_latitude
    !*****************************************************************************************************************
        real(kind=8), intent(in) :: theta, phi, theta_0, phi_0
        real(kind=8), intent(out) :: geo_latitude
        real(kind=8) :: cos_term, sin_term
        
        cos_term = COS(theta) * COS(theta_0)
        sin_term = SIN(theta) * SIN(theta_0) * COS(phi - phi_0)
        
        geo_latitude = ACOS(cos_term + sin_term)        
    
    end subroutine
!==========================================================================================================================

!==========================================================================================================================
    subroutine extract_dipole_chaos(chaos_file, times, coef)
    !*****************************************************************************************************************
    !Extract the coefficients of the dipole (g_10, g_11, h_11) from chaos_file,
    !at all available times.
    !*****************************************************************************************************************   
        character(len=*), intent(in) :: chaos_file
        real(kind=8), intent(inout) :: times(:), coef(:,:)
        integer :: epoches, i, j, n, m
        
        open(10, file=chaos_file, status='old')
            do i = 1, 3
                read(10, *)
            end do
            read(10, *) n, m, epoches
            read(10, *) (times(i), i = 1, epoches)
            do i = 1, 3
                read(10, *) n, m, (coef(j,i), j = 1, epoches)
            end do
        close(10)
    end subroutine
!==========================================================================================================================

!==========================================================================================================================
    subroutine zonal_indexes(max_degree, n_array)
    !*****************************************************************************************************************
    !Generates the indexes of zonal coefficients g_n^0.
    
    !(the same as python , index begin at 0)  eg: max_degree = 3, return [0, 3, 8] in python_index
    !degree = 0, 1, 2  =>   INDEX = 0, 3, 8 means:
    !g_10 => coef_index(1), g_20 => coef_index(4), g_30 => coef_index(9) in fortran array_index 
    
    !:param max_degree: Maximum value of n 
    !:type max_degree: int
    !:return: List of indexes of zonal coefficients
    !:rtype: array int
    !*****************************************************************************************************************
        integer, intent(in) :: max_degree
        integer, intent(out) :: n_array(:)
        integer :: i
        do i = 0, max_degree-1
            n_array(i+1) = i * (i + 2)
        end do
    end subroutine
!==========================================================================================================================

!==========================================================================================================================
    subroutine zonal_mask(max_degree, n_log_array)
    !*****************************************************************************************************************
    !Generates a mask to select indexes of zonal coefficients g_n^(m=0).
    
    !(the same as python , index begin at 0)  eg: max_degree = 3, return
    ![ True False False  True False False False False  True False False False
    !  False False False False False False False False False False False False
    !  False False False False False]
    
    !:param max_degree: Maximum value of n
    !:type max_degree: int
    !:return: List of logical of zonal coefficients  
    !!!! len = 2 * max_degree * (max_degree + 2) attention!!!!
    !:rtype: array logical
    !*****************************************************************************************************************
        integer, intent(in) :: max_degree
        logical, intent(out) :: n_log_array(:)
        integer :: index_(max_degree)
        
        call zonal_indexes(max_degree, index_)
        n_log_array = .false.
        index_ = index_ + 1
        n_log_array(index_) = .true.
    end subroutine
!==========================================================================================================================   

!==========================================================================================================================
    subroutine non_zonal_mask(max_degree, n_non_log_array)
    !*****************************************************************************************************************
    !Generates a mask to select indexes of zonal coefficients g_n^(m!=0).
    
    !(the same as python , index begin at 0)  eg: max_degree = 3, return
    ![False  True  True False  True  True  True  True False  True  True  True
    ! True  True  True  True  True  True  True  True  True  True  True  True
    ! True  True  True  True  True  True]
    
    !:param max_degree: Maximum value of n
    !:type max_degree: int
    !:return: List of logical of non zonal coefficients  
    !!!! len = 2 * max_degree * (max_degree + 2) attention!!!!
    !:rtype: array logical
    !*****************************************************************************************************************
        integer, intent(in) :: max_degree
        logical, intent(out) :: n_non_log_array(:)
        integer :: index_(max_degree)
        
        call zonal_indexes(max_degree, index_)
        n_non_log_array = .true.
        index_ = index_ + 1
        n_non_log_array(index_) = .false.
    end subroutine
!==========================================================================================================================   
    
!==========================================================================================================================   
    subroutine spectral_to_znz_matrix(max_degree, P)
    !*****************************************************************************************************************
    !Generates the permutation matrix P linking U stored in spectral order to U stored in zonal-non_zonal order
    !
    !znz_U = P @ spectral_U
    !
    !:param max_degree: Maximum value of n
    !:type max_degree: int
    !:return: Matrix linking spectral U to zonal-non_zonal U
    !:rtype: np.array (dim: Nu2 x Nu2)
    !!!!   Nuu2 = 2 * max_degree * (max_degree + 2) attention!!!!
    !*****************************************************************************************************************    
        integer, intent(in) :: max_degree
        integer, allocatable, intent(out) :: P(:,:)
        integer  Nu2, i, k
        logical :: n_log_array(2 * max_degree * (max_degree + 2)), n_non_log_array(2 * max_degree * (max_degree + 2))        
        integer :: n_array(2 * max_degree * (max_degree + 2))
        integer, allocatable :: index_(:)
        
        Nu2 = 2 * max_degree * (max_degree + 2)
        allocate(P(Nu2, Nu2))
        call zonal_mask(max_degree, n_log_array)
        call non_zonal_mask(max_degree, n_non_log_array)
        P = 0
        
        do i = 1, Nu2
            n_array(i) = i
        end do        
        
        ! # Link the first members to zonal coefs
        index_ = PACK(n_array, n_log_array)
        do i = 1, max_degree
            P(i, index_(i)) = 1
        end do
        !Link the rest to non zonal coefs
        index_ = PACK(n_array, n_non_log_array)
        do i = max_degree+1, Nu2
            P(i, index_(i-max_degree)) = 1
        end do
    end subroutine
!==========================================================================================================================   
    
!========================================================================================================================== 
    subroutine get_seeds_for_obs(seed, obs_type_to_use)
    !*****************************************************************************************************************
    !Get a seed for each observatory to noise the MF or SV 
    
    !:param seed: seed value of n, obs_type_to_use
    !:type seed: int, type(obs_type)
    !:return: obs_type_to_use
    !:rtype: type(obs_type)
    !*****************************************************************************************************************      
        integer, intent(in) :: seed
        type(obs_type), intent(inout) :: obs_type_to_use(:) 
        integer :: seed_(1)
        integer :: i
        real(kind=8) :: seed_float
        
        seed_(1) = seed
        call random_seed(put = seed_)
        
        do i = 1, size(obs_type_to_use)
            call random_number(seed_float)
            obs_type_to_use(i).seed = INT(seed_float * 50000)
        end do      
    
    end subroutine
!==========================================================================================================================   
    
!==========================================================================================================================     
    
    subroutine read_line_row(path, row, column)
    !*****************************************************************************************************************
    !Reads the columns and rows of the file.
    !
    !:param path: directory where the file is stored
    !:type data_directory: str
    !:param row: row number
    !:type row: float
    !:param line: column number
    !:type line: float
    !:rtype: row, column
    !*****************************************************************************************************************
        character(len=*), intent(in):: path
        integer, intent(inout) :: row, column        
        open(unit=30, file=path)        
        
        read(30, *) row, column
        
        close(30)
    
    end subroutine
    
end module
!==========================================================================================================================   
    
!========================================================================================================================== 
    
    
    
    
    
    
    
    
    
!!!!!!!草稿    

!!==========================================================================================================================
!    subroutine non_zonal_mask(max_degree, n_non_array)
!    !*****************************************************************************************************************
!    !Generates a mask to select indexes of zonal coefficients g_n^(m!=0).
!    
!    !(the same as python , index begin at 0)  eg: max_degree = 3, return [0, 3, 8] in python_index
!    !degree = 0, 1, 2  =>   INDEX = 0, 3, 8 means:
!    !g_10 => coef_index(1), g_20 => coef_index(4), g_30 => coef_index(9) in fortran array_index 
!    
!    !:param max_degree: Maximum value of n
!    !:type max_degree: int
!    !:return: List of indexes of zonal coefficients
!    !:rtype: array int
!    !*****************************************************************************************************************
!        integer, intent(in) :: max_degree
!        integer, intent(out) :: n_non_array(:)
!        integer :: i, j, k
!        
!        k = 1
!        do i = 0, (max_degree + 1)**2 - 2
!            do j = 0, max_degree-1
!                if (i == j * (j + 2)) then
!                    goto 200
!                end if
!            end do
!            n_non_array(k) = i
!            k = k + 1
!        200 continue            
!        end do   
!    end subroutine
!!==========================================================================================================================      