module config
    use utilities
    use reads
    use, intrinsic :: ieee_arithmetic
    use HDF5
    implicit none
    
    type, extends(config_data) :: ComputationConfig
    !*****************************************************************************************************************
    !Defines the configuration of the computation, namely:
    !    - All core state degrees and number of coefs
    !    - Physical constants
    !    - Forecast and Analysis times
    !    - Paths of the prior and observation data
    !    - Output file features
    !*****************************************************************************************************************        
        integer :: do_shear, pca
        real(kind=8) :: dt_a
        real(kind=8), allocatable :: t_forecasts(:), t_analyses(:), t_analyses_full(:)
        character(len=20) :: obs_select
        
    
    contains
        procedure :: init_config
        procedure :: read_from_file
        procedure :: init_bools, init_out, init_physical_constants, init_paths
        procedure :: init_times, init_corestate_inits, init_observations, init_shear
        procedure :: nb_forecasts, nb_analyses, save_hdf5
        procedure :: Nb, Nsv, Nu, Ny, Nu2, Ny2, Nz, Nuz
    
    end type ComputationConfig
    
contains
!==========================================================================================================================
    subroutine init_config(this, do_shear, config_file)
    !*****************************************************************************************************************
    !Loads the config according to parameters. First tries to load config file, then config dict if the config_file is None.
    !
    !:param do_shear: control if we add shear config parameters to config dict
    !:type do_shear: int
    !:param config_file: filename containing the config (see format hereafter)
    !:type config_file: str (default None)
    !
    !config file format:
    !    name type value
    !
    !Example :
    !Lb int 10
    !
    !supported format: int/float/str/bool
    !lines starting by # are not considered.
    !*****************************************************************************************************************
        class(ComputationConfig), intent(inout) :: this
        integer, intent(in) :: do_shear
        character(len=*), intent(in) :: config_file
        
        this.do_shear = do_shear
        call this.read_from_file(TRIM(config_file))
        
        if (this.Lb == 0) then
            write(10,'(A)') "No maximal degree was read for the magnetic field !"
            stop
        end if
        if (this.Lsv == 0) then
            write(10,'(A)') "No maximal degree was read for the secular variation !"
            stop
        end if
        if (this.Lu == 0) then
            write(10,'(A)') "No maximal degree was read for the flow !"
            stop
        end if
        
        call this.init_bools()
        call this.init_out()
        call this.init_physical_constants()
        call this.init_paths()
        call this.init_times()
        call this.init_corestate_inits()
        call this.init_observations()
        
    end subroutine    
!==========================================================================================================================
    function Nb(this) result(res)
        class(ComputationConfig), intent(in) :: this
        integer :: res
        res = this.Lb * (this.Lb + 2)    
    end function
    
    function Nsv(this) result(res)
        class(ComputationConfig), intent(in) :: this
        integer :: res
        res = this.Lsv * (this.Lsv + 2)    
    end function
    
    function Nu(this) result(res)
        class(ComputationConfig), intent(in) :: this
        integer :: res
        res = this.Lu * (this.Lu + 2)    
    end function
    
    function Ny(this) result(res)
        class(ComputationConfig), intent(in) :: this
        integer :: res
        res = this.Ly * (this.Ly + 2)    
    end function
    
    function Nu2(this) result(res)
        class(ComputationConfig), intent(in) :: this
        integer :: res
        res = 2 * this.Nu()    
    end function
    
    function Ny2(this) result(res)
        class(ComputationConfig), intent(in) :: this
        integer :: res
        res = 2 * this.Ny()
    end function
    
    function Nz(this) result(res)
        class(ComputationConfig), intent(in) :: this
        integer :: res
        if (this.pca == 1 ) then
            res = this.N_pca_u + this.Nsv()
        else
            res = this.Nu2() + this.Nsv()
        end if
    end function
        
    function Nuz(this) result(res)
        class(ComputationConfig), intent(in) :: this
        integer :: res
        if (this.pca == 1 ) then
            res = this.N_pca_u
        else
            res = this.Nu2()
        end if
    end function
    
!==========================================================================================================================
    subroutine read_from_file(this, filename)
    !*****************************************************************************************************************
    !""" Reads the config file from a filename and sets the read dictonary as object dict.
    !
    ! :param filename: the file from which we load the data (see format at end)
    ! :type filename: str
    ! """
    !*****************************************************************************************************************
        class(ComputationConfig), intent(inout) :: this
        character(len=*), intent(in) :: filename
        
        call parse_config_file(this, TRIM(filename))    
    end subroutine
!==========================================================================================================================
    
!==========================================================================================================================    
    subroutine init_bools(this)
    !*****************************************************************************************************************
    !"""
    !Init computation booleans.
    !"""
    !*****************************************************************************************************************
        class(ComputationConfig), intent(inout) :: this

        if (this.compute_E == -1) then
            write(10,'(A)') "Compute ER was not set in the config file! Assuming True."
            this.compute_E = 1
        end if
        
        if (TRIM(this.kalman_norm) == 'none') then
            write(10,'(A)') "Kalman_norm was not set in the config file! Assuming l2 norm. string: L2"
            this.kalman_norm = 'l2'
        end if
        
        if (this.N_pca_u == -1) then
            write(10,'(A)') "Number of coefficients for PCA of U was not set in the config file! Assuming no PCA."
            this.N_pca_u = 0
            this.pca = 0
        else if (this.N_pca_u == 0)then
            this.pca = 0
        else
            this.pca = 1
        end if
        
        if (this.remove_spurious == -1.0d0) then
            write(10,'(A)') "Value to discard spurious correlations was not set! Assuming diagonal."
            this.remove_spurious = ieee_value(this.remove_spurious, ieee_positive_inf)
        end if
        
        if (TRIM(this.AR_type) == 'none') then
            write(10,'(A)') "Type of AR process was not set in the config file! Assuming diagonal."
            this.kalman_norm = 'diag'
        end if
        
        if (this.N_pca_u == -1 .AND. TRIM(this.AR_type) == 'AR3') then
            write(10,'(A)') "AR3 process requires N_pca_u to be non-zero"
        end if
        
        if (this.combined_U_ER_forecast == -1) then
            write(10,'(A)') "U ER forecast dependancy was not set in the config file! Assuming independant."
            this.combined_U_ER_forecast = 0
        end if
        
        if (this.combined_U_ER_forecast == 1 .AND. TRIM(this.AR_type) == 'diag') then
            write(10,'(A)') "U ER must be independant for diag forecast"
        end if
        
        !# Checking if we export the data every analysis---------------
        !# Issue #70---------------------------------------------------
        if (this.export_every_analysis == -1) then
            write(10,'(A)') "Data exportation will be done at the end of calculations."
            this.export_every_analysis = 0
        end if
        
        !# Issue #79---------------------------------------------------
        !# Enabling backward analysis or not (defaulting to True)------
        if (this.last_analysis_backward == -1) then
            write(10,'(A)') "Doing last analysis using backward scheme"
            this.last_analysis_backward = 1
        end if        
    end subroutine
!==========================================================================================================================    
    
!==========================================================================================================================    
    subroutine init_out(this)
        class(ComputationConfig), intent(inout) :: this
        
        if (this.out_computed == -1) then
            write(10,'(A)') "Out_computed was not set in the config file ! Assuming computed states saved."
            this.out_computed = 1
        end if
        
        if (this.out_analysis == -1) then
            write(10,'(A)') "Out_analysis was not set in the config file ! Assuming analysis states saved."
            this.out_analysis = 1
        end if
        
        if (this.out_forecast == -1) then
            write(10,'(A)') "Out_forecast was not set in the config file ! Assuming forecast states saved."
            this.out_forecast = 1
        end if
        
        if (this.out_misfits == -1) then
            write(10,'(A)') "Out_misfits was not set in the config file ! Assuming misfits states saved."
            this.out_misfits = 1
        end if
        
        if (TRIM(this.out_format) == 'none') then
            write(10,'(A)') "Out_format was not set in the config file ! Assuming double precision."
            this.out_format = 'float64'
        end if
    
    end subroutine
!==========================================================================================================================    
    
!==========================================================================================================================    
    subroutine init_corestate_inits(this)
        class(ComputationConfig), intent(inout) :: this
        
        if (TRIM(this.core_state_init) == 'none') then
            write(10,'(A)') "Type of CoreState initialisation was not set in the config file! Assuming normal draw around average priors."
            this.core_state_init = 'normal'
        end if
        
        if (TRIM(this.core_state_init) == 'from_file') then
            if (TRIM(this.init_file) == 'none' .OR. this.init_date == -1.0d0) then
                write(10,'(A)') "CoreState initialisation was set to 'from_file' but no file or date were given. Falling back to normal draw initialisation."
                this.core_state_init = 'normal'
            end if
        end if    
    end subroutine
!==========================================================================================================================    
    
!==========================================================================================================================    
    subroutine init_paths(this)
    !*****************************************************************************************************************
    !"""
    !Builds the full path according to the path given in conf or set default.
    !"""
    !*****************************************************************************************************************
        class(ComputationConfig), intent(inout) :: this
        character(len=300) :: paths(3)
        
        if (TRIM(this.prior_dir) == 'none') then
            write(10,'(A)') "Prior directory was not set, assuming 71path."
            this.prior_dir = './data/priors/71path'
            this.prior_type = '71path'
        end if
        
        if (trim(this.prior_type) == '0path') then
            if (TRIM(this.AR_type) .ne. 'diag') then
                write(10,'(A)') 'AR_type must be diag for 0path prior'
            end if
        end if
        
        if (TRIM(this.prior_type) == 'none') then
            write(10,'(A)') 'Prior directory was set ',TRIM(this.prior_dir) ,' without prior type ! Please set the prior type.'
        end if
        
        if (TRIM(this.obs_dir) == 'none') then
            this.obs_dir = './data/observations/COVOBS-x2_maglat'
            this.obs_type = 'COVOBS_hdf5'
        end if    
        
        paths(1) = TRIM(this.prior_dir)
        paths(2) = TRIM(this.obs_dir)
        if (TRIM(this.init_file) /= 'none') then
            paths(3) = TRIM(this.init_file)
        end if
        
        
    end subroutine
!==========================================================================================================================
    
!==========================================================================================================================    
    subroutine init_observations(this)
    !*****************************************************************************************************************
    !"""
    !Builds the full path according to the path given in conf or set default.
    !"""
    !*****************************************************************************************************************
        class(ComputationConfig), intent(inout) :: this
        
        
        if (TRIM(this.obs_type) == 'GO_VO') then
            write(10,'(A)') "default take all obs files in observation directory."
            this.obs_select = 'all'
            if (this.discard_high_lat == -1.0d0) then
                write(10,'(A)') "Suppression of high latitude data was not set, assuming no suppression"
                this.discard_high_lat = 0
            end if
            if (this.SW_err == 'none') then
                write(10,'(A)') 'SW_err was not set, assuming diag error matrix for SWARM'
                this.SW_err = "diag"
            end if          
        end if  
    end subroutine
!==========================================================================================================================
    
!==========================================================================================================================    
    subroutine init_physical_constants(this)
    !*****************************************************************************************************************
    !"""
    !Builds the physical constants of the computation.
    !"""
    !*****************************************************************************************************************
        class(ComputationConfig), intent(inout) :: this
        
        !# Check for theta steps for Legendre polynomial
        if (this.Nth_legendre == -1) then
            write(10,'(A)') "No number of theta steps was read ! Using default value 64."
            this.Nth_legendre = 64
        end if
        
        !# Check the presence of times        
        if (this.TauU == -1) then
            write(10,'(A)') "No characteristic time was read for the core flow ! Using default value 30 yrs(360 mons)."
            this.TauU = 360.0d0
        end if
        
        if (this.TauE == -1) then
            write(10,'(A)') "No characteristic time was read for the subgrid errors ! Using default value 10 yrs(120 mons)."
            this.TauE = 120.0d0
        end if
        
        !# Convert TauU and TauE in decimal values (in years)
        this.TauU = this.TauU / 12.0d0
        this.TauE = this.TauE / 12.0d0
        
    end subroutine
!==========================================================================================================================
    
!==========================================================================================================================    
    subroutine init_times(this)
    !*****************************************************************************************************************
    !""" 
    !Builds the variables linked to the times of computation (forecasts and analysis).
    !"""
    !*****************************************************************************************************************
        class(ComputationConfig), intent(inout) :: this
        real(kind=8) :: t_start_forecast, t_end_forecast
        integer :: i, n
        
        if (this.t_start_analysis == -1.0d0) then
            write(10,'(A)') "No analysis starting time was read !"
            stop
        end if
          
        if (this.t_end_analysis == -1.0d0) then
            write(10,'(A)') "No analysis end time was read !"
            stop
        end if
        
        if (this.dt_f == -1.0d0) then
            write(10,'(A)') "No forecast timestep was read !"
            stop
        end if
        
        if (this.dt_f == 0.0d0) then
            write(10,'(A)') "Forecast timestep cannot be zero !"
            stop
        end if
        
        if (this.dt_a_f_ratio == -1) then
            write(10,'(A)') "No dt_analysis/dt_forecast ratio was read !"
            stop
        end if
        
        this.dt_a = this.dt_f * this.dt_a_f_ratio
        
        if (this.N_dta_start_forecast == -1) then
            this.N_dta_start_forecast = 1
        else if (this.N_dta_start_forecast < 1) then
            write(10,'(A)') "N_dta_start_forecast must be greater or equal to 1"
            stop
        end if
        
        if (this.AR_type == "AR3") then
            if (this.N_dta_end_forecast == -1) then
                this.N_dta_end_forecast = 1
            else
                if (this.N_dta_end_forecast < 1) then
                    write(10,'(A)') "N_dta_end_forecast must be greater or equal to 1"
                    stop
                end if
            end if
        else
            if (this.N_dta_end_forecast == -1) then
                this.N_dta_end_forecast = 0
            end if
        end if
        
        if (this.dt_smoothing == -1.0d0) then
            write(10,'(A)') "No prior dt_smoothing was read ! Using default value 3.2 years."
            this.dt_smoothing = 3.2
        end if
        
        if (this.dt_sampling == -1.0d0) then
            write(10,'(A)') "No prior dt_sampling was read ! Using default value 5.0 years."
            this.dt_sampling = 5.0
        end if
        
        t_start_forecast = this.t_start_analysis - this.dt_a * this.N_dta_start_forecast
        t_end_forecast = this.t_end_analysis + this.dt_a * this.N_dta_end_forecast + this.dt_f
        
        !# Build time arrays
        !# First forecast at t_start + dt_f
        n = ceiling((t_end_forecast - t_start_forecast) / this.dt_f)
        allocate(this.t_forecasts(n))
        i = 0
        do while (t_start_forecast + i*this.dt_f < t_end_forecast)
            this.t_forecasts(i+1) = t_start_forecast + i*this.dt_f
            i = i + 1
        end do
         !# First analysis at t_start
        n = CEILING((this.t_end_analysis + this.dt_f - this.t_start_analysis) / this.dt_a)
        allocate(this.t_analyses(n))
        i = 0
        do while (this.t_start_analysis + i*this.dt_a < this.t_end_analysis + this.dt_f)
            this.t_analyses(i+1) = this.t_start_analysis + i*this.dt_a
            i = i + 1
        end do
        
        if (this.AR_type == "AR3") then
            n = CEILING(((this.t_end_analysis + this.dt_a + this.dt_f) - (this.t_start_analysis - this.dt_a)) / this.dt_a)
            allocate(this.t_analyses_full(n))
            i = 0
            do while (this.t_start_analysis - this.dt_a + i*this.dt_a < this.t_end_analysis + this.dt_a + this.dt_f)
                this.t_analyses_full(i+1) = this.t_start_analysis - this.dt_a + i*this.dt_a
                i = i + 1
            end do
        else
            n = CEILING(((this.t_end_analysis + this.dt_f) - (this.t_start_analysis - this.dt_a)) / this.dt_a)
            allocate(this.t_analyses_full(n))
            i = 0
            do while (this.t_start_analysis - this.dt_a + i*this.dt_a < this.t_end_analysis + this.dt_f)
                this.t_analyses_full(i+1) = this.t_start_analysis - this.dt_a + i*this.dt_a
                i = i + 1
            end do
        end if
    end subroutine
!==========================================================================================================================
    
!==========================================================================================================================    
    subroutine init_shear(this)

        class(ComputationConfig), intent(inout) :: this
        real(kind=8) :: t_start_forecast, t_end_forecast
        integer :: i, n
        
        
    end subroutine
!==========================================================================================================================

!==========================================================================================================================    
    function nb_forecasts(this) result (length)

        class(ComputationConfig), intent(in) :: this
        integer :: length
        
        length = size(this.t_forecasts)        
    end function
!==========================================================================================================================
    
!==========================================================================================================================    
    function nb_analyses(this) result (length)

        class(ComputationConfig), intent(in) :: this
        integer :: length
        
        length = size(this.t_analyses)        
    end function
!==========================================================================================================================

!==========================================================================================================================    
    subroutine save_hdf5(this, hdf5file)

        class(ComputationConfig), intent(inout) :: this
        character(len=*), intent(in) :: hdf5file
        integer :: hdferr, value
        INTEGER(HID_T) :: file_id, attr_id, space_id, root_id
        INTEGER(HSIZE_T), allocatable :: dims(:)
        logical :: attr_exists
        
        ! Initialize FORTRAN interface.
        CALL h5open_f(hdferr)
        
        ! Open file
        CALL h5fopen_f(hdf5file, H5F_ACC_RDWR_F, file_id, hdferr)  
        if (hdferr /= 0) then
            print *, "no HDF5 file"
            stop
        end if
        
        !===============================
        ! check attribute exist
        !===============================
        allocate(dims(1))
        
        !  Lu------------------------------------------------------
        call h5aexists_f(file_id, "Lu", attr_exists, hdferr)
        if (attr_exists) then            
            ! open attribute
            call h5aopen_f(file_id, "Lu", attr_id, hdferr)
            ! write attribute
            value = this.Lu
            call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, value, dims, hdferr)
            call h5aclose_f(attr_id, hdferr)
        else
            ! create dataspace
            call h5screate_f(H5S_SCALAR_F, space_id, hdferr)
            ! create attribute
            call h5acreate_f(file_id, "Lu", H5T_NATIVE_INTEGER, space_id, attr_id, hdferr)
            ! write attribute
            value = this.Lu
            call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, value, dims, hdferr)

            ! ¹Ø±Õ
            call h5aclose_f(attr_id, hdferr)
        end if
        
        !  Lb------------------------------------------------------
        call h5aexists_f(file_id, "Lb", attr_exists, hdferr)
        if (attr_exists) then            
            ! open attribute
            call h5aopen_f(file_id, "Lb", attr_id, hdferr)
            ! write attribute
            value = this.Lb
            call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, value, dims, hdferr)
            call h5aclose_f(attr_id, hdferr)
        else
            ! create dataspace
            call h5screate_f(H5S_SCALAR_F, space_id, hdferr)
            ! create attribute
            call h5acreate_f(file_id, "Lb", H5T_NATIVE_INTEGER, space_id, attr_id, hdferr)
            ! write attribute
            value = this.Lb
            call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, value, dims, hdferr)
        
            ! ¹Ø±Õ
            call h5aclose_f(attr_id, hdferr)
        end if
        
        !  Lsv------------------------------------------------------
        call h5aexists_f(file_id, "Lsv", attr_exists, hdferr)
        if (attr_exists) then            
            ! open attribute
            call h5aopen_f(file_id, "Lsv", attr_id, hdferr)
            ! write attribute
            value = this.Lsv
            call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, value, dims, hdferr)
            call h5aclose_f(attr_id, hdferr)
        else
            ! create dataspace
            call h5screate_f(H5S_SCALAR_F, space_id, hdferr)
            ! create attribute
            call h5acreate_f(file_id, "Lsv", H5T_NATIVE_INTEGER, space_id, attr_id, hdferr)
            ! write attribute
            value = this.Lsv
            call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, value, dims, hdferr)
        
            ! ¹Ø±Õ
            call h5aclose_f(attr_id, hdferr)
        end if
        
        
        !===============================
        ! close file
        !===============================
        call h5fclose_f(file_id, hdferr)

        !===============================
        ! close HDF5
        !===============================
        call h5close_f(hdferr)
        deallocate(dims)
    end subroutine
!==========================================================================================================================
end module
    