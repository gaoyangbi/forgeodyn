module reads
    use utilities   
    use priors
    use iso_fortran_env
    implicit none
    
    
contains
    subroutine extract_realisations(prior_directory, prior_type, dt_sampling, measures, prior_data_obj)
    !*****************************************************************************************************************
    !Return prior data arrays (nb_realisations x nb_coefs) by dynamically getting the reading function from the prior type.
    !
    !:param prior_directory: directory where the prior data is stored
    !:type prior_directory: str
    !:param prior_type: type of the prior data
    !:type prior_type: str
    !:param dt_sampling: smoothing window duration in second
    !:type dt_sampling: float
    !:param measures: sequence of measures to extract, allowed values are a selection of 'times', 'B', 'U' and 'ER'
    !:return: MF, U, ER, times, dt_samp, tag
    !:rtype: lists
    !*****************************************************************************************************************
        character(len=*), intent(in) :: prior_directory, prior_type
        real(kind=8), intent(in) :: dt_sampling
        character(len=*), intent(in) :: measures(:)
        type(prior_data), intent(inout) :: prior_data_obj(:)
        integer  :: i
    
        ! statically get the function name from the prior type
        if (prior_type == '0path') then
            write (10,'(7A)') 'Reading', TRIM(prior_type)//'  ', (trim(measures(i)), i=1,4), 'data for the priors'
            call read_0_path(prior_directory, prior_type, dt_sampling, measures, prior_data_obj)
        else if (prior_type == '50path' .or. prior_type == '71path') then
            write (10,'(7A)') 'Reading', TRIM(prior_type)//'  ', (trim(measures(i)), i=1,4), 'data for the priors'
            call read_50_and_71_path(prior_directory, prior_type, dt_sampling, measures, prior_data_obj)
        else if (prior_type == '100path') then
            write (10,'(7A)') 'Reading', TRIM(prior_type)//'  ', (trim(measures(i)), i=1,4), 'data for the priors'
            call read_100_path(prior_directory, prior_type, dt_sampling, measures, prior_data_obj)
        else
            write (10,'(A)') 'Error: prior type not recognised, allowed values are 0path, 50path, 71path and 100path'
            stop
        end if  
    end subroutine    
!==========================================================================================================================   
    
!==========================================================================================================================
    subroutine read_analysed_states_hdf5(filepath)
    !*****************************************************************************************************************
    !*****************************************************************************************************************
        integer :: status
        character(len=*), intent(in) :: filepath
        !----------------------------------------------------------
        !Ensure the path directory exists; create it if necessary.
        !set output folder   
        inquire(file=trim(filepath), exist=status) 
        ! only intel fortran have directory option
    
        if (.not. status) then
            write (10,'(A)') 'Output directory does not exist, ! Please check the given path.'
        end if
        !----------------------------------------------------------
    
    
    end subroutine
!==========================================================================================================================  
    
!==========================================================================================================================    
    subroutine parse_config_file(config_obj, filename)
    !*****************************************************************************************************************
    !""" Reads the config file and returns a dictionary of the config. Supported types for reading: int, float, str, bool
    !
    ! :param file_name: file name storing the config
    ! :type file_name: str
    ! :return: the configuration read: dict["parameter"] = value
    ! :rtype: dict in Python, here we set the read values as attributes of the config object
    ! """
    !*****************************************************************************************************************
        character(len=*), intent(in) :: filename
        class(config_data), intent(inout) :: config_obj
        integer :: i, ios, pos
        character(len=200) :: key, dtype, value
        character(len=500) :: line
        logical :: status
    
        !----------------------------------------------------------
        !Ensure the path directory exists; create it if necessary.
        !set output folder   
        inquire(file=trim(filename), exist=status) 
        ! only intel fortran have directory option
        if (.not. status) then
            write (10,'(A)') 'Cannot read configuration : no file path was given ! Please check the given path.'
            stop
        end if
        !----------------------------------------------------------
        
        ! read the file into the config_obj------------------------
        open (30, file=TRIM(filename), status='old')

        do
            read(30, '(A)', iostat=ios) line
            if (ios == iostat_end) exit
            if (line(1:1) == '#' .OR. line(1:1) == ' ') cycle ! skip comment lines           
            
            ! find the first space-----------------------
            pos = index(line, ' ')
            if (pos > 0) then
                key = adjustl(line(1:pos-1))  ! È¥µô×ó¿Õ¸ñ
                line = adjustl(line(pos+1:))    ! Ê£Óà×Ö·û´®
            end if

            ! find the second space-----------------------
            pos = index(line, ' ')
            if (pos > 0) then
                dtype = adjustl(line(1:pos-1))
                value = adjustl(line(pos+1:))
            end if

            select case (TRIM(key))
            case ('Lu')
                read(value, *) config_obj.Lu                
            case ('Lb')
                read(value, *) config_obj.Lb                
            case ('Lsv')
                read(value, *) config_obj.Lsv
            case ('Ly')
                read(value, *) config_obj.Ly
            case ('t_start_analysis')
                read(value, *) config_obj.t_start_analysis
            case ('t_end_analysis')
                read(value, *) config_obj.t_end_analysis
            case ('dt_f')
                read(value, *) config_obj.dt_f
            case ('dt_a_f_ratio')
                read(value, *) config_obj.dt_a_f_ratio
            case ('N_dta_start_forecast')
                read(value, *) config_obj.N_dta_start_forecast
            case ('N_dta_end_forecast')
                read(value, *) config_obj.N_dta_end_forecast                
            case ('Nth_legendre')
                read(value, *) config_obj.Nth_legendre
            case ('TauU')
                read(value, *) config_obj.TauU
                if (TRIM(dtype) == 'years') then
                    config_obj.TauU = config_obj.TauU * 12.0      
                end if    
            case ('TauE')
                read(value, *) config_obj.TauE
                if (TRIM(dtype) == 'years') then
                    config_obj.TauE = config_obj.TauE * 12.0      
                end if
            case ('TauG')
                read(value, *) config_obj.TauG
                if (TRIM(dtype) == 'years') then
                    config_obj.TauG = config_obj.TauG * 12.0      
                end if
            case ('prior_dir')
                config_obj.prior_dir = TRIM(value)
            case ('prior_type')
                config_obj.prior_type = TRIM(value)
            case ('dt_smoothing')
                read(value, *) config_obj.dt_smoothing                
            case ('dt_sampling')
                read(value, *) config_obj.dt_sampling
            case ('compute_e')
                read(value, *) config_obj.compute_e
            case ('obs_dir')
                config_obj.obs_dir = TRIM(value)
            case ('obs_type')
                config_obj.obs_type = TRIM(value)
            case ('AR_type')
                config_obj.AR_type = TRIM(value)
            case ('N_pca_u')
                read(value, *) config_obj.N_pca_u
            case ('pca_norm')
                config_obj.pca_norm = TRIM(value)
            case ('core_state_init')
                config_obj.core_state_init = TRIM(value)
            case ('init_file')
                config_obj.init_file = TRIM(value)
            case ('init_date')
                read(value, *) config_obj.init_date
            case ('out_analysis')
                read(value, *) config_obj.out_analysis
            case ('out_forecast')
                read(value, *) config_obj.out_forecast                
            case ('out_computed')
                read(value, *) config_obj.out_computed
            case ('out_misfits')
                read(value, *) config_obj.out_misfits
            case ('remove_spurious')
                read(value, *) config_obj.remove_spurious
            case ('remove_spurious_shear_u')
                read(value, *) config_obj.remove_spurious_shear_u
            case ('remove_spurious_shear_err')
                read(value, *) config_obj.remove_spurious_shear_err
            case ('combined_U_ER_forecast')
                read(value, *) config_obj.combined_U_ER_forecast
            case ('discard_high_lat')
                read(value, *) config_obj.discard_high_lat
            case ('SW_err')
                config_obj.SW_err = TRIM(value)
            case ('out_format')
                config_obj.out_format = TRIM(value)
            case ('prior_dir_shear')
                config_obj.prior_dir_shear = TRIM(value)
            case ('prior_type_shear')
                config_obj.prior_type_shear = TRIM(value)
            case ('kalman_norm')
                config_obj.kalman_norm = TRIM(value)
            case ('export_every_analysis')
                read(value, *) config_obj.export_every_analysis
            case ('last_analysis_backward')
                read(value, *) config_obj.last_analysis_backward
            end  select
        end do
        close(30)
        !----------------------------------------------------------
        
        
    end subroutine
    
end module
    