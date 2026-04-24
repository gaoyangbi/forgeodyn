!  pygeodyn_fortran.f90 
!
!  FUNCTIONS:
!  pygeodyn_fortran - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: pygeodyn_fortran
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program geodyn
    
    use utilities
    use run
    use priors
    use reads
    use config
    use augkf_algo
    implicit none

    ! Variables
    character(len=1000) :: arg !command line arfuments
    integer :: i, num_args, arg_length, status
    integer :: m = 20
    integer :: shear = 0
    integer :: seed = 10
    integer :: d = 2    
    integer :: ios, skip_count
    
    character(len=100) :: conf = 'con_file.txt'
    character(len=100) :: algo = 'augkf'
    character(len=100) :: path = './pygeodyn_results'
    character(len=100) :: cname = 'Current_computation'
    character(len=100) :: l = 'log'
    
    logical :: conf_flag, m_flag, algo_flag, shear_flag, seed_flag, d_flag, path_flag, cname_flag, l_flag
    

    !----------test parts
    type(obs_type) :: obs_type_to_use(3)
    integer :: n(3)
    character(len=20) :: name(3)=['a', 'b', 'c']
    character(len=20) :: name_(3)=['1', '2', '3']
    
    integer*4 n_array(5), P(30,30)
    character(len=20) :: c, tag
    real*8 coef(261,3)
    logical :: n_log_array(70)
    integer :: row, column
    character(len=100) :: measures(4)
    real(kind=8), allocatable :: MF(:,:), U(:,:), ER(:,:), times(:), dt_samp(:)
    type(prior_data) :: aa_test(10)
    type(AugkfAlgo) :: Augkf_test
    type(ComputationConfig) :: config_test
    logical :: exis
    !--------------------------
    
    
    
    ! Body of pygeodyn_fortran
    conf_flag = .false.
    m_flag    = .false.
    algo_flag = .false.
    shear_flag= .false.
    seed_flag = .false.
    d_flag    = .false.
    path_flag = .false.
    cname_flag= .false.
    l_flag    = .false.
    !------------------------------get the command_argument
    num_args = command_argument_count()
     
    if (num_args == 0) then
        call print_usage()
        stop
    end if  

    
    !----------------------read the num_args
    do i = 1, num_args
        if (skip_count > 0) then
            skip_count = skip_count - 1
            cycle   ! 
        end if
        
        call get_command_argument(i, arg, arg_length, status)

        select case (trim(arg))
            
            case ("-h", "--help")
                call print_help()
                stop
                
            case ("-v", "--version")
                call print_version()
                stop
                
            case ("-conf")
                if (i == num_args) then
                    print *, "Error: -conf requires a value"
                    stop
                end if
                call get_command_argument(i+1, arg)
                read(arg, '(A)', iostat=ios) conf
                if (ios /= 0) then
                    print *, "Error: invalid integer for -conf"
                    stop
                end if
                skip_count = 1
                
            case ("-m")
                if (i == num_args) then
                    print *, "Error: -m requires a value"
                    stop
                end if
                call get_command_argument(i+1, arg)
                read(arg, *, iostat=ios) m
                if (ios /= 0) then
                    print *, "Error: invalid integer for -m"
                    stop
                end if
                skip_count = 1
                
                
            case ("-algo")
                if (i == num_args) then
                    print *, "Error: -algo requires a value"
                    stop
                end if
                call get_command_argument(i+1, arg)
                read(arg, *, iostat=ios) algo
                if (ios /= 0) then
                    print *, "Error: invalid integer for -algo"
                    stop
                end if
                skip_count = 1
                
            case ("-shear")
                if (i == num_args) then
                    print *, "Error: -shear requires a value"
                    stop 
                end if
                call get_command_argument(i+1, arg)
                read(arg, *, iostat=ios) shear
                if (ios /= 0) then
                    print *, "Error: invalid value for -shear"
                    stop
                end if
                skip_count = 1
                
            case ("-seed")
                if (i == num_args) then
                    print *, "Error: -seed requires a value"
                    stop 
                end if
                call get_command_argument(i+1, arg)
                read(arg, *, iostat=ios) seed
                if (ios /= 0) then
                    print *, "Error: invalid value for -seed"
                    stop 
                end if
                skip_count = 1
                
            case ("-d")
                if (i == num_args) then
                    print *, "Error: -d requires a value"
                    stop
                end if
                call get_command_argument(i+1, arg)
                read(arg, *, iostat=ios) d
                if (ios /= 0) then
                    print *, "Error: invalid value for -d"
                    stop
                end if
                skip_count = 1
                
            case ("-path")
                if (i == num_args) then
                    print *, "Error: -path requires a value"
                    stop
                end if
                call get_command_argument(i+1, arg)
                path = TRIM(arg)
                if (ios /= 0) then
                    print *, "Error: invalid value for -path"
                    stop
                end if
                skip_count = 1
                
            case ("-cname")
                if (i == num_args) then
                    print *, "Error: -cname requires a value"
                    stop
                end if
                call get_command_argument(i+1, arg)
                read(arg, *, iostat=ios) cname
                if (ios /= 0) then
                    print *, "Error: invalid value for -cname"
                    stop
                end if
                skip_count = 1
                
            case ("-l")
                if (i == num_args) then
                    print *, "Error: -l requires a value"
                    stop
                end if
                call get_command_argument(i+1, arg)
                read(arg, *, iostat=ios) l
                if (ios /= 0) then
                    print *, "Error: invalid value for -l"
                    stop
                end if
                skip_count = 1
                
            case default
                print *, "Unknown option: ", trim(arg)
                stop
                
        end select
        
    end do
    
    
    !----------------------code draft
    !do i = 1, num_args
    !    call get_command_argument(i, arg, arg_length, status)        
    !    if (status == 0) then
    !        !  -h  --help
    !        if (trim(arg) == "-h" .or. trim(arg) == "--help") then
    !            call print_help()
    !            cycle
    !        !  -v  --version
    !        else if (trim(arg) == "-v" .or. trim(arg) == "--version") then
    !            call print_version()
    !            cycle
    !        else if (trim(arg) == "-conf") then
    !            conf_flag = .true.
    !            cycle
    !        else if (trim(arg) == "-m") then
    !            m_flag = .true.
    !            cycle
    !        else if (trim(arg) == "-algo") then
    !            algo_flag = .true.
    !            cycle
    !        else if (trim(arg) == "-shear") then
    !            shear_flag = .true.
    !            cycle
    !        else if (trim(arg) == "-seed") then
    !            seed_flag = .true.
    !            cycle
    !        else if (trim(arg) == "-d") then
    !            d_flag = .true.
    !            cycle
    !        else if (trim(arg) == "-path") then
    !            path_flag = .true.
    !            cycle
    !        else if (trim(arg) == "-cname") then
    !            cname_flag = .true.
    !            cycle
    !        else if (trim(arg) == "-l") then
    !            l_flag = .true.
    !            cycle
    !        else
    !            if (conf_flag) then
    !                conf = trim(arg)
    !                conf_flag = .false.
    !            else if (m_flag) then
    !                read(arg, *) m
    !                m_flag = .false. 
    !            else if (algo_flag) then
    !                algo = trim(arg)
    !                algo_flag = .false.    
    !            else if (shear_flag) then
    !                read(arg, *) shear
    !                shear_flag = .false. 
    !            else if (seed_flag) then
    !                read(arg, *) seed
    !                seed_flag = .false.
    !            else if (d_flag) then
    !                read(arg, *) d
    !                d_flag = .false. 
    !            else if (path_flag) then
    !                path = trim(arg)
    !                path_flag = .false.    
    !            else if (cname_flag) then
    !                cname = trim(arg)
    !                cname_flag = .false. 
    !            else if (l_flag) then
    !                l = trim(arg)
    !                l_flag = .false.
    !            else
    !                print *, "Error: Unknown option ", trim(arg)
    !                call exit()
    !            endif
    !        end if
    !    else
    !        print *, "Error retrieving argument ", i
    !        call exit()
    !    end if
    !    
    !end do    

    call algorithm(path, cname, conf, m, shear, seed, l, d, algo)
    
    
    !
    !print *, conf
    !print *, m
    !print *, algo
    !print *, shear
    !print *, seed
    !print *, d
    !print *, path
    !print *, cname
    !print *, l
    
    
    
    !----------test parts
    print *, "test parts--------------------------------------------"
    
    measures = ['times', 'B', 'U', 'ER']
    
    call extract_realisations('d:\vs\program_fortran\pygeodyn_fortran\pygeodyn_fortran\data\priors\100path', '100path', 0.5d0, measures, aa_test)
    !print *, aa_test.dim_times, aa_test.dim_MF_1, aa_test.dim_MF_2, aa_test.dim_U_1, aa_test.dim_U_2, aa_test.dim_ER_1, aa_test.dim_ER_2 
    !print *, aa_test.MF(1, :)
    !print *, aa_test.tag
    !print *, aa_test.dt_samp
    print *, size(aa_test(2).U, 1)
    !write (*,'(4A)') (trim(measures(i))// '  ', i=1,4) 
    
    !
    !n = 10
    !obs_type_to_use.seed = 20
    !obs_type_to_use.num_obs = name_
    !obs_type_to_use.measure_name = name
    !seed = 56
    !print *, obs_type_to_use.seed
    !call get_seeds_for_obs(seed, obs_type_to_use)
    !print *, obs_type_to_use.seed

    call config_test.init_config(0, 'D:/VS/program_Fortran/pygeodyn_fortran/pygeodyn_fortran/code_use.conf')
    call config_test.save_hdf5('D:\VS\program_Fortran\pygeodyn_fortran\test.hdf5')
    print *, 21 / 2
    print *, config_test.Nuz()
    

    
    print *, "test parts--------------------------------------------"
    !----------------------
    end program geodyn

    
