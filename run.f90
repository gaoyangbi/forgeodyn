module run
    use corestate
    use augkf_algo
    use config
    implicit none
    private
    !to improve performance in parrallel
    integer :: MKL_NUM_THREADS = 1
    integer :: NUMEXPR_NUM_THREADS = 1
    integer :: OMP_NUM_THREADS = 1
    
    public :: algorithm
    
    
contains
    
    function choose_algorithm(algo_name, config_file, nb_models, global_seed, attributed_models, do_shear) result(return_obj)
        implicit none
        integer, intent(in) :: nb_models, do_shear, global_seed
        integer, intent(in) :: attributed_models(:)
        character(len=100), intent(in) :: config_file, algo_name 
        integer :: return_obj
        
        if (trim(algo_name) == 'augkf') then
            write(10,'(A)') 'Using augmented state Kalman filter algorithm'
            return_obj = 1
        else
            print *, 'Algorithm not supported: ', trim(algo_name)
            stop
        end if
        
    
    end function choose_algorithm

    subroutine algorithm(output_path, computation_name, config_file, nb_models, do_shear, seed, log_file, logging_level, algo_name)
    !***************************************************************************************************************"""
    !Runs the chosen algorithm : takes care of running the forecasts, the analysis, of logging and of saving the data
    !
    !:param output_path: path where the data should be saved
    !:type output_path: str
    !:param computation_name: will create a folder of this name to store the output files
    !:type computation_name: str
    !:param config_file: path to the configuration file
    !:type config_file: str
    !:param nb_models: number of realisations/models to consider
    !:type nb_models: int
    !:param do_shear: control parameter of the shear computation
    !:type do_shear: int
    !:param seed: seed to use for theCs.Corestateel 
    !:param algo_name: name of the algorithm. Supported algorithms is 'augkf'
    !:type algo_name: str
    !:return: CoreStates containing all the results of the computation, the forecasts and analysis.
    !:rtype: CoreState, CoreState, CoreState
    !*****************************************************************************************************************"""
        use mpi
        use config
        use augkf_algo
        implicit none      
        integer, intent(in) :: nb_models, do_shear, seed, logging_level
        character(len=*), intent(in) :: output_path, computation_name, config_file, log_file, algo_name    
        logical status
        character(len=100) :: log_path, str
        integer :: i, val 
        real(kind=8) :: seed_float
        integer, allocatable :: attributed_models(:)
        
        
        ! mpi variables---------------------
        integer :: comm, rank, nb_proc, ierr
        integer :: begin_time, end_time, elapsed_time, count_rate
        integer, allocatable :: process_seeds(:)
        integer :: seed_put(1)
        integer :: pseed
        logical :: flag 
        !-----------------------------------
        
        ! test part-------------------------
        type(AugkfAlgo) :: Augkf_test
        type(ComputationConfig) :: config_test
        type(cov_prior_type) :: cov_prior
        type(set_prior_type) :: avg_prior
        ! test part-------------------------
        
         
        ! type ComputationConfig in config.f90-----------------
        type(ComputationConfig) :: com_config
        !------------------------------------------------------
        
        
        !print *, trim(config_file), trim(algo_name), trim(output_path), trim(computation_name), trim(log_file)
        !print *, nb_models, do_shear, seed, logging_level
    
        
        ! Initialize MPI-------------------------------------------
        call MPI_Initialized(flag, ierr)
        if (.not. flag) then
            call MPI_Init(ierr)
        end if
        
        comm = MPI_COMM_WORLD
        call MPI_Comm_size(comm, nb_proc, ierr)
        call MPI_Comm_rank(comm, rank, ierr)
        !print *, 'MPI initialized: ', flag, ierr
        !print *, 'Number of processes: ', nb_proc, ' Rank: ', rank
        !print *, first_process(),compute_shear()        
        !----------------------------------------------------------
        
        !----------------------------------------------------------
        !Ensure the path directory exists; create it if necessary.
        !set output folder   
        inquire(directory=trim(output_path)//'/'//trim(computation_name), exist=status) 
        ! only intel fortran have directory option
    
        if (.not. status) then
            call execute_command_line('mkdir "'//trim(output_path)//'/'//trim(computation_name)//'"', exitstat=status)            
        end if
        !----------------------------------------------------------
        
        
        !Set log---------------------------------------------------
        write(str, '(I0)') rank
        log_path = TRIM(output_path) // '/' // TRIM(computation_name) // '/' // TRIM(log_file) // TRIM(str) // '.txt'
        print *, 'Logs will be saved in: ', log_path
        !----------------------------------------------------------
        
        
        !Start time------------------------------------------------
        call system_clock(begin_time, count_rate)
        !----------------------------------------------------------
        
        ! INITIALISATION OF ALGO ===================================================
        
        ! set seed-------------------------------------------------
        open(unit=10, file=log_path, status='unknown')
        allocate(process_seeds(nb_proc))
        process_seeds = 0
        if (first_process()) then
            write(10,'(A, I0)') 'Global seed = ', seed
            
            ! Generate a seed for each MPI processes from the global seed
            ! If the same seed was given to each process, then the same random numbers would appear 
            ! on each process
            seed_put(1) = seed
            call RANDOM_SEED(put=seed_put)
            
            ! generate nb_proc random seeds for each process
            do i = 1, nb_proc
                call RANDOM_NUMBER(seed_float)
                process_seeds(i) = int(seed_float * 50000)
            end do
            
        end if
        !----------------------------------------------------------
        
        ! Broadcast the process seed for each MPI process----------
        call MPI_Bcast(seed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(process_seeds, nb_proc, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        pseed = process_seeds(rank+1) ! rank starts at 0, fortran array at 1
        write(10,'(A, I0, A, I0)') 'Process : ',rank,', seed = ', pseed
        
        ! Initialise the random state using the process seed
        seed_put(1) = pseed
        call RANDOM_SEED(put=seed_put)
        !----------------------------------------------------------
        
        
        ! Build the attribution of realisations if several processes in parallel----
        if (MOD(nb_models, nb_proc) > rank) then 
            allocate(attributed_models(nb_models/nb_proc + 1))
        else
            allocate(attributed_models(nb_models/nb_proc))
        end if
         
        val = rank
        do i = 1, SIZE(attributed_models)
            attributed_models(i) = val
            val = val + nb_proc
        end do        

        write(10,'(A,I0,A)', advance='no') 'Process ', rank, ' will process realisations: '
        do i = 1, size(attributed_models)
            write(10,'(I0,1X)', advance='no') attributed_models(i)
        end do 
        write(10,*) ''
        !---------------------------------------------------------------------------
        
        ! Set algo
        !print *,  choose_algorithm(algo_name, config_file, nb_models, pseed, attributed_models, do_shear)
        !call com_config.init_config(do_shear, config_file)
        print *, attributed_models
        !===========================================================================
        
        
        !test parts---------------------------------------------------------------------
        print *, "test parts run--------------------------------------------"
        call config_test.init_config(0, 'D:/VS/program_Fortran/pygeodyn_fortran/pygeodyn_fortran/code_use.conf')
        call config_test.save_hdf5('D:\VS\program_Fortran\pygeodyn_fortran\test.hdf5')
        call Augkf_test.init_AugkfAlgo(config_test, 500, 500, attributed_models)
        call Augkf_test.extract_prior_and_covariances(avg_prior, cov_prior)
        print *, cov_prior.A(2,:)
        print *, SIZE(cov_prior.A, 1), SIZE(cov_prior.A, 2)
        print *, "test parts run--------------------------------------------"
        !---------------------------------------------------------------------------
        
        
        
        !----------------------------------------------------------
        call MPI_FINALIZE (ierr)
        deallocate(process_seeds, attributed_models)
        close(10)
        
    contains
        logical function first_process()
            first_process = (rank == 0)
        end function first_process
        
        logical function compute_shear()
            compute_shear = (do_shear == 1)
        end function compute_shear
        
        
    end subroutine algorithm
end module
    