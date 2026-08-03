module run
    use corestate
    use augkf_algo
    use config
    use computer
    use observations
    use analyser
    use writes
    implicit none
    private
    !to improve performance in parrallel
    integer :: MKL_NUM_THREADS = 1
    integer :: NUMEXPR_NUM_THREADS = 1
    integer :: OMP_NUM_THREADS = 1
    
    public :: algorithm
    
    
    contains
!==========================================================================================================================
    subroutine send_receive(corestate, comm, rank, attributed_models, corestate_gather)
    !***************************************************************************************************************"""
    !Gather corestate to rank 0 (and broadcast to all ranks if do_bcast=True) 
    !without comm.Gather that had issues with large matrices
    !
    !:param corestate: Corestate to be gathered
    !:type corestate: Cs.Corestate
    !:param comm: MPI communicator
    !:type comm: MPI.comm
    !:param rank: process rank
    !:type rank: int
    !:type attributed_models: 1D numpy array
    !*****************************************************************************************************************"""
        class(CoreState_type), intent(in) :: corestate
        integer, intent(in) :: comm, rank
        integer, intent(in) :: attributed_models(:)
        class(CoreState_type), intent(out) :: corestate_gather
        integer :: ierr, i, j, nprocs, n_rea, n_t, n_coef
        
        call MPI_COMM_SIZE(comm, nprocs, ierr)
        
        
        if (rank == 0) then
            allocate(corestate_gather.measures_, source=corestate.measures_)
            allocate(corestate_gather.max_degrees_, source=corestate.max_degrees_)         
            do j = 1, size(corestate_gather.measures_, 1)
                deallocate(corestate_gather.measures_(j).measure_data)
                n_rea = SIZE(corestate.measures_(j).measure_data, 1) * nprocs
                n_t   = SIZE(corestate.measures_(j).measure_data, 2)
                n_coef = SIZE(corestate.measures_(j).measure_data, 3)
                allocate(corestate_gather.measures_(j).measure_data(n_rea, n_t, n_coef))
            end do
        end if
        
        
    
    end subroutine
!==========================================================================================================================
    
!==========================================================================================================================
    subroutine gather_states(corestate, attributed_models, comm, rank, do_bcast, corestate_gather)
    !***************************************************************************************************************"""
    !Gather corestate to rank 0 (and broadcast to all ranks if do_bcast=True)
    !
    !:param corestate: Corestate to be gathered
    !:type corestate: Cs.Corestate
    !:param attributed_models: models handled by rank process 
    !:type attributed_models: 1D numpy array
    !:param comm: MPI communicator
    !:type comm: MPI.comm
    !:param rank: process rank
    !:type rank: int
    !:param do_bcast: Controls whether the gathered corestate is broadcasted to all process
    !:type do_bcast: boolean
    !*****************************************************************************************************************"""
        class(CoreState_type), intent(in) :: corestate
        integer, intent(in) :: attributed_models(:)
        integer, intent(in) :: comm, rank
        class(CoreState_type), allocatable, intent(out) :: corestate_gather
        logical :: do_bcast
        integer :: ierr, i, j, nprocs, n_rea, n_t, n_coef
        integer :: local_idx, global_idx
        integer, allocatable :: recv_models(:,:)
        real(kind=8), allocatable :: sendbuf(:,:,:)
        real(kind=8), allocatable :: recvbuf(:,:,:,:)
        
        !#synchronyze all processes
        !call MPI_BARRIER(comm, ierr)
        !#gather from all cores to rank 0
        ! In this program, we havn't write the gather code in the order of the realisations but in the order of the measures. 
        ! we have solved this problem! 2026.6.23
        
        !#for each measure in corestate
        call MPI_Comm_size(comm, nprocs, ierr)
        allocate(corestate_gather, source=corestate)
        do i = 1, size(corestate.measures_, 1)
            !#synchronyze all processes
            call MPI_BARRIER(comm, ierr)
            !#gather from all cores to rank 0
            n_rea = SIZE(corestate.measures_(i).measure_data, 1)
            n_t = SIZE(corestate.measures_(i).measure_data, 2)
            n_coef = SIZE(corestate.measures_(i).measure_data, 3)
            
            allocate(sendbuf, source=corestate.measures_(i).measure_data)
            if (rank == 0) then
                allocate(recvbuf(n_rea, n_t, n_coef, nprocs))
                allocate(recv_models(n_rea, nprocs))
            end if
            
            call MPI_GATHER(sendbuf, n_rea*n_t*n_coef, MPI_DOUBLE_PRECISION, &
                  recvbuf, n_rea*n_t*n_coef, MPI_DOUBLE_PRECISION, &
                  0, MPI_COMM_WORLD, ierr)
            
            call MPI_GATHER(attributed_models, n_rea, MPI_INTEGER, &
                  recv_models, n_rea, MPI_INTEGER, &
                  0, MPI_COMM_WORLD, ierr)
            
            deallocate(corestate_gather.measures_(i).measure_data)
            allocate(corestate_gather.measures_(i).measure_data(n_rea*nprocs, n_t, n_coef))
            
            if (rank == 0) then
                do j = 1, nprocs
                    
                    do local_idx = 1, n_rea

                        global_idx = recv_models(local_idx,j) + 1

                        corestate_gather.measures_(i).measure_data(global_idx,:,:) = recvbuf(local_idx,:,:,j)

                    end do                    
                    !corestate_gather.measures_(i).measure_data(n_rea*(j-1)+1:n_rea*j,:,:) = recvbuf(:,:,:,j)
                end do
            end if
            
            !# Broadcast the gathered corestate to all processes if do_bcast is true
            if (do_bcast) then
                !#synchronyze all processes
                call MPI_BARRIER(comm, ierr)
                call MPI_BCAST(corestate_gather.measures_(i).measure_data, n_rea*n_t*n_coef*nprocs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)                
            end if
            
            deallocate(sendbuf)
            if (ALLOCATED(recvbuf)) deallocate(recvbuf)
            if (ALLOCATED(recv_models)) deallocate(recv_models)
        end do
        !print *, "gathering done on rank ", rank, corestate.measures_(1).measure_data(5,1,1), corestate_gather.measures_(1).measure_data(10,1,1), corestate.measures_(1).key
    end subroutine
!==========================================================================================================================
    
!==========================================================================================================================
    subroutine choose_algorithm(algo_name, config_file, nb_models, global_seed, attributed_models, do_shear, return_obj)
        integer, intent(in) :: nb_models, do_shear, global_seed
        integer, intent(in) :: attributed_models(:)
        character(len=*), intent(in) :: config_file, algo_name 
        type(AugkfAlgo), intent(out) :: return_obj
        type(ComputationConfig) :: config
        
        if (trim(algo_name) == 'augkf') then
            write(10,'(A)') 'Using augmented state Kalman filter algorithm'
            write(* ,'(A)') 'Using augmented state Kalman filter algorithm'
            call config.init_config(do_shear, config_file)
            call create_augkf_algo(config, nb_models, global_seed, attributed_models, return_obj)
        else
            print *, 'Algorithm not supported: ', trim(algo_name)
            stop
        end if        
    
    end subroutine
!==========================================================================================================================
    
!==========================================================================================================================
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
        logical status, log_ 
        character(len=300) :: log_path, str, folder_name, output_file
        integer :: i, val, i_t, i_analysis , idx_max, j, ratio
        real(kind=8) :: seed_float, process_rstate, t
        integer, allocatable :: attributed_models(:)
        type(AugkfAlgo) :: algo
        type(CoreState_type) :: computed_states, forecast_states, analysed_states, misfits
        type(CoreState_type) :: computed_states_slice
        class(CoreState_type), allocatable :: analysed_states_slice, forecast_states_slice, gather_computed_states
        class(CoreState_type), allocatable :: computed_states_window
        class(CoreState_type), allocatable :: analysed_states_gather, forecast_states_gather, computed_states_gather
        REAL(kind=8), allocatable :: Z_AR3(:,:,:)
        integer :: R
        
        ! mpi variables---------------------
        integer :: comm, rank, nb_proc, ierr
        integer :: begin_time, end_time, elapsed_time, count_rate
        integer, allocatable :: process_seeds(:)
        integer :: seed_put(1)
        integer :: pseed
        logical :: flag, do_backward_analysis
        integer :: window_start, window_end
        !-----------------------------------
        
        ! test part-------------------------
        
        ! test part-------------------------
        

        
        
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
        folder_name = trim(output_path)//'/'//trim(computation_name)
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
        write(*,'(A, I0, A, I0)') 'Process : ',rank,', seed = ', pseed
        
        ! Initialise the random state using the process seed
        seed_put(1) = pseed
        call RANDOM_SEED(put=seed_put)
        call RANDOM_NUMBER(process_rstate)
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
        call choose_algorithm(algo_name, config_file, nb_models, seed, attributed_models, do_shear, algo)
        !===========================================================================
        
        !# Initialization of states is done in each process using the global seed
        !# Each process has its own attributed models so that there is less transfer of arrays
        call algo.init_corestates(seed, computed_states, forecast_states, analysed_states, misfits, Z_AR3)
        
        if (first_process()) then
            write(10,'(A, *(F7.2,1X))') "Forecast will be performed at following times: ", algo.config.t_forecasts
            write(*,'(A, *(F7.2,1X))') "Forecast will be performed at following times: ", algo.config.t_forecasts
            
            write(10,'(A, *(F7.2,1X))') "Analyses will be performed at following times: ", algo.config.t_analyses
            write(*,'(A, *(F7.2,1X))') "Analyses will be performed at following times: ", algo.config.t_analyses
        end if
        
        !# LAUNCH THE ALGORITHM
        !# INIT
        i_t = 0
        i_analysis = 0
        ratio = algo.config.dt_a_f_ratio
        idx_max = size(algo.config.t_forecasts, 1) - 1
        
        ! # Loop over all time indices
        ! # init the core states types
        allocate(computed_states_slice.measures_, source=computed_states.measures_)
        allocate(computed_states_slice.max_degrees_, source=computed_states.max_degrees_)         
        do j = 1, size(computed_states_slice.measures_, 1)
            deallocate(computed_states_slice.measures_(j).measure_data)
            allocate(computed_states_slice.measures_(j).measure_data(SIZE(computed_states.measures_(j).measure_data, 1), 1, SIZE(computed_states.measures_(j).measure_data, 3)))
        end do
        allocate(forecast_states_slice, source=computed_states_slice)
        allocate(analysed_states_slice, source=computed_states_slice)
        
        do while (i_t < idx_max)
            t = algo.config.t_forecasts(i_t+1)
            
            !# PREPARE
            !# Check if obs data is available and if next analysis on mf and/or sv is performed
            if (trim(algo.config.AR_type) == 'AR3') then
                call algo.analyser_3.check_if_analysis_data(algo.config, i_analysis, .False.)
            else 
                call algo.analyser_1.check_if_analysis_data(algo.config, i_analysis, .False.)
            end if
            
            if (algo.config.last_analysis_backward == 1) then
                if (ABS(t - algo.config.t_analyses(SIZE(algo.config.t_analyses, 1) - 1)) < algo.config.dt_f / 2.0d0) then
                    write(10,'(A)') "Preparing states for backward analysis scheme."
                    write(*,'(A)') "Preparing states for backward analysis scheme."
                    
                    if (trim(algo.config.AR_type) == 'AR3') then
                        call algo.analyser_3.check_if_analysis_data(algo.config, i_analysis, .True.)
                    else 
                        call algo.analyser_1.check_if_analysis_data(algo.config, i_analysis, .True.)
                    end if
                end if
            end if
            
            !# Adapt forecast range R to eventual analysis and AR type
            if ((TRIM(algo.config.AR_type) == "AR3")) then
                if (algo.analyser_3.sv_analysis() .or. algo.analyser_3.mf_analysis()) then
                    R = 2*ratio
                else
                    R = ratio
                end if
            else
                R = ratio
            end if
            
            if (i_t + R > idx_max) then
                R = idx_max - i_t
            end if
            
            !# FORECAST
            do i = 1, R
                !# Increment i_t
                i_t = i_t +1
                t = algo.config.t_forecasts(i_t+1)

                ! # slice computed states for the current time step
                do j = 1, size(computed_states_slice.measures_, 1)
                    computed_states_slice.measures_(j).measure_data = computed_states.measures_(j).measure_data(:, i_t:i_t, :)
                end do
                
                ! # Compute forecast 
                if (trim(algo.config.AR_type) == 'AR3') then
                    call algo.forecaster_3.parallel_forecast_step_AR3( &
                        algo.config,             &
                        algo.nb_realisations,    &
                        algo.attributed_models,  &
                        algo.pcaU_operator,      &
                        algo.avg_prior,          &
                        algo.cov_prior,          &
                        computed_states_slice,   &
                        Z_AR3,                   &
                        seed,                    &
                        i_t,                     &
                        forecast_states_slice    &
                    )
                else
                    call algo.forecaster_1.parallel_forecast_step_AR1(algo.config, algo.nb_realisations, algo.attributed_models, algo.pcaU_operator, algo.avg_prior, algo.cov_prior, computed_states_slice, seed, i_t, forecast_states_slice)
                end if
                
                ! # Update the computed core_state array with the forecast result
                do j = 1, size(computed_states_slice.measures_, 1)
                    computed_states.measures_(j).measure_data(:, (i_t+1):(i_t+1), :) = forecast_states_slice.measures_(j).measure_data
                    forecast_states.measures_(j).measure_data(:, (i_t+1):(i_t+1), :) = forecast_states_slice.measures_(j).measure_data
                end do
            end do
            
            !# ANALYSIS
            if ((trim(algo.config.AR_type) == 'AR3')) then
                if (algo.analyser_3.sv_analysis() .or. algo.analyser_3.mf_analysis()) then
                    ! # Set i_t back from t_a+ to t_a
                    i_t = i_t - ratio
                    t = algo.config.t_forecasts(i_t+1)
                end if
            end if
            
            !# If at least mf or sv analysis
            if (trim(algo.config.AR_type) == 'AR3') then
                log_ = (algo.analyser_3.mf_analysis()) .OR. (algo.analyser_3.sv_analysis())
            else
                log_ = (algo.analyser_1.mf_analysis()) .OR. (algo.analyser_1.sv_analysis())
            end if
            
            if (log_) then
                write(10,'(A, i4, A, F7.2, A)') "Starting analysis #", i_analysis+1, " at time ", t, "..."
                write(*,'(A, i4, A, F7.2, A)') "Starting analysis #", i_analysis+1, " at time ", t, "..."
                
                !# Compute analysis
                if (trim(algo.config.AR_type) == 'AR3') then
                    !# Gather computed_states to get all realisations for analysis
                    !# Issue #79
                    !# Doing backward analysis step for the last analysis
                    !# Avoiding biased forecast
                    do_backward_analysis = .false.
                    if (algo.config.last_analysis_backward == 1) then
                        if (ABS(t - algo.config.t_analyses(SIZE(algo.config.t_analyses))) &
                            < algo.config.dt_f / 2.0d0) then
                            do_backward_analysis = .true.
                            write (10,'(A)') 'Performing an analysis using backward scheme.'
                            write (*,'(A)') 'Performing an analysis using backward scheme.'
                        end if
                    end if

                    if (do_backward_analysis) then
                        window_start = i_t - 2*ratio + 1
                        window_end = i_t + 1
                    else
                        window_start = i_t - ratio + 1
                        window_end = i_t + ratio + 1
                    end if

                    if (window_start < 1 .or. &
                        window_end > SIZE(computed_states.measures_(1).measure_data, 2)) then
                        write (10,'(A,2I8)') 'Invalid AR3 analysis window: ', window_start, window_end
                        write (*,'(A,2I8)') 'Invalid AR3 analysis window: ', window_start, window_end
                        stop
                    end if

                    if (ALLOCATED(computed_states_window)) deallocate(computed_states_window)
                    allocate(computed_states_window, source=computed_states_slice)
                    do j = 1, SIZE(computed_states_window.measures_)
                        deallocate(computed_states_window.measures_(j).measure_data)
                        allocate( &
                            computed_states_window.measures_(j).measure_data( &
                                SIZE(attributed_models),                    &
                                2*ratio+1,                                 &
                                SIZE(computed_states.measures_(j).measure_data, 3) &
                            ),                                             &
                            source=computed_states.measures_(j).measure_data( &
                                :, window_start:window_end, :               &
                            )                                              &
                        )
                    end do

                    call gather_states( &
                        computed_states_window,       &
                        algo.attributed_models,       &
                        comm,                         &
                        rank,                         &
                        .true.,                       &
                        gather_computed_states        &
                    )
                    call algo.analyser_3.analysis_step_AR3( &
                        gather_computed_states,       &
                        Z_AR3,                        &
                        do_backward_analysis,         &
                        algo.config,                  &
                        algo.nb_realisations,         &
                        algo.attributed_models,       &
                        algo.pcaU_operator,           &
                        algo.avg_prior,               &
                        analysed_states_slice         &
                    )

                    do j = 1, SIZE(analysed_states_slice.measures_)
                        analysed_states.measures_(j).measure_data( &
                            :, i_analysis+1:i_analysis+1, :) =      &
                            analysed_states_slice.measures_(j).measure_data
                    end do
                else
                    do j = 1, size(computed_states_slice.measures_, 1)
                        computed_states_slice.measures_(j).measure_data = computed_states.measures_(j).measure_data(:, (i_t+1):(i_t+1), :)
                    end do
                    call gather_states(computed_states_slice, algo.attributed_models, comm, rank, .true., gather_computed_states)
                    call algo.analyser_1.analysis_step(gather_computed_states, algo.config, algo.nb_realisations, algo.attributed_models, algo.pcaU_operator, algo.avg_prior, analysed_states_slice)
                    
                    do j = 1, size(analysed_states_slice.measures_, 1)
                        analysed_states.measures_(j).measure_data(:, (i_analysis+1):(i_analysis+1), :) = analysed_states_slice.measures_(j).measure_data
                    end do
                end if
                
                if (first_process()) then
                    ! write the analysed states to the analysed_states array
                end if
                
                !# Update the computed core_state array with the analysis result (overwrites forecast result)
                do j = 1, size(analysed_states_slice.measures_, 1)
                    computed_states.measures_(j).measure_data(:, (i_t+1):(i_t+1), :) = analysed_states.measures_(j).measure_data(:, (i_analysis+1):(i_analysis+1), :)
                end do
                
                !# Update the misfits
                !# misfits_data(1).key = 'MF'  misfits_data(2).key = 'SV'
                if (trim(algo.config.AR_type) == 'AR3') then
                    misfits.measures_(1).measure_data(:, (i_analysis+1):(i_analysis+1), :) = algo.analyser_3.current_misfits(1)
                    misfits.measures_(2).measure_data(:, (i_analysis+1):(i_analysis+1), :) = algo.analyser_3.current_misfits(2)
                else
                    misfits.measures_(1).measure_data(:, (i_analysis+1):(i_analysis+1), :) = algo.analyser_1.current_misfits(1)
                    misfits.measures_(2).measure_data(:, (i_analysis+1):(i_analysis+1), :) = algo.analyser_1.current_misfits(2)
                end if
                
                ! # logging
                write(10,'(A, i4, A)') "Analysis #", i_analysis+1, " finished !"
                write(*,'(A, i4, A)') "Analysis #", i_analysis+1, " finished !"     
                !print *, 'misfits', algo.analyser_1.mf_X(1).mat(1,1), algo.analyser_1.mf_X(1).mat(1,1)
                
                ! # Issue #70
                ! # Exporting data  
                if (algo.config.export_every_analysis == 1) then
                    write(*,'(A)') "Exporting..."
                    !# TODO: Implement HDF5 output for computed results and misfits.
                    !# Currently, the HDF5 file creation and configuration saving are implemented,
                    !# but numerical results have not yet been written to the output file.
                    if (first_process()) then
                        ! # Creating file 
                    end if      
                end if
                
            end if
            
            if (i_analysis < algo.config.nb_analyses()) then
                if (ABS(t - algo.config.t_analyses(i_analysis+1)) < algo.config.dt_f / 2.0d0) then !# if analysis time
                    ! # Update the analysed core_state array with the computed core_state
                    do j = 1, size(analysed_states_slice.measures_, 1)
                        analysed_states.measures_(j).measure_data(:, (i_analysis+1):(i_analysis+1), :) = computed_states.measures_(j).measure_data(:, (i_t+1):(i_t+1), :)
                    end do
                    ! # Increment i_analysis
                    i_analysis = i_analysis + 1
                end if
                
            end if                        
        end do
        ! # END OF ALGO TIME LOOP
        
        !# START SHEAR COMPUTATION
        
        
        
        
        !# END SHEAR COMPUTATION
        
        !# gather computed_states, analysed_states, forecast_states to rank 0 for saving
        call gather_states(computed_states, algo.attributed_models, comm, rank, .true., computed_states_gather)
        call gather_states(analysed_states, algo.attributed_models, comm, rank, .true., analysed_states_gather)
        call gather_states(forecast_states, algo.attributed_models, comm, rank, .true., forecast_states_gather)
        
        if (first_process()) then
            ! # Create FileHandlers for saving
            output_file = trim(folder_name)//'/'//trim(computation_name)//'.hdf5'
            call algo.config.save_hdf5(TRIM(output_file))
            
            ! # Save requested output groups
            if (algo.config.out_computed == 1) then
                call write_hdf5( &
                    computed_states_gather, algo.config.t_forecasts, output_file, "computed", "float64", &
                    [character(len=8) :: "Z", "S", "dUdt", "dERdt", "d2Udt2", "d2ERdt2"] &
                )
            end if
            if (algo.config.out_forecast == 1) then
                call write_hdf5( &
                    forecast_states_gather, algo.config.t_forecasts, output_file, "forecast", "float64", &
                    [character(len=8) :: "Z", "S", "dUdt", "dERdt", "d2Udt2", "d2ERdt2"] &
                )
            end if
            if (algo.config.out_analysis == 1) then
                call write_hdf5( &
                    analysed_states_gather, algo.config.t_analyses, output_file, "analysed", "float64", &
                    [character(len=1) :: "Z"] &
                )
            end if
            if (algo.config.out_misfits == 1) then
                call write_hdf5(misfits, algo.config.t_analyses, output_file, "misfits", "float64")
            end if
        end if
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
!==========================================================================================================================
    
!==========================================================================================================================
end module
    
