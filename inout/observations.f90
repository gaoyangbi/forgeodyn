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
    contains
        procedure :: init_Observation
    end type

contains
!==========================================================================================================================    
    subroutine init_Observation(self)
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
        
    end subroutine    
!==========================================================================================================================
    
!==========================================================================================================================    
    subroutine find_obs_analysis_match()
    !"""
    !return index if obs date found in analysis times +-dt_forecast/2 
    !"""
        
        
    end subroutine    
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
    subroutine build_chaos_hdf5_observations(cfg, nb_realisations, measure_type, seed)
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
        integer, intent(in) :: nb_realisations
        character(len=*), intent(in) :: measure_type
        integer, intent(in) :: seed
        character(len=200) :: datadir, dataset_name
        integer :: max_degree, nb_coefs
                
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
    end subroutine    
!==========================================================================================================================
end module
    