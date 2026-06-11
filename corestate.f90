module corestate
    use common
    use utilities
    use config
    use analyser
    use pca
    use blas95
    use lapack95
    use f95_precision
    use omp_lib
    implicit none
    
    type :: CoreState_type
        !*****************************************************************************************************************
        !"""
        !CoreState object. It stores the measure datas in a dict of ndarrays.
        !"""
        !*****************************************************************************************************************        
        type(key_measures), allocatable :: measures_(:)
        type(key_max_degrees), allocatable :: max_degrees_(:)
    contains
        procedure :: init_CoreState, addMeasure
        procedure :: initialise_from_file, initialise_from_noised_priors
        procedure :: cs_Lb, cs_Nb, cs_Lu, cs_Nu2, cs_Lsv, cs_Nsv
    end type
    
    
contains
!==========================================================================================================================  
    subroutine init_CoreState(self, init_measures)
    !*****************************************************************************************************************
    !"""
    !Initiates the CoreState. Initial measures can be given in the init_measures arg:
    !    - Either with only the data (max_degree will be inferred from the data)
    !        Ex: CoreState({SV: np.zeros(224), ...})
    !    - Or by giving also the max_degree as a second member of a 2-tuple/list (the first being the data):
    !        Ex: CoreState({SV: [np.zeros(224), 14], ...})
    !
    !:param init_measures: The measures to add in a dict with members can be np.ndarray or 2-tuple/lists.
    !:type init_measures: dict or None
    !"""
    !*****************************************************************************************************************
        class(CoreState_type), intent(inout) :: self    
        class(key_measures), intent(in) :: init_measures(:)
        integer :: i
        
        allocate(self.measures_, source=init_measures)
        allocate(self.max_degrees_(SIZE(init_measures)))
        do i = 1, SIZE(init_measures)
            self.max_degrees_(i).key = init_measures(i).key
            call self.addMeasure(init_measures(i).key, init_measures(i).measure_data, i)
        end do
        
        
    end subroutine    
!==========================================================================================================================
    
!==========================================================================================================================
    function cs_Lb(self) result(res)
        class(CoreState_type), intent(inout) :: self
        integer :: res
        res = self.max_degrees_(1).max_d
    end function
    
    function cs_Nb(self) result(res)
        class(CoreState_type), intent(inout) :: self
        integer :: res
        res = self.cs_Lb() * (self.cs_Lb() + 2)
    end function
    
    function cs_Lu(self) result(res)
        class(CoreState_type), intent(inout) :: self
        integer :: res
        res = self.max_degrees_(2).max_d
    end function
    
    function cs_Nu2(self) result(res)
        class(CoreState_type), intent(inout) :: self
        integer :: res
        res = self.cs_Lu() * (self.cs_Lu() + 2) * 2
    end function
    
    function cs_Lsv(self) result(res)
        class(CoreState_type), intent(inout) :: self
        integer :: res
        res = self.max_degrees_(3).max_d
    end function
    
    function cs_Nsv(self) result(res)
        class(CoreState_type), intent(inout) :: self
        integer :: res
        res = self.cs_Lsv() * (self.cs_Lsv() + 2)
    end function
!==========================================================================================================================
    
!==========================================================================================================================  
    subroutine addMeasure(self, meas_id, meas_data, idx)
    !*****************************************************************************************************************
    !"""
    !Adds a measure to the CoreState.
    !
    !:param meas_id: name of the measure. Used as key of dict for internal storing.
    !:type meas_id: str
    !:param meas_data: data of the measure.
    !:type meas_data: np.ndarray or list
    !:param idx: the location of key.
    !:type idx: int
    !"""
    !*****************************************************************************************************************
        class(CoreState_type), intent(inout) :: self        
        character(len=*), intent(in) :: meas_id
        real(kind=8), intent(in) :: meas_data(:,:,:)
        integer :: idx
        integer :: nb_coeffs, computed_Lmax
        
        
        !# Try to infer the max degree from last dimension of data N=L(L+2)
        if (trim(meas_id) .ne. 'Z') then
            nb_coeffs = size(meas_data, 3)
            !# If a measure derived from U or S, the equation is N/2 = L(L+2)
            if (('U' == trim(meas_id)) .or. ('dUdt' == trim(meas_id)) .or. ('d2Udt2' == trim(meas_id)) .or. ('S' == trim(meas_id))) then
                if (mod(nb_coeffs, 2) /= 0) then
                    stop
                end if
                nb_coeffs = nb_coeffs / 2
            end if
            computed_Lmax = int(sqrt(real(nb_coeffs + 1)) - 1.0)
            self.max_degrees_(idx).max_d = computed_Lmax
        else
            self.max_degrees_(idx).max_d = -1
        end if        
    end subroutine    
!==========================================================================================================================
    
!==========================================================================================================================  
    subroutine initialise_from_noised_priors(self, random_state, algo_config, algo_attributed_models, algo_avg_prior, algo_cov_prior, algo_analyser, algo_pcaU, Z_AR3)
    !*****************************************************************************************************************
    !"""
    !Initialise the core state Z MF U ER SV at t=0 from the CoreState in a file at a given date.
    !    
    !:param algo: algo instance
    !:type algo: Augkf.algo
    !:param file_path: path of the hdf5 file containing the computed states to use for initialisation
    !:type file_path: str
    !:param date: date of the CoreState to use for the initialisation
    !:type date: datetime64
    !:return: nothing. Simply update self.
    !"""
    !*****************************************************************************************************************
        class(CoreState_type), intent(inout) :: self
        real(kind=8), intent(in) :: random_state
        class(ComputationConfig), intent(in) :: algo_config
        integer, intent(in) :: algo_attributed_models(:)
        class(set_prior_type), intent(in) :: algo_avg_prior
        class(cov_prior_type), intent(in) :: algo_cov_prior
        class(NormedPCAOperator), intent(in) :: algo_pcaU
        class(GenericComputer), intent(in) :: algo_analyser
        real(kind=8), allocatable, intent(out) :: Z_AR3(:,:,:)
        real(kind=8), allocatable :: avg_b(:,:), L_bb(:,:), L_zz(:,:), L_dzdz(:,:), L_d2zd2z(:,:)
        character(len=10) :: AR_type
        integer :: Nz, info, i, j, i_idx
        real(kind=8) :: dt_f
        real(kind=8), allocatable :: w_b(:), w_z(:), w_dz(:,:), w_d2z(:,:)
        real(kind=8), allocatable :: dZ(:,:), d2Z(:,:), AbT(:,:)
        type(input_core_state_type) :: CoreState_temp
        
        !# THE CORESTATE HAS A NUMBER OF REALISATIONS THAT MATCHES THE NUMBER OF ATTRIBUTED MODELS
        !
        !#Averages
        allocate(avg_b, source=algo_avg_prior.MF)
        AR_type = trim(algo_config.AR_type)
        Nz = algo_config.Nz()
        dt_f = algo_config.dt_f
        
        !# Lower Cholesky matrices
        allocate(L_bb, source=algo_cov_prior.B_B)
        allocate(L_zz, source=algo_cov_prior.Z_Z)
        call potrf(L_bb, 'L', info)
        call potrf(L_zz, 'L', info)
        
        if (trim(AR_type) == "AR3") then
            allocate(L_dzdz, source=algo_cov_prior.dZ_dZ)
            allocate(L_d2zd2z, source=algo_cov_prior.d2Z_d2Z)
            call potrf(L_dzdz, 'L', info)
            call potrf(L_d2zd2z, 'L', info)
        end if
           
        !$omp parallel do schedule(static)
        do j = 2, size(L_bb,1)
            L_bb(1:j-1, j) = 0.0d0
        end do
        !$omp end parallel do
        
        !$omp parallel do schedule(static)
        do j = 2, size(L_zz,1)
            L_zz(1:j-1, j) = 0.0d0
            if (trim(AR_type) == "AR3") then
                L_dzdz(1:j-1, j) = 0.0d0
                L_d2zd2z(1:j-1, j) = 0.0d0
            end if           
        end do
        !$omp end parallel do
        
        !# Set random draw
        call random_seed(put=[INT(random_state*50000)])

        !# Loop over attributed models
        CoreState_temp.Lsv = self.cs_Lsv()
        CoreState_temp.Lu = self.cs_Lu()
        CoreState_temp.Lb = self.cs_Lb()
        CoreState_temp.Nsv = self.cs_Nsv()
        CoreState_temp.Nu2 = self.cs_Nu2()
        CoreState_temp.Nb = self.cs_Nb()
        allocate(CoreState_temp.B, source = RESHAPE(avg_b, [SIZE(avg_b)]))
        do i_idx = 1, SIZE(algo_attributed_models)
            !# Set normal draw
            call randn_vec(w_b, algo_config.Nb())
            call randn_vec(w_z, Nz)
            
            !# Initialise B part of core state by normal distrib N(mean_b, sigma_b)
            self.measures_(1).measure_data(i_idx, 1, :) = RESHAPE(avg_b, [SIZE(avg_b)]) + matmul(L_bb, w_b)
            !# Initialise Z part of core state by normal distrib N(0, sigma_z) 
            self.measures_(5).measure_data(i_idx, 1, :) = matmul(L_zz, w_z)
            !# Z to U ER
            call Z_to_U_ER1(algo_config, algo_avg_prior, algo_pcaU, self.measures_(5).measure_data(i_idx, 1, :), self.measures_(2).measure_data(i_idx, 1, :), self.measures_(4).measure_data(i_idx, 1, :))
            !# SV = A(b) U + ER
            CoreState_temp.B = self.measures_(1).measure_data(i_idx, 1, :)
            call algo_analyser.compute_Ab(CoreState_temp, AbT)
            self.measures_(3).measure_data(i_idx, 1, :) = matmul(AbT, self.measures_(2).measure_data(i_idx, 1, :)) + self.measures_(4).measure_data(i_idx, 1, :)
        end do
        
        allocate(Z_AR3(SIZE(algo_attributed_models), 3, Nz), source=0.0d0)
        if (trim(AR_type) == "AR3") then
            !# Compute Z_AR3
            call randn_mat(w_dz, SIZE(algo_attributed_models), Nz)
            call randn_mat(w_d2z, SIZE(algo_attributed_models), Nz)
            !# Set dZ d2Z normal distribution
            allocate(dZ(SIZE(w_dz, 1), SIZE(L_dzdz, 1)), source=0.0d0)
            allocate(d2Z(SIZE(w_d2z, 1), SIZE(L_d2zd2z, 1)), source=0.0d0)
            dZ = matmul(w_dz, TRANSPOSE(L_dzdz))
            d2Z = matmul(w_d2z, TRANSPOSE(L_d2zd2z))
            !# Taylor series of order 2
            Z_AR3(:, 3, :) = self.measures_(5).measure_data(:, 1, :)
            Z_AR3(:, 2, :) = self.measures_(5).measure_data(:, 1, :) - dt_f * dZ  + dt_f**2 / 2.0d0 * d2Z
            Z_AR3(:, 1, :) = self.measures_(5).measure_data(:, 1, :) - (2.0d0*dt_f) * dZ  + (2.0d0*dt_f)**2 / 2.0d0 * d2Z
        end if
    end subroutine    
!==========================================================================================================================
    
!==========================================================================================================================  
    subroutine Z_to_U_ER1(algo_config, algo_avg_prior, algo_pcaU, Z, U, ER)
    !*****************************************************************************************************************
    !"""
    !Compute U ER from Z 
    !
    !:param Z: Augmented state Z
    !:type Z: numpy array, dim Ncoef (1D) or Nreal x Ncoef (2D)
    !:param dimension: dimension of Z
    !:type dimension: int 
    !:return: U and ER states
    !:rtype: 1D or 2D U and ER states
    !"""
    !*****************************************************************************************************************
        class(ComputationConfig), intent(in) :: algo_config
        class(set_prior_type), intent(in) :: algo_avg_prior
        class(NormedPCAOperator), intent(in) :: algo_pcaU
        real(kind=8), intent(in) :: Z(:)
        real(kind=8), intent(inout) :: U(:), ER(:)
        real(kind=8), allocatable :: U_(:,:), Z_(:,:)
        
        allocate(Z_, source = RESHAPE(Z, [1, SIZE(Z)]))
        if (algo_config.pca == 1) then
            ER = Z(algo_config.N_pca_u+1:SIZE(Z)) + RESHAPE(algo_avg_prior.ER, [SIZE(algo_avg_prior.ER)])
            call algo_pcaU.inverse_transform(U_, Z_(:, 1:algo_config.N_pca_u))
            U = RESHAPE(U_, [SIZE(U_)])
        else
            ER = Z(algo_config.Nu2()+1:SIZE(Z)) + RESHAPE(algo_avg_prior.ER, [SIZE(algo_avg_prior.ER)])
            U  = Z(1:algo_config.Nu2()) + RESHAPE(algo_avg_prior.U, [SIZE(algo_avg_prior.U)])
        end if
    end subroutine    
!==========================================================================================================================
    
!==========================================================================================================================  
    subroutine initialise_from_file(self, algo_config, algo_attributed_models, Z_AR3)
    !*****************************************************************************************************************
    !"""
    !Initialise the core state Z MF U ER SV at t=0 from the CoreState in a file at a given date.
    !    
    !:param algo: algo instance
    !:type algo: Augkf.algo
    !:param file_path: path of the hdf5 file containing the computed states to use for initialisation
    !:type file_path: str
    !:param date: date of the CoreState to use for the initialisation
    !:type date: datetime64
    !:return: nothing. Simply update self.
    !"""
    !*****************************************************************************************************************
        class(CoreState_type), intent(inout) :: self        
        class(ComputationConfig), intent(in) :: algo_config
        integer, intent(in) :: algo_attributed_models(:)
        real(kind=8), allocatable, intent(out) :: Z_AR3(:,:,:)
                
        print *, "subroutine initialise_from_file will be added in a future version."
        stop             
    end subroutine    
!==========================================================================================================================
    
!==========================================================================================================================  
    subroutine coef_print(core_state_1D, n_coef, to_print)
    !*****************************************************************************************************************
    !"""
    !Convenience function to print the coef of all core state quantities of a certain index. Note that the index should therefore be lower than the length of the smallest quantity.
    !
    !:param core_state_1D: 1D Core state with all quantities
    !:type core_state_1D: corestates.CoreState
    !:param n_coef: index of the coef to print
    !:type n_coef: int
    !:returns: a string with the 'n_coef'-th coefficients of all measures of the input Corestate
    !:rtype: str
    !"""
    !*****************************************************************************************************************
        class(CoreState_type), intent(in) :: core_state_1D
        integer, intent(in) :: n_coef
        character(len=1000), intent(out) :: to_print
        character(len=100) :: num_str
        integer :: i
                
        to_print = ''
        do i = 1, size(core_state_1D.measures_)
            write(num_str,'(F0.6)') core_state_1D.measures_(i).measure_data(1, 1, n_coef)
            to_print = trim(to_print) // trim(core_state_1D.measures_(i).key) // ': ' // trim(adjustl(num_str)) // ' |  '
        end do
       
    end subroutine    
!==========================================================================================================================
    
    

    
end module
    