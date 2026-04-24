module pca
!==================================
!pca_module
!In python, index begin at 0
!In fortran, index begin at 1
!==================================
    use blas95
    use lapack95
    use f95_precision
    use utilities
    use config
    implicit none
    public
    
    type :: NormedPCAOperator
    !"""
    !Object handling the transforms with normalisation of the PCA of U.
    !Generates the sklearn operator (used to fit the normed data) and the norm matrix used for normalisation.
    !Implements the fit_transform, transform and inverse_transform methods with normalisation taken into account.
    !"""
        real(kind=8), allocatable :: norm_matrix(:,:), inverse_norm_matrix(:,:)
        real(kind=8), allocatable :: sklearn_operator_components(:,:), sklearn_operator_mean(:,:)
        class(ComputationConfig), allocatable :: config

    contains
        procedure :: init_PCAOperator, fit, fit_
        procedure :: S_u, inv_S_u, U0, n_components
        procedure :: transform, inverse_transform, transform_deriv, inverse_transform_deriv
    end type
    
    
    type, extends(NormedPCAOperator) :: NormedNzPCAOperator
        
    contains
        procedure :: init_NzPCAOperator
    end type
    
    
    
contains
!========================================================================================================================== 
    subroutine init_PCAOperator(self, config)
    !***************************************************************************************************************
    !"""
    !Object handling the transforms with normalisation of the PCA of U.
    !Generates the sklearn operator (used to fit the normed data) and the norm matrix used for normalisation.
    !Implements the fit_transform, transform and inverse_transform methods with normalisation taken into account.
    !"""
    !****************************************************************************************************************
        class(NormedPCAOperator), intent(inout) :: self
        class(ComputationConfig), intent(in) :: config
        self.config = config
        call compute_normalisation_matrix(config.Lu, config.pca_norm, self.norm_matrix, self.inverse_norm_matrix)        
        !# Assuming that the sklearn_operator was used to do the PCA on the normed data Û = NM*U (NM = norm_matrix):
        !# we have U_pca = SKC*(Û - Û0) with SKC = sklearn_operator.components_ and Û0 = sklearn_operator.mean_
        !# so Û = SKC^T * U_pca + Û0
        !# => U = NM^(-1) SKC^T * U_pca + NM^(-1) * Û0
        !# => U = S_u * U_pca + U0
        !# => U_pca = inv_S_u * (U - U0)
        !# /!\ S_u is no longer orthogonal due to normalisation so the pseudo inverse need to be stored
    end subroutine
!==========================================================================================================================  

!========================================================================================================================== 
    subroutine init_NzPCAOperator(self, config)
    !***************************************************************************************************************
    !"""
    !Object handling the transforms with normalisation of the PCA of U.
    !Generates the sklearn operator (used to fit the normed data) and the norm matrix used for normalisation.
    !Implements the fit_transform, transform and inverse_transform methods with normalisation taken into account.
    !"""
    !****************************************************************************************************************
        class(NormedNzPCAOperator), intent(inout) :: self
        class(ComputationConfig), intent(in) :: config
        integer, allocatable :: P_U_Ua(:, :)
        real(kind=8), allocatable :: matrix_1(:,:), matrix_2(:,:)
        
        call self.init_PCAOperator(config)
        ! # Project into zU,nzU space to recover only the nz part of norm_matrix
        call spectral_to_znz_matrix(self.config.Lu, P_U_Ua)
        
        allocate(matrix_1(SIZE(self.norm_matrix,1), SIZE(self.norm_matrix,2)))
        allocate(matrix_2(SIZE(self.norm_matrix,1), SIZE(self.norm_matrix,2)))
        
        matrix_1 = MATMUL(MATMUL(P_U_Ua, self.norm_matrix), TRANSPOSE(P_U_Ua))
        matrix_2 = MATMUL(MATMUL(P_U_Ua, self.inverse_norm_matrix), TRANSPOSE(P_U_Ua))
        if (ALLOCATED(self.norm_matrix)) deallocate(self.norm_matrix)
        if (ALLOCATED(self.inverse_norm_matrix)) deallocate(self.inverse_norm_matrix)
        
        allocate(self.norm_matrix(SIZE(matrix_1,1)-self.config.Lu, SIZE(matrix_1,2)-self.config.Lu))
        allocate(self.inverse_norm_matrix(SIZE(matrix_1,1)-self.config.Lu, SIZE(matrix_1,2)-self.config.Lu))

        self.norm_matrix = matrix_1(self.config.Lu+1:, self.config.Lu+1:)
        self.inverse_norm_matrix = matrix_2(self.config.Lu+1:, self.config.Lu+1:)
    end subroutine
!==========================================================================================================================  
    
!==========================================================================================================================
    subroutine fit(self, matrix)
        class(NormedPCAOperator), intent(inout) :: self
        real(kind=8), intent(in) :: matrix(:,:)

        call self.fit_(MATMUL(matrix, self.norm_matrix))
    end subroutine    
!==========================================================================================================================  
    
!==========================================================================================================================  
    subroutine fit_(self, realisations_U)
        class(NormedPCAOperator), intent(inout) :: self
        real(kind=8), intent(in) :: realisations_U(:,:)
        real(kind=8), allocatable :: mean_X(:, :), l(:, :)
        real(kind=8), allocatable :: Y(:, :)
        real(kind=8) :: alpha
        real(kind=8), allocatable :: s(:)
        real(kind=8), allocatable :: u(:, :), vt(:, :)
        integer :: i, idx(1)
        
        allocate(mean_X(1, SIZE(realisations_U, 2)))
        allocate(l(1, SIZE(realisations_U, 1)))
        allocate(Y(SIZE(realisations_U, 1), SIZE(realisations_U, 2)))
        allocate(s(MIN(SIZE(realisations_U, 1), SIZE(realisations_U, 2))))
        allocate(u(SIZE(realisations_U, 1), SIZE(realisations_U, 1)))
        allocate(vt(SIZE(realisations_U, 2), SIZE(realisations_U, 2)))
        !-----------------------juping remove mean value
        alpha = 1.0d0 / real(SIZE(realisations_U, 1), kind=8)
        l     = 1.0d0
        mean_X = 0.0d0
        Y     = realisations_U
        call gemm(l, realisations_U, mean_X, 'N', 'N', alpha)
        call gemm(l*(-1.0d0), mean_X, Y, 'T', 'N', 1.0d0, 1.0d0)
        
        !-----------------------SVD

        call gesvd(Y, s, u, vt)
        
        if (ALLOCATED(self.sklearn_operator_components)) deallocate(self.sklearn_operator_components)
        if (ALLOCATED(self.sklearn_operator_mean)) deallocate(self.sklearn_operator_mean)
        
        allocate(self.sklearn_operator_components(self.config.N_pca_u, SIZE(realisations_U, 2)))
        allocate(self.sklearn_operator_mean(1, SIZE(realisations_U, 2)))
        
        !do i = 1, self.config.N_pca_u
        !    u(:,i) = u(:,i) * s(i)
        !end do
        !
        self.sklearn_operator_components = vt(1:self.config.N_pca_u,:)
        self.sklearn_operator_mean = mean_X
        
        do i = 1, self.config.N_pca_u
            idx = maxloc(abs(self.sklearn_operator_components(i,:)))
            if (self.sklearn_operator_components(i,idx(1)) <=0 ) self.sklearn_operator_components(i, :) = self.sklearn_operator_components(i, :) * (-1.0d0)
        end do       

    end subroutine
!==========================================================================================================================

!==========================================================================================================================   
    subroutine S_u(self, result_matrix)
        class(NormedPCAOperator), intent(in) :: self
        real(kind=8), allocatable, intent(out) :: result_matrix(:,:)
        
        if (ALLOCATED(result_matrix)) deallocate(result_matrix)
        allocate(result_matrix(SIZE(self.inverse_norm_matrix, 1),self.config.N_pca_u))
        result_matrix = 0.0d0
        call gemm(self.inverse_norm_matrix, self.sklearn_operator_components, result_matrix, 'N', 'T')
    end subroutine
    
    subroutine inv_S_u(self, result_matrix)
        class(NormedPCAOperator), intent(in) :: self
        real(kind=8), allocatable, intent(out) :: result_matrix(:,:)
        
        if (ALLOCATED(result_matrix)) deallocate(result_matrix)
        allocate(result_matrix(SIZE(self.sklearn_operator_components, 1), SIZE(self.norm_matrix, 2)))
        result_matrix = 0.0d0
        call gemm(self.sklearn_operator_components, self.norm_matrix, result_matrix, 'N', 'N')
    end subroutine
    
    subroutine U0(self, result_matrix)
        class(NormedPCAOperator), intent(in) :: self
        real(kind=8), allocatable, intent(out) :: result_matrix(:,:)
        
        if (ALLOCATED(result_matrix)) deallocate(result_matrix)
        allocate(result_matrix(SIZE(self.sklearn_operator_mean, 1), SIZE(self.inverse_norm_matrix, 1)))
        result_matrix = 0.0d0
        call gemm(self.sklearn_operator_mean, self.inverse_norm_matrix, result_matrix, 'N', 'T')
    end subroutine
    
    subroutine n_components(self, result_)
        class(NormedPCAOperator), intent(in) :: self
        integer, intent(out):: result_
        
        result_ = self.config.N_pca_u
    end subroutine
    
    subroutine transform(self, result_matrix, U)
        class(NormedPCAOperator), intent(in) :: self
        real(kind=8), allocatable, intent(out) :: result_matrix(:,:)
        real(kind=8), intent(in) :: U(:,:) ! row vector
        real(kind=8), allocatable :: matrix_inv_S_u(:,:), matrix_U0(:,:), matrix(:,:)
        integer :: i
        
        if (ALLOCATED(result_matrix)) deallocate(result_matrix)
        call self.inv_S_u(matrix_inv_S_u)
        call self.U0(matrix_U0)
        allocate(result_matrix(SIZE(U, 1), SIZE(matrix_inv_S_u, 1)))
        allocate(matrix(SIZE(U, 1), SIZE(U, 2)))
        do concurrent (i = 1: SIZE(U, 1))
             matrix(i,:) = U(i,:) - matrix_U0(1,:)
        end do
        result_matrix = 0.0d0
        call gemm(matrix, matrix_inv_S_u, result_matrix, 'N', 'T')
    end subroutine
    
    subroutine inverse_transform(self, result_matrix, pcaU)
        class(NormedPCAOperator), intent(in) :: self
        real(kind=8), allocatable, intent(out) :: result_matrix(:,:)
        real(kind=8), intent(in) :: pcaU(:,:) ! row vector
        real(kind=8), allocatable :: matrix_S_u(:,:), matrix_U0(:,:)
        integer :: i
        
        if (ALLOCATED(result_matrix)) deallocate(result_matrix)
        call self.S_u(matrix_S_u)
        call self.U0(matrix_U0)
        allocate(result_matrix(SIZE(pcaU, 1), SIZE(matrix_S_u, 1)))
        result_matrix = 0.0d0
        call gemm(pcaU, matrix_S_u, result_matrix, 'N', 'T')
        
        do concurrent (i = 1: SIZE(result_matrix, 1))
             result_matrix(i,:) = result_matrix(i,:) + matrix_U0(1, :)
        end do
        
    end subroutine
    
    subroutine transform_deriv(self, result_matrix, dU)
        class(NormedPCAOperator), intent(in) :: self
        real(kind=8), allocatable, intent(out) :: result_matrix(:,:)
        real(kind=8), intent(in) :: dU(:,:) ! row vector
        real(kind=8), allocatable :: matrix_inv_S_u(:,:)
        
        if (ALLOCATED(result_matrix)) deallocate(result_matrix)
        call self.inv_S_u(matrix_inv_S_u)
        allocate(result_matrix(SIZE(dU, 1), SIZE(matrix_inv_S_u, 1)))
        result_matrix = 0.0d0
        call gemm(dU, matrix_inv_S_u, result_matrix, 'N', 'T')
    end subroutine
    
    subroutine inverse_transform_deriv(self, result_matrix, pcadU)
        class(NormedPCAOperator), intent(in) :: self
        real(kind=8), allocatable, intent(out) :: result_matrix(:,:)
        real(kind=8), intent(in) :: pcadU(:,:)
        real(kind=8), allocatable :: matrix_S_u(:,:)
        
        if (ALLOCATED(result_matrix)) deallocate(result_matrix)
        call self.S_u(matrix_S_u)
        allocate(result_matrix(SIZE(pcadU, 1), SIZE(matrix_S_u, 1)))
        result_matrix = 0.0d0
        call gemm(pcadU, matrix_S_u, result_matrix, 'N', 'T')
    end subroutine
!==========================================================================================================================  
    
!==========================================================================================================================  
    subroutine compute_normalisation_matrix(Lu, norm_type, norm_matrix, inverse_norm_matrix)
    !***************************************************************************************************************
    !"""
    !Computes the normalisation matrix for the PCA according to the asked norm type.
    !Default is None: the normalisation matrix (and its inverse) is then the identity.
    !
    !:param Lu: max_degree of core flow U
    !:type Lu: int
    !:param norm_type: Normalisation type
    !:type norm_type: "energy" or None
    !:return: Normalisation matrix and its inverse
    !:rtype: np.array (Nu2 x Nu2), np.array (Nu2 x Nu2)
    !"""
    !****************************************************************************************************************
        integer, intent(in) :: Lu
        character(len=*), intent(in) :: norm_type
        real(kind=8), allocatable, intent(out) :: norm_matrix(:,:), inverse_norm_matrix(:,:)
        integer :: i 
        
        if (ALLOCATED(norm_matrix)) deallocate(norm_matrix)
        if (ALLOCATED(inverse_norm_matrix)) deallocate(inverse_norm_matrix)
        if (TRIM(norm_type) == 'energy') then
            write(10, '(A)') "Computing the normalisation matrix for the PCA with energy norm..."
            call compute_energy_normalisation_matrix(Lu, norm_matrix, inverse_norm_matrix)
        else
            !# No normalisation
            write(10, '(A)') "Computing the normalisation matrix for the PCA without any normalisation..."
            allocate(norm_matrix(2*Lu*(Lu+2), 2*Lu*(Lu+2)))
            allocate(inverse_norm_matrix(2*Lu*(Lu+2), 2*Lu*(Lu+2)))
            
            norm_matrix = 0.0d0
            inverse_norm_matrix = 0.0d0
            do concurrent (i = 1:2*Lu*(Lu+2))
                norm_matrix(i, i) = 1.0d0
                inverse_norm_matrix(i, i) = 1.0d0
            end do
        end if
    
    end subroutine
!========================================================================================================================== 
    
!==========================================================================================================================  
    subroutine compute_energy_normalisation_matrix(Lu, norm_matrix, inverse_norm_matrix)
    !***************************************************************************************************************"""
    !"""
    !Computes the matrix normalising coefficients with respect to energy.
    !That is: norm_tnm = sqrt(n*(n + 1)/(2 * n + 1))*tnm
    !
    !:param Lu: max_degree of core flow U
    !:type Lu: int
    !:return: Energy normalisation matrix and its inverse
    !:rtype: np.array (Nu2 x Nu2), np.array (Nu2 x Nu2)
    !"""
    !*****************************************************************************************************************"""
        integer, intent(in) :: Lu
        real(kind=8), allocatable, intent(out) :: norm_matrix(:,:), inverse_norm_matrix(:,:)
        integer :: i, Nu, min_index, next_index, j
        real(kind=8) :: norm_factor
        
        if (ALLOCATED(norm_matrix)) deallocate(norm_matrix)
        if (ALLOCATED(inverse_norm_matrix)) deallocate(inverse_norm_matrix)
        Nu = Lu*(Lu+2)
        allocate(norm_matrix(2*Nu, 2*Nu))
        allocate(inverse_norm_matrix(2*Nu, 2*Nu))
        norm_matrix = 0.0d0
        inverse_norm_matrix = 0.0d0
        min_index = 0
        
        do i = 1, Lu
            next_index = i * (i + 2)
            
            norm_factor = sqrt(real(i * (i + 1), kind=8) / real(2 * i + 1, kind=8))
            do concurrent (j = min_index+1: next_index)
                norm_matrix(j, j) = norm_factor
                norm_matrix(j + Nu, j + Nu) = norm_factor
                inverse_norm_matrix(j, j) = 1.0d0 / norm_factor
                inverse_norm_matrix(j + Nu, j + Nu) = 1.0d0 / norm_factor
            end do
            
            min_index = next_index
            
        end do

    end subroutine
!========================================================================================================================== 
    
end module