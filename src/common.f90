module common
    use utilities
    use deterministic_rng, only: normal_vector, RNG_FORECAST
    use linear_interpolation_module
    use blas95
    use lapack95
    use f95_precision
    use legendre
    use omp_lib
    use mpi
    implicit none
    
    
contains
!========================================================================================================================== 
    subroutine sample_timed_quantity(times, X_quantity, sampling_dt, result_1, result_2)
    !*****************************************************************************************************************
    !"""
    !Samples a timed quantity according to sampling_dt.
    !First, generates sampling times from the times array and then evaluates the X quantity at these times using interpolation.
    !
    !:param times: times of X quantity
    !:type times: 1D numpy.array (dim: Ntimes)
    !:param X_quantity: array containing values of the quantity X at times
    !:type X_quantity: 2D numpy.array (dim: Ntimes x Ncoefs)
    !:param sampling_dt: time step to use for the sampling (in years)
    !:type sampling_dt: float
    !:return: 2-tuple containing the sampling times and the quantity X evaluated at these sampling times
    !:rtype: 1D numpy.array (dim: Nsamples), 2D numpy.array (dim: Nsamples x Ncoefs)
    !"""
    !*****************************************************************************************************************
        real(kind=8), intent(in) :: times(:)
        real(kind=8), intent(in) :: X_quantity(:,:)
        real(kind=8), intent(in) :: sampling_dt
        real(kind=8), allocatable, intent(out) :: result_1(:), result_2(:,:)
        integer :: nb_coeffs, nb_samples, i, j, iflag
        type(linear_interp_1d) :: s1        
        real(kind=8) :: fval
        
        nb_coeffs = SIZE(X_quantity, 2)
        nb_samples = CEILING((MAXVAL(times) - MINVAL(times)) / sampling_dt)
        allocate(result_1(nb_samples))
        allocate(result_2(nb_samples, nb_coeffs))
        result_1 = 0.0d0
        result_2 = 0.0d0
        i = 0
        do while (MINVAL(times) + i*sampling_dt < MAXVAL(times))
            result_1(i+1) = MINVAL(times) + i*sampling_dt
            i = i + 1
        end do
        
        do i = 1, nb_coeffs
            call s1.initialize(times,X_quantity(:,i),iflag)
            
            do j = 1, nb_samples
                call s1.evaluate(result_1(j),fval)
                result_2(j,i) = fval
            end do
            
            call s1.destroy()
        end do       

    end subroutine    
!==========================================================================================================================
    
!========================================================================================================================== 
    subroutine compute_legendre_polys(tmax, Lb, Lu, Lsv, legendre_polys)
    !*****************************************************************************************************************
    !"""
    !Computes the coefficients of Legendre polynomials, first derivative and second derivative of B, U and SV using a wrapped Fortran function.
    !
    !:param tmax: Number of angles for computation
    !:type tmax: int
    !:param Lb: Max degree of magnetic field
    !:type Lb: int
    !:param Lu: Max degree of core flow
    !:type Lu: int
    !:param Lsv: Max degree of secular variation
    !:type Lsv: int
    !:return: A dictionary containing the angles used for the computations and the Legendre coefs of B, U and SV.
    !:rtype: dict
    !"""
    !*****************************************************************************************************************
        integer, intent(in) :: tmax, Lb, Lu, Lsv
        class(legendre_polys_type), intent(out) :: legendre_polys
        integer :: LLb, LLu, LLsv, i
        real(kind=8) :: x_i
        real(kind=8), allocatable :: lp_b(:,:), d_lp_b(:,:), d2_lp_b(:,:)
        real(kind=8), allocatable :: lp_u(:,:), d_lp_u(:,:), d2_lp_u(:,:)
        real(kind=8), allocatable :: lp_sv(:,:), d_lp_sv(:,:), d2_lp_sv(:,:)
        real(kind=8), allocatable :: gauss_points(:), gauss_weights(:), gauss_thetas(:)
        
        !# Compute number of coefs
        LLb = ((Lb + 1) * (Lb + 2)) / 2
        LLu = ((Lu + 1) * (Lu + 2)) / 2
        LLsv = ((Lsv + 1) * (Lsv + 2)) / 2
        
        !# Init of arrays storing the legendre polynomials, their derivative and second derivative
        !# For magnetic field
        allocate(lp_b(LLb, tmax), source=0.0d0)
        allocate(d_lp_b(LLb, tmax), source=0.0d0)
        allocate(d2_lp_b(LLb, tmax), source=0.0d0)
        
        !# For core flow
        allocate(lp_u(LLu, tmax), source=0.0d0)
        allocate(d_lp_u(LLu, tmax), source=0.0d0)
        allocate(d2_lp_u(LLu, tmax), source=0.0d0)
        
        !# For secular variation
        allocate(lp_sv(LLsv, tmax), source=0.0d0)
        allocate(d_lp_sv(LLsv, tmax), source=0.0d0)
        allocate(d2_lp_sv(LLsv, tmax), source=0.0d0)
        
        !# Compute Gauss-Legendre quadrature : returns tmax gauss_points in [-1, -1] with their associated weights
        call gauleg(-1.0d0, 1.0d0, gauss_points, gauss_weights, tmax)
        
        !# For each gauss point x_i, compute legendre associated functions at x_i
        do i = 1, SIZE(gauss_points)
            x_i = gauss_points(i)
            !# Last arg of plmbar2 is 1 to have Schmidt quasi-normalisation
            call plmbar2(lp_b(:,i), d_lp_b(:,i), d2_lp_b(:,i), x_i, Lb, 1)
            call plmbar2(lp_u(:,i), d_lp_u(:,i), d2_lp_u(:,i), x_i, Lu, 1)
            call plmbar2(lp_sv(:,i), d_lp_sv(:,i), d2_lp_sv(:,i), x_i, Lsv, 1)
        end do
        
        !# Convert the values used for the polynomial computation in angles
        allocate(gauss_thetas, source=gauss_points)
        gauss_thetas = acos(gauss_points)
        
        allocate(legendre_polys.thetas, source=gauss_thetas)
        allocate(legendre_polys.weights, source=gauss_weights)
        
        allocate(legendre_polys.MF(3, LLb, tmax))
        legendre_polys.MF(1, :, :) = lp_b
        legendre_polys.MF(2, :, :) = d_lp_b
        legendre_polys.MF(3, :, :) = d2_lp_b
        
        allocate(legendre_polys.U(3, LLu, tmax))
        legendre_polys.U(1, :, :) = lp_u
        legendre_polys.U(2, :, :) = d_lp_u
        legendre_polys.U(3, :, :) = d2_lp_u
        
        allocate(legendre_polys.SV(3, LLsv, tmax))
        legendre_polys.SV(1, :, :) = lp_sv
        legendre_polys.SV(2, :, :) = d_lp_sv
        legendre_polys.SV(3, :, :) = d2_lp_sv
    end subroutine    
!==========================================================================================================================
    
!========================================================================================================================== 
    subroutine compute_Kalman_gain_matrix(P, H, R, use_cholesky, Kalman_gain)
    !*****************************************************************************************************************
    !"""
    !Computes the Kalman gain matrix from the correlation matrix, observation operator and observation error matrix:
    !    K = P*trans(H)*inv(H*P*trans(H)+R)
    !The inversion can be carried out using Cholesky decomposition.
    !
    !:param P: Correlation matrix
    !:type P: 2D numpy.ndarray (dim: Nstate x Nstate)
    !:param H: Observation operator
    !:type H: 2D numpy.ndarray (dim: Nobs x Nstate)
    !:param R: Observation error matrix
    !:type R: 2D numpy.ndarray (dim: Nobs x Nobs)
    !:param use_cholesky: if True, matrix inversion is done using Cholesky decomposition. Uses scipy.linalg.inv otherwise.
    !:type use_cholesky: bool
    !:return: Kalman gain matrix K = P*trans(H)*inv(H*P*trans(H)+R)
    !:rtype: 2D numpy.ndarray (dim: Nstate x Nobs)
    !"""
    !*****************************************************************************************************************
        real(kind=8), intent(in) :: P(:,:), H(:,:), R(:,:)
        logical, intent(in) :: use_cholesky
        real(kind=8), allocatable, intent(out) :: Kalman_gain(:,:)
        real(kind=8), allocatable :: PHT(:,:), matrix_to_invert(:,:), inverted_matrix(:,:)
        real(kind=8), allocatable :: S(:,:), KT(:,:), ST(:,:)
        integer :: info
        
        
        !# Compute PHT = P*transpose(H) that is needed twice
        allocate(PHT, source=MATMUL(P, TRANSPOSE(H)))
        
        !# Compute the Kalman gain K = P*trans(H)*inv(H*P*trans(H)+R)
        allocate(matrix_to_invert, source=MATMUL(H, PHT) + R)
        call potrf(matrix_to_invert, 'L', info) 
        !!! atention: potrf modifies the input matrix, so we need to make a copy of it if we want to keep it for later use
        
        if (use_cholesky .AND. (info == 0)) then
            allocate(KT, source=transpose(PHT))
            call potrs(matrix_to_invert, KT, 'L',info)
            allocate(Kalman_gain, source=transpose(KT))
        else
            allocate(KT, source=transpose(PHT))
            allocate(ST, source=transpose(MATMUL(H, PHT) + R))
            call gesv(ST, KT, info=info)
            allocate(Kalman_gain, source=transpose(KT))
        end if
        
    end subroutine    
!==========================================================================================================================
    
!========================================================================================================================== 
    subroutine compute_Kalman_huber(b_f, y, P_eig_val, P_eig_vec, H, Rbb, max_steps, b_kplus1)
    !*****************************************************************************************************************
    !"""
    !Solve iteratively the least square problem with a Huber norm, the equation is written (cf Walker and jackson 2000):
    !b^k+1 = b^f + (P^{-1/2} + H^T R^{-1/2} W R^{-1/2} H)^-1 ( H^T R^{-1/2} W R^{-1/2}) (y - H b^f)
    !
    !:param max_steps: maximum number of realisations
    !:type max_steps: int 
    !"""
    !*****************************************************************************************************************
        real(kind=8), intent(in) :: b_f(:), y(:), P_eig_val(:), P_eig_vec(:,:), H(:,:), Rbb(:,:)
        integer, intent(in) :: max_steps
        real(kind=8), intent(inout) :: b_kplus1(:)
        real(kind=8), allocatable :: P_diag(:,:), Pbb_inv(:,:)
        real(kind=8), allocatable :: sq_R(:), R_H(:,:), H_R(:,:), W(:), ones_min(:), y_minus_H_b(:)
        real(kind=8), allocatable :: A(:,:), B(:), diag(:,:), epsilon(:), b_k(:)
        integer :: info, i
        real(kind=8) :: c, cutoff
        
        ! # calculate P_diag   (1/lamda)
        allocate(P_diag, source=P_eig_vec)
        P_diag = 0.0d0
        do i = 1, SIZE(P_diag, 1)
            P_diag(i,i) = 1.0d0 / P_eig_val(i)
        end do
        Pbb_inv = MATMUL(P_eig_vec, MATMUL(P_diag, TRANSPOSE(P_eig_vec)))
        
        ! # calculate Rbb^{-1/2}
        allocate(sq_R(SIZE(Rbb, 1)))
        sq_R = 0.0d0
        do i = 1, SIZE(Rbb, 1)
            sq_R(i) = 1.0d0 / SQRT(Rbb(i,i))
        end do
        
        ! Calculate R_H = Rbb^{-1/2} * H.
        ! Since Rbb^{-1/2} is diagonal, this is implemented by scaling each row of H
        ! instead of using MATMUL.
        allocate(R_H(SIZE(H, 1), SIZE(H, 2)))
        R_H = 0.0d0
        do i = 1, SIZE(H, 1)
            R_H(i,:) = sq_R(i) * H(i,:)
        end do
        allocate(H_R, source=TRANSPOSE(R_H))
        
        ! Initialize Huber weights to 1 (equal weights for all observations)
        allocate(W(SIZE(Rbb,1)))
        W = 1.0d0
        allocate(ones_min(SIZE(Rbb,1)))
        ones_min = 1.0d0
        
        ! parameter and y_minus_H_b
        c = 1.5
        cutoff = 1.0d-8
        allocate(y_minus_H_b, source=y-MATMUL(H, b_f))
        
        ! calcualte b_kplus1
        allocate(A, source=Pbb_inv)
        A = 0.0d0
        allocate(B, source=y_minus_H_b)
        B = 0.0d0
        allocate(diag(SIZE(W,1), SIZE(W,1)))
        diag = 0.0d0
        b_kplus1 = 0.0d0
        allocate(b_k, source=b_f)
        b_k = 0.0d0
        allocate(epsilon, source=y)
        do i = 1, max_steps
            call diag_(W, diag)
            A = Pbb_inv + MATMUL(H_R, MATMUL(diag, R_H))
            call diag_(W * sq_R, diag)
            B = MATMUL(H_R, MATMUL(diag, y_minus_H_b))
            call gesv(A, B, info=info)
            if (info /= 0) then
                print *, "Error in solving the linear system, info = ", info
                stop
            end if
            b_kplus1 = b_f + B
            ! # if difference between last iteration and this one is less than 2% of the deviation, break the loop
            if (i>1) then
                if (maxval(abs(b_kplus1 - b_k)) < 1.0d-4 * maxval(abs(b_k))) exit
            end if
            b_k = b_kplus1
            epsilon = abs(y - matmul(H, b_kplus1)) * sq_R
            epsilon = max(epsilon, cutoff)
            !# if residual is bigger than one, weight is replaced by one
            W = MIN(c / epsilon, ones_min)
            !# conditions to break the loop to make it faster
            !# then nothing will change in the loop, LS solution because no outliers
            if (all(abs(W - 1.0d0) < 1.0d-12)) exit
        end do        
    end subroutine    
!==========================================================================================================================
    
!========================================================================================================================== 
    subroutine compute_Kalman_huber_parameter_basis(b_f, y, HPzzHT, PzzHT, H, Rbb, norm, n_steps, b_kplus1)
    !*****************************************************************************************************************
    !"""
    !very similar to compute_Kalman_huber but written in the parameter basis.
    !Solve iteratively the least square problem with a Huber norm, the equation is written (cf Walker and jackson 2000):
    !b^k+1 = b^f + (P^{-1/2} + H^T R^{-1/2} W R^{-1/2} H)^-1 (H^T R^{-1/2} W R^{-1/2}) (y - H b^f)
    !
    !:param max_steps: maximum number of realisations
    !:type max_steps: int (dim: N_real x N_coefs_A)
    !"""
    !*****************************************************************************************************************
        real(kind=8), intent(in) :: b_f(:), y(:), HPzzHT(:,:), PzzHT(:,:), H(:,:), Rbb(:,:)
        integer, intent(in) :: n_steps
        character(len=*), intent(in) :: norm
        real(kind=8), intent(inout) :: b_kplus1(:)
        real(kind=8), allocatable :: sq_R(:), W(:), ones_min(:), y_minus_H_b(:)
        real(kind=8), allocatable :: R_k(:,:), A(:,:), B(:), epsilon(:), prev_W(:)
        integer :: info, i, j, max_steps
        real(kind=8) :: c, cutoff
        
        ! # Calculate sqrt(diag(Rbb))
        allocate(sq_R(SIZE(Rbb, 1)))
        sq_R = 0.0d0
        do i = 1, SIZE(Rbb, 1)
            sq_R(i) = SQRT(Rbb(i,i))
        end do
        
        ! Initialize Huber weights to 1 (equal weights for all observations)
        allocate(W(SIZE(Rbb,1)))
        W = 1.0d0
        allocate(ones_min(SIZE(Rbb,1)))
        ones_min = 1.0d0
        
        ! parameter and y_minus_H_b
        c = 1.5
        cutoff = 1.0d-8
        allocate(y_minus_H_b, source=y-MATMUL(H, b_f))
        
        max_steps = n_steps
        if (trim(norm) == 'L2' .OR. trim(norm) == 'l2') then
            max_steps = 1 !# first step just computes least square, one iteration in the loop is fine
        end if
        
        !# initialize matrix
        allocate(R_k(SIZE(sq_R), SIZE(sq_R)))
        allocate(A, source = HPzzHT)
        allocate(B, source=y_minus_H_b)
        allocate(epsilon, source=y)
        allocate(prev_W, source=W)
        
        do i = 1, max_steps
            R_k = 0.0d0
            do j = 1, SIZE(R_k, 1)
                R_k(j, j) = sq_R(j) * sq_R(j) / W(j)
            end do
            A = HPzzHT + R_k
            B = y_minus_H_b
            call gesv(A, B, info=info)
            
            if (info /= 0) then
                print *, "Error in solving the linear system, info = ", info
                stop
            end if
            b_kplus1 = b_f + MATMUL(PzzHT, B)
            
            !# if difference between last iteration and this one is less than 2% of the deviation, break the loop
            epsilon = abs(y - MATMUL(H, b_kplus1)) / sq_R
            epsilon = max(epsilon, cutoff)
            
            prev_W = W
            W = MIN(c / epsilon, ones_min) !# if residual is bigger than c, weight is replaced by one
            !# conditions to break the loop to make it faster
            if (all(W == ones_min)) exit
            
            ! Python-compatible stopping criterion for comparison
            if (i > 1 .and. maxval(abs(W - prev_W)) < 1.0e-4 * maxval(W)) exit
            !if (i > 1 .and. maxval(W - prev_W) < 1.0e-4 * maxval(W)) exit
        end do
    end subroutine    
!==========================================================================================================================
    
!========================================================================================================================== 
    subroutine compute_misfit(realisations, mean, inverse_weight_matrix, result_)
    !*****************************************************************************************************************
    !"""
    !Computes the weighted observation-space misfit (Y - HX) using the Mahalanobis distance.
    !The misfit is a weighted norm according to the weight_matrice (Mahalanobis distance).
    !:warning: The inverse weight matrix must be supplied
    !
    !:param realisations: Array of realisations (dim: nb_real x nb_coefs)
    !:type realisations: np.array
    !:param mean: Mean that will be subtracted to the realisations array (dim: nb_coefs)
    !:type mean: np.array
    !:param inverse_weight_matrix: Inverse of the weight matrix (WM^{-1}, dim: nb_coefs x nb_coefs)
    !:type inverse_weight_matrix: np.array
    !:return: sqrt(1/(nb_real(nb_coefs-1))*sum[(realisations - mean)^T WM^{-1} (realisations-mean)])
    !:rtype: float
    !"""
    !*****************************************************************************************************************
        real(kind=8), intent(in) :: realisations(:,:), mean(:,:), inverse_weight_matrix(:,:)
        real(kind=8), intent(inout) :: result_
        real(kind=8), allocatable :: distance(:,:)
        integer :: nb_realisations, nb_coefs, i
        real(kind=8) :: non_norm_squared_norm
        
        
        nb_realisations = size(realisations, 1)
        nb_coefs = size(realisations, 2)
        
        !# Regular distance D (this step will fail if the dimensions are not matching)
        allocate(distance, source=realisations - mean)
        
        !# Compute the sum over realisations (assumed axis=0) of D^T * WM^{-1} * D that is the squared weighted norm.
        non_norm_squared_norm = 0.0d0
        do i = 1, nb_realisations
            non_norm_squared_norm = non_norm_squared_norm + DOT_PRODUCT(distance(i,:), MATMUL(inverse_weight_matrix, distance(i,:)))
        end do
        
        result_ = SQRT(non_norm_squared_norm / (nb_realisations * (nb_coefs - 1)))
    end subroutine    
!==========================================================================================================================
    
!========================================================================================================================== 
    subroutine get_BLUE(X, Y, P, H, R, K, cholesky, BLUE_X)
    !*****************************************************************************************************************
    !"""
    !Returns the Best Linear Unbiased Estimate (BLUE) for X from the observations Y and the covariances matrices.
    !
    !:param X: current state vector
    !:type X: 1D numpy.ndarray (dim: Nstate)
    !:param Y: observation vector
    !:type Y: 1D numpy.ndarray (dim: Nobs)
    !:param P: covariance of forecast priors
    !:type P: 2D numpy.ndarray (dim: Nstate x Nstate)
    !:param H: Observation operator
    !:type H: 2D numpy.ndarray (dim: Nobs x Nstate)
    !:param R: covariance of observations errors
    !:type R: 2D numpy.ndarray (dim: Nobs x Nobs)
    !:param K: Kalman gain matrix. Optionnal: if not supplied, will be computed from other matrices.
    !:type K: 2D numpy.ndarray (dim: Nstate x Nobs)
    !:param cholesky: if True, uses Cholesky decomposition to compute matrices' inverses
    !:type cholesky: bool
    !:return: Best Linear Unbiased Estimate of X
    !:rtype: 1D numpy.ndarray (dim: Nstate)
    !"""
    !*****************************************************************************************************************
        real(kind=8), intent(in) :: X(:), Y(:), P(:,:), H(:,:), R(:,:), K(:,:)
        logical, intent(in) :: cholesky
        real(kind=8), intent(inout) :: BLUE_X(:)
        real(kind=8), allocatable :: D(:), matrix_to_invert(:,:), inverted_matrix(:,:)
        real(kind=8), allocatable :: S(:,:), KT(:,:), ST(:,:)
        integer :: info
        
        !# Compute the innovation vector first
        allocate(D, source=Y - MATMUL(H, X))
        
        !# Compute Kalman gain matrix if not provided
        ! not needed here since K is provided as input
        
        BLUE_X = X + MATMUL(K, D)     
    end subroutine    
!==========================================================================================================================
    
!==========================================================================================================================    
    subroutine compute_diag_AR1_coefs(cov_U, cov_ER, Tau_U, Tau_E, A, Chol)
    !*****************************************************************************************************************
    !"""
    !Computes the matrices for a diagonal AR-1 process
    !
    !:param correlation_matrix: Correlation matrix of the AR-1 quantity
    !:type correlation_matrix: np.array (dim: N x N)
    !:param ar_timestep: Timestep of the AR-1 process in years (also called Tau)
    !:type ar_timestep: float
    !:return: The drift matrix and Cholesky lower matrix of the correlation matrix
    !:rtype: np.array (dim: N x N), np.array (dim: N x N)
    !"""
    !*****************************************************************************************************************
        real(kind=8), intent(in) :: cov_U(:,:), cov_ER(:,:)
        REAL(kind=8), intent(in) :: Tau_U, Tau_E
        real(kind=8), allocatable, intent(out) :: A(:,:), Chol(:,:)
        integer :: Nu, Ner, Nz, i, j
        real(kind=8), allocatable :: matrix_tmp1(:,:), matrix_tmp2(:,:)
        real(kind=8), allocatable :: Chol_U(:,:), Chol_ER(:,:)
        
        Nu = SIZE(cov_U, 1)
        Ner = SIZE(cov_ER, 1)
        Nz = Nu + Ner
        allocate(A(Nz, Nz), source=0.0d0)
        do concurrent (i = 1:Nz)
            if (i <= Nu) then
                A(i, i) = 1.0d0 / Tau_U
            else
                A(i, i) = 1.0d0 / Tau_E
            end if
        end do
        
        allocate(matrix_tmp1, source=cov_U)
        allocate(matrix_tmp2, source=cov_ER)
        call potrf(matrix_tmp1, 'L')
        do j=1, SIZE(matrix_tmp1, 1)
            do concurrent (i=1:j-1)
                matrix_tmp1(i,j) = 0.0d0
            end do
        end do
        
        call potrf(matrix_tmp2, 'L')
        do j=1, SIZE(matrix_tmp2, 1)
            do concurrent (i=1:j-1)
                matrix_tmp2(i,j) = 0.0d0
            end do
        end do
        
        allocate(Chol_U, source=matrix_tmp1)
        allocate(Chol_ER, source=matrix_tmp2)
        Chol_U = SQRT(2.0d0 / Tau_U) * Chol_U
        Chol_ER = SQRT(2.0d0 / Tau_E) * Chol_ER
        call block_diag(Chol_U, Chol_ER, Chol)
    end subroutine
!==========================================================================================================================
    
!==========================================================================================================================
    subroutine prep_AR_matrix(Z, dt_samp, Nt, x0_avg, dX_avg, d2X_avg, d3X_avg)
    !*****************************************************************************************************************
    !"""
    !Compute time derivatives
    !
    !:Z param: state array
    !:Z type: ndarray (Ntimes x Ncoeffs)
    !:dt_samp param: time step
    !:dt_samp type: float
    !:Nt param: average window size
    !:Nt type: int
    !:X return: X
    !:X rtype: ndarray (Ntimes x Ncoeffs)
    !:dX return: first time derivative X
    !:dX rtype: ndarray (Ntimes x Ncoeffs)
    !:d2X return: second time derivative X
    !:d2X rtype: ndarray (Ntimes x Ncoeffs)
    !:d3X return: third time derivative X
    !:d3X rtype: ndarray (Ntimes x Ncoeffs)
    !"""
    !*****************************************************************************************************************
        real(kind=8), intent(in) :: Z(:,:)
        real(kind=8), intent(in) :: dt_samp
        integer, intent(in) :: Nt
        real(kind=8), allocatable :: x0(:,:), dX(:,:), d2X(:,:), d3X(:,:)
        real(kind=8), allocatable :: x1(:,:), x2(:,:), x3(:,:)
        real(kind=8), allocatable, intent(out) ::x0_avg(:,:), dX_avg(:,:), d2X_avg(:,:), d3X_avg(:,:)
        
        if (SIZE(Z, 1) < 4) then
            write (10,'(A)') "First dimension of Z must be at least 4."
        end if
        
        allocate(x0(SIZE(Z,1)-3,  SIZE(Z,2)))
        allocate(dX(SIZE(Z,1)-3,  SIZE(Z,2)), x1(SIZE(Z,1)-3,  SIZE(Z,2)))
        allocate(d2X(SIZE(Z,1)-3, SIZE(Z,2)), x2(SIZE(Z,1)-3,  SIZE(Z,2)))
        allocate(d3X(SIZE(Z,1)-3, SIZE(Z,2)), x3(SIZE(Z,1)-3,  SIZE(Z,2)))
        
        x0 = Z(1:SIZE(Z,1)-3,:)
        x1 = Z(2:SIZE(Z,1)-2,:)
        x2 = Z(3:SIZE(Z,1)-1,:)
        x3 = Z(4:SIZE(Z,1),:)
        !# compute dX d2X d3X derivatives
        dX  = (x1 - x0) / dt_samp
        d2X = (x2 - 2.0d0*x1 + x0) / (dt_samp**2)
        d3X = (x3 - 3.0d0*x2 + 3.0d0*x1 - x0) / (dt_samp**3)
        
        !# compute average (not divided by smoothing window total weigth !!!)
        !# requires to be divided by the smoothing window total weigth to get a proper average
        call compute_average(x0, Nt, 'valid', x0_avg)
        call compute_average(dX, Nt, 'valid', dX_avg)
        call compute_average(d2X, Nt, 'valid', d2X_avg)
        call compute_average(d3X, Nt, 'valid', d3X_avg)
        
    end subroutine
!==========================================================================================================================
    
!==========================================================================================================================
    subroutine compute_AR_coefs_avg(container, AR_type, A, B, C,Chol)
    !*****************************************************************************************************************
    !"""
    !Compute the dense AR coefs:
    !A,Chol if AR1
    !A,B,C,Chol if AR3
    !
    !:param container: list containing the stored variables [X, dX, dt_sampling] for every prior.
    !:type container: list
    !:param AR_type: type of AR process
    !:type AR_type: str
    !"""
    !*****************************************************************************************************************
        class(container_type), intent(in) :: container(:)
        character(len=*), intent(in) :: AR_type
        real(kind=8), allocatable, intent(out) :: A(:,:), B(:,:), C(:,:), Chol(:,:)
        integer :: i, NZ, NT_, Nt, j
        real(kind=8), allocatable :: tmp1(:,:), tmp2(:,:), X(:,:), dX(:,:), XXT(:,:), dXXT(:,:), tmp2_inv(:,:)
        real(kind=8), allocatable :: AA(:,:)
        real(kind=8) :: Dt, Mt
        real(kind=8), allocatable :: black_(:), W(:,:), S(:,:)
        class(WWT), allocatable :: WWT(:)
        
        if (trim(AR_type) == "AR1") then
            NZ = SIZE(container(1).z_T, 1)
            allocate(tmp1(NZ, NZ), source=0.0d0)
            allocate(tmp2(NZ, NZ), source=0.0d0)
        else if (trim(AR_type) == "AR3") then
            NZ = SIZE(container(1).z_T, 1) / 3
            allocate(tmp1(NZ*3, NZ*3), source=0.0d0)
            allocate(tmp2(NZ*3, NZ*3), source=0.0d0)
        end if
        
        Dt = container(1).dt_prior
        
        ! # compute AA
        NT_ = 0
        do i = 1, SIZE(container)
            if (ALLOCATED(X)) deallocate(X)
            if (ALLOCATED(dX)) deallocate(dX)
            if (ALLOCATED(XXT)) deallocate(XXT)
            if (ALLOCATED(dXXT)) deallocate(dXXT)
            
            allocate(X, source=container(i).z_T)
            allocate(dX, source=container(i).dz_T)
            allocate(XXT(SIZE(X, 1), SIZE(X, 1)))
            allocate(dXXT(SIZE(dX, 1), SIZE(X, 1)))
            Nt = container(i).Nt
            call blackman(Nt, black_) 
            Mt = dot(black_, black_)
            NT_ = NT_ + NINT(real(SIZE(X, 2))/Nt)
            XXT = MATMUL(X, TRANSPOSE(X)) / (Mt*Nt)
            dXXT = MATMUL(dX, TRANSPOSE(X)) / (Mt*Nt)
            tmp1 = tmp1 + dXXT
            tmp2 = tmp2 + XXT            
        end do
        
        allocate(AA, source=tmp1)
        call max_inv(tmp2, tmp2_inv)
        AA = MATMUL(tmp1, tmp2_inv) * (-1.0d0 / Dt)
        
        !# compute WWT
        allocate(WWT(SIZE(container)))
        do i = 1, SIZE(container)
            if (ALLOCATED(X)) deallocate(X)
            if (ALLOCATED(dX)) deallocate(dX)
            if (ALLOCATED(W)) deallocate(W)
            
            allocate(X, source=container(i).z_T)
            allocate(dX, source=container(i).dz_T)
            allocate(W, source=dX)
            Dt = container(i).dt_prior
            Nt = container(i).Nt
            call blackman(Nt, black_) 
            Mt = dot(black_, black_)
            W = (dX + Dt * MATMUL(AA, X)) / SQRT(Dt)
            allocate(WWT(i).matrix(SIZE(W, 1), SIZE(W, 1)))
            WWT(i).matrix = MATMUL(W, TRANSPOSE(W)) / (Mt*Nt)
        end do
        
        allocate(S(SIZE(W, 1), SIZE(W, 1)), source=0.0d0)
        if (trim(AR_type) == "AR1") then
            do i = 1, SIZE(container)
                S = S + WWT(i).matrix
            end do
            S = S / NT_
            allocate(Chol, source=S)
            allocate(A, source=AA)
            call potrf(Chol, 'L')
            do j=1, SIZE(Chol, 1)
                do concurrent (i=1:j-1)
                    Chol(i,j) = 0.0d0
                end do
            end do
            A = TRANSPOSE(AA)
        else if (trim(AR_type) == "AR3") then
            do i = 1, SIZE(container)
                S = S + WWT(i).matrix
            end do
            S = S / NT_
            allocate(Chol(NZ, NZ))
            allocate(A(NZ, NZ), B(NZ, NZ), C(NZ, NZ))
            Chol = S(2*NZ+1:3*NZ, 2*NZ+1:3*NZ)
            call potrf(Chol, 'L')
            do j=1, SIZE(Chol, 1)
                do concurrent (i=1:j-1)
                    Chol(i,j) = 0.0d0
                end do
            end do
            C = AA(2*NZ+1:3*NZ, 1:NZ)
            B = AA(2*NZ+1:3*NZ, NZ+1:2*NZ)
            A = AA(2*NZ+1:3*NZ, 2*NZ+1:3*NZ)
            C = TRANSPOSE(C)
            B = TRANSPOSE(B)
            A = TRANSPOSE(A)
        end if
    end subroutine
!==========================================================================================================================
    
!==========================================================================================================================
    subroutine compute_AR1_coefs_forecast(A, Chol, dt_forecast, Ncoef, A_forecast, Chol_forecast)
    !*****************************************************************************************************************
    !"""
    !Compute the forecast coeffcients A, and Chol for AR-1 process
    !"""
    !*****************************************************************************************************************
        real(kind=8), intent(in) :: A(:,:), Chol(:,:)
        real(kind=8), intent(in) :: dt_forecast
        integer, intent(in) :: Ncoef
        real(kind=8), allocatable :: Id(:,:)
        real(kind=8), allocatable, intent(out) :: A_forecast(:,:), Chol_forecast(:,:)
        
        !# matrices pour AR1 sch¨¦ma d¨¦centr¨¦
        allocate(Id(Ncoef, Ncoef))
        call dlaset('A', Ncoef, Ncoef, 0.0d0, 1.0d0, Id, Ncoef)
        allocate(A_forecast, source=A)
        allocate(Chol_forecast, source=Chol)
        
        A_forecast = dt_forecast * A - Id
        Chol_forecast = SQRT(dt_forecast) * Chol
    end subroutine
!==========================================================================================================================
    
!==========================================================================================================================
    subroutine compute_AR3_coefs_forecast(A, B, C, Chol, dt_forecast, Ncoef, A_forecast, B_forecast, C_forecast, Chol_forecast)
    !*****************************************************************************************************************
    !"""
    !Compute the forecast coeffcients A, B, C and Chol for AR-3 process
    !"""
    !*****************************************************************************************************************
        real(kind=8), intent(in) :: A(:,:), B(:,:), C(:,:), Chol(:,:)
        real(kind=8), intent(in) :: dt_forecast
        integer, intent(in) :: Ncoef
        real(kind=8), allocatable :: Id(:,:)
        real(kind=8), allocatable, intent(out) :: A_forecast(:,:), B_forecast(:,:), C_forecast(:,:) ,Chol_forecast(:,:)
        
        !# matrices pour AR3 sch¨¦ma d¨¦centr¨¦
        allocate(Id(Ncoef, Ncoef))
        call dlaset('A', Ncoef, Ncoef, 0.0d0, 1.0d0, Id, Ncoef)
        allocate(A_forecast, source=A)
        allocate(B_forecast, source=B)
        allocate(C_forecast, source=C)
        allocate(Chol_forecast, source=Chol)
        
        A_forecast = dt_forecast * A - 3.0d0 * Id
        B_forecast = (dt_forecast**2) * B - 2.0d0 * dt_forecast * A + 3.0d0 * Id
        C_forecast = (dt_forecast**3) * C - (dt_forecast**2) * B + dt_forecast * A - Id
        Chol_forecast = (dt_forecast**(5.0d0/2.0d0)) * Chol
    end subroutine
!==========================================================================================================================
    
!==========================================================================================================================
    subroutine compute_average(X, Nt, mode, result_)
    !# Perform a convolution of X by a blackman smoothing window of size Nt (NOT A PROPER AVERAGE)    
        real(kind=8), intent(in) :: X(:,:)
        integer, intent(in) :: Nt
        character(len=*), intent(in) :: mode
        real(kind=8), allocatable, intent(out) :: result_(:,:)
        real(kind=8), allocatable :: X_avg(:,:)
        integer :: i, j
        real(kind=8), allocatable :: blackman_w(:)
        !integer :: comm, rank, nb_proc, ierr
        
        allocate(X_avg(SIZE(X, 1)-Nt+1, SIZE(X, 2)))
        allocate(blackman_w(Nt))
        allocate(result_(SIZE(X, 1)-Nt+1, SIZE(X, 2)))
        X_avg = 0.0d0
        call blackman(Nt, blackman_w)
        
        !comm = MPI_COMM_WORLD
        !call MPI_Comm_size(comm, nb_proc, ierr)
        !call MPI_Comm_rank(comm, rank, ierr)       
        
        if (TRIM(mode) == 'valid') then
            !$omp parallel do collapse(2)
            do i = 1, SIZE(X, 2)
                do j = 1, SIZE(X_avg, 1)
                    X_avg(j, i) = DOT_PRODUCT(X(j:j+Nt-1,i), blackman_w)
                    !if (rank==1) print *, "rank=", rank, "threads=", omp_get_max_threads()
                end do
            end do
            !$omp end parallel do
        end if
        result_ = X_avg        
    end subroutine
!==========================================================================================================================
    
!==========================================================================================================================
    subroutine ar1_process(X, A, Chol, global_seed, global_i_real, i_t, check_Cholesky, AR1_result)
        !"""
        !Applies an Auto-Regressive process of order 1 to the augmented state Z.
        !
        !X1 = - X0 @ A + Chol @ normal_noise
        !
        !:param X: quantity on which the AR-1 process will be applied (X0 above)
        !:type X: 1D numpy.ndarray (dim: Nz)
        !:param A: AR-1 operator (can be diagonal or dense)
        !:type A: 2D numpy.ndarray (dim: Nz x Nz)
        !:param global_seed: global random seed
        !:type global_seed: int
        !:param global_i_real: global model realisation index
        !:type global_i_real: int
        !:param i_t: time index
        !:type i_t: int
        !:param check_Cholesky: If True (default), checks that the Cholesky is lower triangular. Setting to False may enhance performance.
        !:type: bool
        !:return: vector containing the result of the AR-1 process (X1 above)
        !:rtype: 1D numpy.ndarray (dim: Nz)
        !"""
        real(kind=8), intent(in) :: X(:)
        real(kind=8), intent(in) :: A(:,:), Chol(:,:)
        real(kind=8), allocatable, intent(out) :: AR1_result(:)
        integer, intent(in) :: global_seed, global_i_real, i_t
        logical, intent(in) :: check_Cholesky
        integer :: Ncoef, i, j
        real(kind=8), allocatable :: normal_noise(:), Chol_test(:,:)
        
        Ncoef = SIZE(X, 1)
        if ((SIZE(A, 1) .ne. Ncoef) .or. (SIZE(A, 2) .ne. Ncoef)) then
            write (10,'(A, i4, i4, A, i4, i4, A)')  'A matrix in AR-1 should have dimensions (',Ncoef,Ncoef,'). Got', SIZE(A, 1),SIZE(A, 1),' instead.'
            write (*,'(A, i4, i4, A, i4, i4, A)')  'A matrix in AR-1 should have dimensions (',Ncoef,Ncoef,'). Got', SIZE(A, 1),SIZE(A, 1),' instead.'
            stop
        end if
        allocate(Chol_test, source=Chol)
        do i = 1, SIZE(Chol_test, 1)
            do j = i+1, SIZE(Chol_test, 2)
                Chol_test(i,j) = 0.0d0
            end do
        end do
        if (check_Cholesky)then
            if (ANY(ABS(Chol_test-Chol) > 1.0d-12)) then
                write (10,'(A)')  'Cholesky matrix supplied in AR-1 is not a lower triangular matrix ! Did you supply the upper Cholesky matrix instead ?'
                write (*,'(A)') 'Cholesky matrix supplied in AR-1 is not a lower triangular matrix ! Did you supply the upper Cholesky matrix instead ?'
                stop
            end if
        end if     
        
        !# Returns samples of the same size as qty from normal distribution with zero mean and unit variance N(0,1)
        allocate(normal_noise(Ncoef))
        call normal_vector(global_seed, RNG_FORECAST, global_i_real, i_t, normal_noise)
        
        !build a scaled noise
        !# Compute qty at t+1 from qty at t and the scaled noise (Euler-Maruyama)
        allocate(AR1_result(Ncoef), source=0.0d0)
        AR1_result = -MATMUL(X, A) + MATMUL(Chol, normal_noise)
    end subroutine
!==========================================================================================================================
    
!==========================================================================================================================
    subroutine ar3_process(X, A, B, C, Chol, global_seed, global_i_real, i_t, check_Cholesky, AR3_result)
    !"""
    !Applies an Auto-Regressive process of order 3 to the augmented state Z.
    !
    !X3 = - X2 @ A - X1 @ B - X0 @ C + Chol @ normal_noise
    !
    !:param X: quantity on which the AR-3 process will be applied (X0, X1, X2 above)
    !:type X: 2D numpy.ndarray (dim: 3 x Nz)
    !:param A: AR-3 operator (must be dense)
    !:type A: 2D numpy.ndarray (dim: Nz x Nz)
    !:param B: AR-3 operator (must be dense)
    !:type B: 2D numpy.ndarray (dim: Nz x Nz)
    !:param C: AR-3 operator (must be dense)
    !:type C: 2D numpy.ndarray (dim: Nz x Nz)
    !:param global_seed: global random seed
    !:type global_seed: int
    !:param global_i_real: global model realisation index
    !:type global_i_real: int
    !:param i_t: time index
    !:type i_t: int
    !:param check_Cholesky: If True, checks that the Cholesky matrix is lower triangular.
    !:type check_Cholesky: bool
    !:return: vector containing the result of the AR-3 process
    !:rtype: 1D numpy.ndarray (dim: Nz)
    !"""
        real(kind=8), intent(in) :: X(:,:)
        real(kind=8), intent(in) :: A(:,:), B(:,:), C(:,:), Chol(:,:)
        real(kind=8), allocatable, intent(out) :: AR3_result(:)
        integer, intent(in) :: global_seed, global_i_real, i_t
        logical, intent(in) :: check_Cholesky

        integer :: Ncoef, i, j
        real(kind=8), allocatable :: normal_noise(:), Chol_test(:,:)

        Ncoef = SIZE(X, 2)

        !# X must contain the three previous AR states
        if (SIZE(X, 1) .ne. 3) then
            write (10,'(A, i4, A)') &
                'First dimension of input quantity in AR-3 process should be 3. Got ', &
                SIZE(X, 1), ' instead.'
            write (*,'(A, i4, A)') &
                'First dimension of input quantity in AR-3 process should be 3. Got ', &
                SIZE(X, 1), ' instead.'
            stop
        end if

        !# Check dimensions of A
        if ((SIZE(A, 1) .ne. Ncoef) .or. (SIZE(A, 2) .ne. Ncoef)) then
            write (10,'(A, i4, A, i4, A, i4, A, i4, A)') &
                'A matrix in AR-3 should have dimensions (', Ncoef, ',', Ncoef, &
                '). Got (', SIZE(A, 1), ',', SIZE(A, 2), ') instead.'
            write (*,'(A, i4, A, i4, A, i4, A, i4, A)') &
                'A matrix in AR-3 should have dimensions (', Ncoef, ',', Ncoef, &
                '). Got (', SIZE(A, 1), ',', SIZE(A, 2), ') instead.'
            stop
        end if

        !# Check dimensions of B
        if ((SIZE(B, 1) .ne. Ncoef) .or. (SIZE(B, 2) .ne. Ncoef)) then
            write (10,'(A, i4, A, i4, A, i4, A, i4, A)') &
                'B matrix in AR-3 should have dimensions (', Ncoef, ',', Ncoef, &
                '). Got (', SIZE(B, 1), ',', SIZE(B, 2), ') instead.'
            write (*,'(A, i4, A, i4, A, i4, A, i4, A)') &
                'B matrix in AR-3 should have dimensions (', Ncoef, ',', Ncoef, &
                '). Got (', SIZE(B, 1), ',', SIZE(B, 2), ') instead.'
            stop
        end if

        !# Check dimensions of C
        if ((SIZE(C, 1) .ne. Ncoef) .or. (SIZE(C, 2) .ne. Ncoef)) then
            write (10,'(A, i4, A, i4, A, i4, A, i4, A)') &
                'C matrix in AR-3 should have dimensions (', Ncoef, ',', Ncoef, &
                '). Got (', SIZE(C, 1), ',', SIZE(C, 2), ') instead.'
            write (*,'(A, i4, A, i4, A, i4, A, i4, A)') &
                'C matrix in AR-3 should have dimensions (', Ncoef, ',', Ncoef, &
                '). Got (', SIZE(C, 1), ',', SIZE(C, 2), ') instead.'
            stop
        end if

        !# Check dimensions of Chol
        if ((SIZE(Chol, 1) .ne. Ncoef) .or. (SIZE(Chol, 2) .ne. Ncoef)) then
            write (10,'(A, i4, A, i4, A, i4, A, i4, A)') &
                'Cholesky matrix in AR-3 should have dimensions (', Ncoef, ',', Ncoef, &
                '). Got (', SIZE(Chol, 1), ',', SIZE(Chol, 2), ') instead.'
            write (*,'(A, i4, A, i4, A, i4, A, i4, A)') &
                'Cholesky matrix in AR-3 should have dimensions (', Ncoef, ',', Ncoef, &
                '). Got (', SIZE(Chol, 1), ',', SIZE(Chol, 2), ') instead.'
            stop
        end if

        !# Check that Chol is lower triangular
        if (check_Cholesky) then
            allocate(Chol_test, source=Chol)

            do i = 1, SIZE(Chol_test, 1)
                do j = i+1, SIZE(Chol_test, 2)
                    Chol_test(i,j) = 0.0d0
                end do
            end do

            if (ANY(ABS(Chol_test-Chol) > 1.0d-12)) then
                write (10,'(A)') &
                    'Cholesky matrix supplied in AR-3 is not a lower triangular matrix !'
                write (*,'(A)') &
                    'Cholesky matrix supplied in AR-3 is not a lower triangular matrix !'
                stop
            end if
        end if

        !# Returns samples from the normal distribution N(0,1)
        allocate(normal_noise(Ncoef))
        call normal_vector( &
            global_seed,       &
            RNG_FORECAST,      &
            global_i_real,     &
            i_t,               &
            normal_noise       &
        )

        !# Compute X3 from X0, X1 and X2 with the scaled noise
        allocate(AR3_result(Ncoef), source=0.0d0)

        AR3_result = -MATMUL(X(3,:), A) &
                     -MATMUL(X(2,:), B) &
                     -MATMUL(X(1,:), C) &
                   + MATMUL(Chol, normal_noise)      
    end subroutine
!==========================================================================================================================
    
!==========================================================================================================================
    subroutine cov(x, result_)
        ! 1/(N-1) (xt*x)
        real(kind=8), intent(in) :: x(:,:)
        real(kind=8), allocatable, intent(out) :: result_(:,:)
        real(kind=8), allocatable :: mean_X(:, :), l(:, :)
        real(kind=8), allocatable :: Y(:, :)
        real(kind=8) :: alpha
        
        allocate(result_(SIZE(x, 2), SIZE(x, 2)))
        allocate(mean_X(1, SIZE(x, 2)))
        allocate(l(1, SIZE(x, 1)))
        allocate(Y, source=x)
        !-----------------------juping remove mean value
        alpha = 1.0d0 / (real(SIZE(x, 1), kind=8))
        l     = 1.0d0
        mean_X = 0.0d0
        Y     = x
        call gemm(l, x, mean_X, 'N', 'N', alpha)
        call gemm(l*(-1.0d0), mean_X, Y, 'T', 'N', 1.0d0, 1.0d0)
        
        alpha = 1.0d0 / (real(SIZE(x, 1), kind=8) - 1.0d0)
        call gemm(Y, Y, result_, 'T', 'N', alpha)
    end subroutine
!==========================================================================================================================
end module
    