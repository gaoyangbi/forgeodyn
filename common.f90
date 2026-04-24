module common
    use utilities
    use linear_interpolation_module
    use blas95
    use lapack95
    use f95_precision
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
    subroutine compute_average(X, Nt, mode, result_)
    !# Perform a convolution of X by a blackman smoothing window of size Nt (NOT A PROPER AVERAGE)
        real(kind=8), intent(in) :: X(:,:)
        integer, intent(in) :: Nt
        character(len=*), intent(in) :: mode
        real(kind=8), allocatable, intent(out) :: result_(:,:)
        real(kind=8), allocatable :: X_avg(:,:)
        integer :: i, j
        real(kind=8), allocatable :: blackman_w(:)
        
        allocate(X_avg(SIZE(X, 1)-Nt+1, SIZE(X, 2)))
        allocate(blackman_w(Nt))
        allocate(result_(SIZE(X, 1)-Nt+1, SIZE(X, 2)))
        X_avg = 0.0d0
        call blackman(Nt, blackman_w)
        
        if (TRIM(mode) == 'valid') then
            do concurrent (i = 1: SIZE(X, 2))
                do concurrent (j = 1: SIZE(X_avg, 1))
                    X_avg(j, i) = DOT_PRODUCT(X(j:j+Nt-1,i), blackman_w)
                end do
            end do
        end if
        result_ = X_avg        
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
    