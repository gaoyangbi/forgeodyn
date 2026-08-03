module deterministic_rng

    use, intrinsic :: iso_fortran_env, only: int64, real64
    use mkl_vsl_type
    use mkl_vsl

    implicit none
    private

    ! Different scientific uses must have different stream identifiers.
    integer, parameter, public :: RNG_INIT_B    = 1
    integer, parameter, public :: RNG_INIT_Z    = 2
    integer, parameter, public :: RNG_FORECAST  = 3
    integer, parameter, public :: RNG_INIT_DZ   = 4
    integer, parameter, public :: RNG_INIT_D2Z  = 5

    ! Reserve a large counter block for every
    ! (stream_id, global_i_real, i_t) tuple.
    integer(int64), parameter :: BLOCK_STRIDE = 1048576_int64

    public :: normal_vector

contains

    subroutine normal_vector( &
        global_seed,       &
        stream_id,         &
        global_i_real,     &
        i_t,               &
        normal_noise       &
    )

        integer, intent(in) :: global_seed
        integer, intent(in) :: stream_id
        integer, intent(in) :: global_i_real
        integer, intent(in) :: i_t
        real(real64), intent(out) :: normal_noise(:)

        type(VSL_STREAM_STATE) :: rng_stream

        integer :: status
        integer :: seed32
        integer(int64) :: key
        integer(int64) :: offset

        if (stream_id < 0) then
            error stop "normal_vector: stream_id must be non-negative"
        end if

        if (global_i_real < 0) then
            error stop "normal_vector: global_i_real must be non-negative"
        end if

        if (i_t < 0) then
            error stop "normal_vector: i_t must be non-negative"
        end if

        if (size(normal_noise) > BLOCK_STRIDE) then
            error stop "normal_vector: noise vector exceeds RNG block"
        end if

        ! MKL vslNewStream accepts a default 32-bit integer seed.
        ! Convert any global seed to the positive interval
        ! [1, 2147483646].
        seed32 = int( &
            modulo( &
                int(global_seed, int64), &
                2147483646_int64 &
            ) + 1_int64 &
        )

        ! Map the tuple
        ! (stream_id, global_i_real, i_t)
        ! to one deterministic non-negative integer.
        key = cantor_pair( &
            cantor_pair( &
                int(stream_id, int64), &
                int(global_i_real, int64) &
            ), &
            int(i_t, int64) &
        )

        offset = key * BLOCK_STRIDE

        status = vslNewStream( &
            rng_stream,                &
            VSL_BRNG_PHILOX4X32X10,    &
            seed32                     &
        )
        call check_vsl_status(status, "vslNewStream")

        status = vslSkipAheadStream(rng_stream, offset)
        call check_vsl_status(status, "vslSkipAheadStream")

        status = vdRngGaussian( &
            VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, &
            rng_stream,                         &
            size(normal_noise),                  &
            normal_noise,                       &
            0.0_real64,                         &
            1.0_real64                          &
        )
        call check_vsl_status(status, "vdRngGaussian")

        status = vslDeleteStream(rng_stream)
        call check_vsl_status(status, "vslDeleteStream")

    end subroutine normal_vector


    pure function cantor_pair(a, b) result(paired_value)

        integer(int64), intent(in) :: a
        integer(int64), intent(in) :: b
        integer(int64) :: paired_value
        integer(int64) :: sum_ab

        sum_ab = a + b
        paired_value = sum_ab * (sum_ab + 1_int64) / 2_int64 + b

    end function cantor_pair


    subroutine check_vsl_status(status, routine_name)

        integer, intent(in) :: status
        character(len=*), intent(in) :: routine_name

        if (status /= VSL_STATUS_OK) then
            write(*, '(A, A, A, I0)') &
                "MKL VSL error in ", trim(routine_name), &
                ", status = ", status

            error stop "MKL VSL random number generation failed"
        end if

    end subroutine check_vsl_status

end module deterministic_rng
