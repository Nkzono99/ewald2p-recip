module ewald2p_reciprocal_state
    use ewald2p_recip_types, only: dp, ip, ewald2p_recip_plan_type, ewald2p_recip_state_type
    use ewald2p_recip_utils, only: ERR_INVALID_OPTION, ERR_MEMORY_ALLOCATION, ERR_NON_NEUTRAL_SYSTEM, &
                                   ERR_PLAN_NOT_BUILT, ERR_SIZE_MISMATCH, charge_scale, clear_error, &
                                   set_error, total_charge_value
    use ewald2p_recip_wrap, only: wrap_periodic
    implicit none
    private

    public :: update_recip_state
    public :: destroy_recip_state
    public :: check_neutrality

contains

    subroutine check_neutrality(src_q, neutrality_tol, is_neutral, ierr, message, total_charge)
        real(dp), intent(in) :: src_q(:)
        real(dp), intent(in) :: neutrality_tol
        logical, intent(out) :: is_neutral
        integer(ip), intent(out) :: ierr
        character(len=*), intent(out), optional :: message
        real(dp), intent(out), optional :: total_charge

        real(dp) :: qsum
        real(dp) :: qscale

        call clear_error(ierr, message)

        if (neutrality_tol < 0.0_dp) then
            is_neutral = .false.
            if (present(total_charge)) total_charge = 0.0_dp
            call set_error(ierr, ERR_INVALID_OPTION, "Invalid neutrality tolerance: must be nonnegative.", message)
            return
        end if

        qsum = total_charge_value(src_q)
        qscale = charge_scale(src_q)

        if (present(total_charge)) total_charge = qsum
        is_neutral = abs(qsum) <= neutrality_tol * qscale
    end subroutine check_neutrality

    subroutine update_recip_state(plan, state, src_pos, src_q, ierr, message)
        type(ewald2p_recip_plan_type), intent(in) :: plan
        type(ewald2p_recip_state_type), intent(inout) :: state
        real(dp), intent(in) :: src_pos(:, :)
        real(dp), intent(in) :: src_q(:)
        integer(ip), intent(out) :: ierr
        character(len=*), intent(out), optional :: message

        integer(ip) :: nsrc
        integer(ip) :: ik
        integer(ip) :: isrc
        integer :: stat
        logical :: is_neutral
        real(dp) :: phase
        real(dp) :: qsum
        complex(dp) :: moment_plus
        complex(dp) :: moment_minus
        character(len=128) :: local_message

        call clear_error(ierr, message)

        if (.not. plan%is_built) then
            call set_error(ierr, ERR_PLAN_NOT_BUILT, "Plan is not built.", message)
            return
        end if

        if (size(src_pos, 1) /= 3) then
            call set_error(ierr, ERR_SIZE_MISMATCH, "src_pos must have shape (3, nsrc).", message)
            return
        end if

        nsrc = int(size(src_pos, 2), kind=ip)
        if (size(src_q) /= nsrc) then
            call set_error(ierr, ERR_SIZE_MISMATCH, "src_q size must match src_pos second dimension.", message)
            return
        end if

        call check_neutrality(src_q, plan%neutrality_tol, is_neutral, ierr, message, total_charge=qsum)
        if (ierr /= 0_ip) return

        if (.not. is_neutral) then
            write (local_message, '(a,es24.16)') "Non-neutral system: total charge = ", qsum
            call set_error(ierr, ERR_NON_NEUTRAL_SYSTEM, trim(local_message), message)
            return
        end if

        call destroy_recip_state(state)

        allocate (state%src_pos(3, nsrc), state%src_q(nsrc), state%s_plus(plan%nk), state%s_minus(plan%nk), stat=stat)
        if (stat /= 0) then
            call set_error(ierr, ERR_MEMORY_ALLOCATION, "Failed to allocate reciprocal state arrays.", message)
            call destroy_recip_state(state)
            return
        end if

        do isrc = 1, nsrc
            state%src_pos(1, isrc) = wrap_periodic(src_pos(1, isrc), plan%lx)
            state%src_pos(2, isrc) = wrap_periodic(src_pos(2, isrc), plan%ly)
            state%src_pos(3, isrc) = src_pos(3, isrc)
            state%src_q(isrc) = src_q(isrc)
        end do

        do ik = 1, plan%nk
            moment_plus = cmplx(0.0_dp, 0.0_dp, kind=dp)
            moment_minus = cmplx(0.0_dp, 0.0_dp, kind=dp)

            do isrc = 1, nsrc
                phase = plan%kx(ik) * state%src_pos(1, isrc) + plan%ky(ik) * state%src_pos(2, isrc)
                moment_plus = moment_plus + state%src_q(isrc) * cmplx(cos(phase), -sin(phase), kind=dp)
                moment_minus = moment_minus + state%src_q(isrc) * cmplx(cos(phase), sin(phase), kind=dp)
            end do

            state%s_plus(ik) = moment_plus
            state%s_minus(ik) = moment_minus
        end do

        state%nsrc = nsrc
        state%is_ready = .true.
    end subroutine update_recip_state

    subroutine destroy_recip_state(state)
        type(ewald2p_recip_state_type), intent(inout) :: state

        if (allocated(state%src_pos)) deallocate (state%src_pos)
        if (allocated(state%src_q)) deallocate (state%src_q)
        if (allocated(state%s_plus)) deallocate (state%s_plus)
        if (allocated(state%s_minus)) deallocate (state%s_minus)

        state%is_ready = .false.
        state%nsrc = 0_ip
    end subroutine destroy_recip_state
end module ewald2p_reciprocal_state
