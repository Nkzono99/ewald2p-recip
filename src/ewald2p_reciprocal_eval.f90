module ewald2p_reciprocal_eval
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use ewald2p_recip_kspace, only: nonzero_mode_kernel
    use ewald2p_recip_types, only: dp, ip, ewald2p_recip_plan_type, ewald2p_recip_state_type
    use ewald2p_recip_utils, only: ERR_INTERNAL_NUMERICAL, ERR_PLAN_NOT_BUILT, ERR_SIZE_MISMATCH, &
                                   ERR_STATE_NOT_READY, coincidence_tol, clear_error, pi_dp, set_error
    use ewald2p_recip_wrap, only: wrap_xy_point
    use ewald2p_recip_zero_mode, only: eval_zero_mode_contribution, self_potential_scale
    implicit none
    private

    public :: eval_recip_point
    public :: eval_recip_points

contains

    subroutine eval_recip_point(plan, state, r, e, phi, ierr, message)
        type(ewald2p_recip_plan_type), intent(in) :: plan
        type(ewald2p_recip_state_type), intent(in) :: state
        real(dp), intent(in) :: r(3)
        real(dp), intent(out) :: e(3)
        real(dp), intent(out), optional :: phi
        integer(ip), intent(out) :: ierr
        character(len=*), intent(out), optional :: message

        real(dp) :: phi_local

        call clear_error(ierr, message)

        if (.not. plan%is_built) then
            e = 0.0_dp
            if (present(phi)) phi = 0.0_dp
            call set_error(ierr, ERR_PLAN_NOT_BUILT, "Plan is not built.", message)
            return
        end if

        if (.not. state%is_ready) then
            e = 0.0_dp
            if (present(phi)) phi = 0.0_dp
            call set_error(ierr, ERR_STATE_NOT_READY, "State is not ready.", message)
            return
        end if

        call evaluate_single_target(plan, state, r, e, phi_local, present(phi), ierr, message)
        if (ierr /= 0_ip) return

        if (present(phi)) phi = phi_local
    end subroutine eval_recip_point

    subroutine eval_recip_points(plan, state, target_pos, e, phi, ierr, message)
        type(ewald2p_recip_plan_type), intent(in) :: plan
        type(ewald2p_recip_state_type), intent(in) :: state
        real(dp), intent(in) :: target_pos(:, :)
        real(dp), intent(out) :: e(:, :)
        real(dp), intent(out), optional :: phi(:)
        integer(ip), intent(out) :: ierr
        character(len=*), intent(out), optional :: message

        integer(ip) :: ntgt
        integer(ip) :: itgt
        real(dp) :: phi_local

        call clear_error(ierr, message)

        if (.not. plan%is_built) then
            e = 0.0_dp
            if (present(phi)) phi = 0.0_dp
            call set_error(ierr, ERR_PLAN_NOT_BUILT, "Plan is not built.", message)
            return
        end if

        if (.not. state%is_ready) then
            e = 0.0_dp
            if (present(phi)) phi = 0.0_dp
            call set_error(ierr, ERR_STATE_NOT_READY, "State is not ready.", message)
            return
        end if

        if (size(target_pos, 1) /= 3) then
            e = 0.0_dp
            if (present(phi)) phi = 0.0_dp
            call set_error(ierr, ERR_SIZE_MISMATCH, "target_pos must have shape (3, ntgt).", message)
            return
        end if

        ntgt = int(size(target_pos, 2), kind=ip)

        if (size(e, 1) /= 3 .or. size(e, 2) /= ntgt) then
            e = 0.0_dp
            if (present(phi)) phi = 0.0_dp
            call set_error(ierr, ERR_SIZE_MISMATCH, "e must have shape (3, ntgt).", message)
            return
        end if

        if (present(phi)) then
            if (size(phi) /= ntgt) then
                e = 0.0_dp
                phi = 0.0_dp
                call set_error(ierr, ERR_SIZE_MISMATCH, "phi must have shape (ntgt).", message)
                return
            end if
        end if

        do itgt = 1, ntgt
            call evaluate_single_target(plan, state, target_pos(:, itgt), e(:, itgt), phi_local, present(phi), ierr, message)
            if (ierr /= 0_ip) then
                if (present(phi)) phi(itgt:) = 0.0_dp
                if (itgt < ntgt) e(:, itgt:ntgt) = 0.0_dp
                return
            end if

            if (present(phi)) phi(itgt) = phi_local
        end do
    end subroutine eval_recip_points

    subroutine evaluate_single_target(plan, state, r_in, e, phi, need_phi, ierr, message)
        type(ewald2p_recip_plan_type), intent(in) :: plan
        type(ewald2p_recip_state_type), intent(in) :: state
        real(dp), intent(in) :: r_in(3)
        real(dp), intent(out) :: e(3)
        real(dp), intent(out) :: phi
        logical, intent(in) :: need_phi
        integer(ip), intent(out) :: ierr
        character(len=*), intent(out), optional :: message

        integer(ip) :: ik
        integer(ip) :: isrc
        real(dp) :: area
        real(dp) :: prefactor
        real(dp) :: phase
        real(dp) :: cos_phase
        real(dp) :: sin_phase
        real(dp) :: dz
        real(dp) :: kernel_sum
        real(dp) :: kernel_diff
        real(dp) :: phi0
        real(dp) :: ez0
        real(dp) :: r(3)

        call clear_error(ierr, message)

        call wrap_xy_point(plan%lx, plan%ly, r_in, r)

        area = plan%lx * plan%ly
        prefactor = pi_dp() / area

        e = 0.0_dp
        phi = 0.0_dp

        do ik = 1, plan%nk
            do isrc = 1, state%nsrc
                dz = r(3) - state%src_pos(3, isrc)
                call nonzero_mode_kernel(plan%xi, plan%kabs(ik), dz, kernel_sum, kernel_diff)

                phase = plan%kx(ik) * (r(1) - state%src_pos(1, isrc)) + &
                        plan%ky(ik) * (r(2) - state%src_pos(2, isrc))
                cos_phase = cos(phase)
                sin_phase = sin(phase)

                if (need_phi) then
                    phi = phi + prefactor * state%src_q(isrc) * kernel_sum * cos_phase / plan%kabs(ik)
                end if

                e(1) = e(1) + prefactor * state%src_q(isrc) * (plan%kx(ik) / plan%kabs(ik)) * kernel_sum * sin_phase
                e(2) = e(2) + prefactor * state%src_q(isrc) * (plan%ky(ik) / plan%kabs(ik)) * kernel_sum * sin_phase
                e(3) = e(3) - prefactor * state%src_q(isrc) * kernel_diff * cos_phase
            end do
        end do

        call eval_zero_mode_contribution(plan%lx, plan%ly, plan%xi, state%src_pos(3, :), state%src_q, r(3), phi0, ez0)

        if (need_phi) phi = phi + phi0
        e(3) = e(3) + ez0

        if (need_phi) phi = phi + self_correction_if_needed(plan, state, r)

        if (.not. all(ieee_is_finite(e))) then
            call set_error(ierr, ERR_INTERNAL_NUMERICAL, "Non-finite field encountered during evaluation.", message)
            return
        end if

        if (need_phi) then
            if (.not. ieee_is_finite(phi)) then
                call set_error(ierr, ERR_INTERNAL_NUMERICAL, "Non-finite potential encountered during evaluation.", message)
                return
            end if
        end if
    end subroutine evaluate_single_target

    pure real(dp) function self_correction_if_needed(plan, state, r)
        type(ewald2p_recip_plan_type), intent(in) :: plan
        type(ewald2p_recip_state_type), intent(in) :: state
        real(dp), intent(in) :: r(3)

        integer(ip) :: isrc
        real(dp) :: tol

        self_correction_if_needed = 0.0_dp
        tol = coincidence_tol(plan%lx, plan%ly, r(3))

        do isrc = 1, state%nsrc
            if (abs(r(1) - state%src_pos(1, isrc)) <= tol .and. &
                abs(r(2) - state%src_pos(2, isrc)) <= tol .and. &
                abs(r(3) - state%src_pos(3, isrc)) <= tol) then
                self_correction_if_needed = self_correction_if_needed + self_potential_scale(plan%xi) * state%src_q(isrc)
            end if
        end do
    end function self_correction_if_needed
end module ewald2p_reciprocal_eval
