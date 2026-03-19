program test_ewald2p_reciprocal
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use ewald2p_reciprocal, only: build_recip_plan, check_neutrality, destroy_recip_plan, destroy_recip_state, &
                                  dp, ewald2p_recip_options_type, ewald2p_recip_plan_type, &
                                  ewald2p_recip_state_type, eval_recip_point, eval_recip_points, ip, &
                                  update_recip_state
    implicit none

    call test_build_and_neutrality()
    call test_periodicity_and_batch()
    call test_zero_mode_only()
    call test_self_policy()
    call test_cutoff_convergence()

contains

    subroutine test_build_and_neutrality()
        type(ewald2p_recip_options_type) :: options
        type(ewald2p_recip_plan_type) :: plan
        type(ewald2p_recip_state_type) :: state
        real(dp) :: src_pos_single(3, 1)
        real(dp) :: src_q_single(1)
        real(dp) :: src_pos(3, 2)
        real(dp) :: src_q(2)
        real(dp) :: r(3)
        real(dp) :: e(3)
        real(dp) :: phi
        integer(ip) :: ierr
        logical :: is_neutral
        character(len=256) :: message

        call fill_options(options, 1.5_dp, 2.0_dp, 1.1_dp, 2_ip, 2_ip)

        call build_recip_plan(plan, options, ierr, message)
        call require_success(ierr, message, "build_recip_plan basic")
        call assert_true(plan%is_built, "plan should be built")
        call assert_true(plan%nk > 0_ip, "plan should contain nonzero modes")

        src_pos_single(:, 1) = [0.1_dp, -0.2_dp, -0.3_dp]
        src_q_single = [1.0_dp]
        call update_recip_state(plan, state, src_pos_single, src_q_single, ierr, message)
        call assert_true(ierr == 5_ip, "single source must be rejected as non-neutral")

        src_pos(:, 1) = [0.1_dp, -0.2_dp, -0.3_dp]
        src_pos(:, 2) = [1.0_dp, 1.6_dp, 0.4_dp]
        src_q = [1.0_dp, -1.0_dp]

        call check_neutrality(src_q, options%neutrality_tol, is_neutral, ierr, message)
        call require_success(ierr, message, "check_neutrality basic")
        call assert_true(is_neutral, "two-source system should be neutral")

        call update_recip_state(plan, state, src_pos, src_q, ierr, message)
        call require_success(ierr, message, "update_recip_state basic")
        call assert_true(state%is_ready, "state should be ready after update")

        r = [0.25_dp, 0.75_dp, 0.2_dp]
        call eval_recip_point(plan, state, r, e, phi, ierr, message)
        call require_success(ierr, message, "eval_recip_point basic")
        call assert_finite_vector(e, "field should be finite")
        call assert_finite_scalar(phi, "potential should be finite")

        call destroy_recip_state(state)
        call destroy_recip_plan(plan)
    end subroutine test_build_and_neutrality

    subroutine test_periodicity_and_batch()
        type(ewald2p_recip_options_type) :: options
        type(ewald2p_recip_plan_type) :: plan
        type(ewald2p_recip_state_type) :: state
        real(dp) :: src_pos(3, 3)
        real(dp) :: src_q(3)
        real(dp) :: r(3)
        real(dp) :: r_x(3)
        real(dp) :: r_y(3)
        real(dp) :: e_ref(3)
        real(dp) :: e_x(3)
        real(dp) :: e_y(3)
        real(dp) :: phi_ref
        real(dp) :: phi_x
        real(dp) :: phi_y
        real(dp) :: targets(3, 3)
        real(dp) :: e_batch(3, 3)
        real(dp) :: phi_batch(3)
        real(dp) :: e_single(3)
        real(dp) :: phi_single
        integer(ip) :: ierr
        integer(ip) :: itgt
        character(len=256) :: message

        call fill_options(options, 1.4_dp, 1.8_dp, 1.3_dp, 3_ip, 3_ip)

        src_pos(:, 1) = [0.15_dp, 0.20_dp, -0.40_dp]
        src_pos(:, 2) = [0.70_dp, 1.10_dp, 0.10_dp]
        src_pos(:, 3) = [-0.30_dp, 0.55_dp, 0.60_dp]
        src_q = [1.0_dp, -0.4_dp, -0.6_dp]

        call build_recip_plan(plan, options, ierr, message)
        call require_success(ierr, message, "build_recip_plan periodicity")
        call update_recip_state(plan, state, src_pos, src_q, ierr, message)
        call require_success(ierr, message, "update_recip_state periodicity")

        r = [0.37_dp, -0.81_dp, 0.25_dp]
        r_x = r + [options%lx, 0.0_dp, 0.0_dp]
        r_y = r + [0.0_dp, options%ly, 0.0_dp]

        call eval_recip_point(plan, state, r, e_ref, phi_ref, ierr, message)
        call require_success(ierr, message, "eval_recip_point reference periodicity")
        call eval_recip_point(plan, state, r_x, e_x, phi_x, ierr, message)
        call require_success(ierr, message, "eval_recip_point x-periodicity")
        call eval_recip_point(plan, state, r_y, e_y, phi_y, ierr, message)
        call require_success(ierr, message, "eval_recip_point y-periodicity")

        call assert_close_vector(e_x, e_ref, 1.0e-11_dp, "x periodicity field")
        call assert_close_vector(e_y, e_ref, 1.0e-11_dp, "y periodicity field")
        call assert_close_scalar(phi_x, phi_ref, 1.0e-11_dp, "x periodicity potential")
        call assert_close_scalar(phi_y, phi_ref, 1.0e-11_dp, "y periodicity potential")

        targets(:, 1) = r
        targets(:, 2) = r_x
        targets(:, 3) = [0.88_dp, 1.25_dp, -0.35_dp]

        call eval_recip_points(plan, state, targets, e_batch, phi_batch, ierr, message)
        call require_success(ierr, message, "eval_recip_points batch")

        do itgt = 1, 3
            call eval_recip_point(plan, state, targets(:, itgt), e_single, phi_single, ierr, message)
            call require_success(ierr, message, "eval_recip_point batch comparison")
            call assert_close_vector(e_batch(:, itgt), e_single, 1.0e-12_dp, "point/batch field consistency")
            call assert_close_scalar(phi_batch(itgt), phi_single, 1.0e-12_dp, "point/batch potential consistency")
        end do

        call destroy_recip_state(state)
        call destroy_recip_plan(plan)
    end subroutine test_periodicity_and_batch

    subroutine test_zero_mode_only()
        type(ewald2p_recip_options_type) :: options
        type(ewald2p_recip_plan_type) :: plan
        type(ewald2p_recip_state_type) :: state
        real(dp) :: src_pos(3, 2)
        real(dp) :: src_q(2)
        real(dp) :: r(3)
        real(dp) :: e(3)
        real(dp) :: phi
        real(dp) :: phi_ref
        real(dp) :: ez_ref
        integer(ip) :: ierr
        character(len=256) :: message

        call fill_options(options, 1.2_dp, 1.6_dp, 1.7_dp, 0_ip, 0_ip)

        src_pos(:, 1) = [0.2_dp, 0.4_dp, -0.5_dp]
        src_pos(:, 2) = [0.9_dp, 1.1_dp, 0.3_dp]
        src_q = [1.0_dp, -1.0_dp]

        call build_recip_plan(plan, options, ierr, message)
        call require_success(ierr, message, "build_recip_plan zero mode")
        call assert_true(plan%nk == 0_ip, "h_max=l_max=0 should give nk=0")

        call update_recip_state(plan, state, src_pos, src_q, ierr, message)
        call require_success(ierr, message, "update_recip_state zero mode")

        r = [0.45_dp, 0.35_dp, 0.15_dp]
        call eval_recip_point(plan, state, r, e, phi, ierr, message)
        call require_success(ierr, message, "eval_recip_point zero mode only")

        call reference_zero_mode(options%lx, options%ly, options%xi, state%src_pos(3, :), state%src_q, r(3), phi_ref, ez_ref)

        call assert_close_scalar(e(1), 0.0_dp, 1.0e-13_dp, "zero mode ex should vanish")
        call assert_close_scalar(e(2), 0.0_dp, 1.0e-13_dp, "zero mode ey should vanish")
        call assert_close_scalar(e(3), ez_ref, 1.0e-12_dp, "zero mode ez")
        call assert_close_scalar(phi, phi_ref, 1.0e-12_dp, "zero mode phi")

        call destroy_recip_state(state)
        call destroy_recip_plan(plan)
    end subroutine test_zero_mode_only

    subroutine test_self_policy()
        type(ewald2p_recip_options_type) :: options
        type(ewald2p_recip_plan_type) :: plan
        type(ewald2p_recip_state_type) :: state
        real(dp) :: src_pos(3, 2)
        real(dp) :: src_q(2)
        real(dp) :: e(3)
        real(dp) :: phi
        real(dp) :: e_ref(3)
        real(dp) :: phi_ref
        integer(ip) :: ierr
        character(len=256) :: message

        call fill_options(options, 1.8_dp, 1.4_dp, 1.5_dp, 4_ip, 4_ip)

        src_pos(:, 1) = [0.20_dp, 0.30_dp, 0.10_dp]
        src_pos(:, 2) = [1.10_dp, 1.00_dp, -0.45_dp]
        src_q = [1.0_dp, -1.0_dp]

        call build_recip_plan(plan, options, ierr, message)
        call require_success(ierr, message, "build_recip_plan self policy")
        call update_recip_state(plan, state, src_pos, src_q, ierr, message)
        call require_success(ierr, message, "update_recip_state self policy")

        call eval_recip_point(plan, state, src_pos(:, 1), e, phi, ierr, message)
        call require_success(ierr, message, "eval_recip_point self policy")

        call reference_eval_no_self(plan, state, src_pos(:, 1), e_ref, phi_ref)

        call assert_close_vector(e, e_ref, 1.0e-12_dp, "self policy should not alter field")
        call assert_close_scalar(phi, phi_ref + self_scale_local(options%xi) * src_q(1), 1.0e-12_dp, &
                                 "self policy should subtract standard potential self term")

        call destroy_recip_state(state)
        call destroy_recip_plan(plan)
    end subroutine test_self_policy

    subroutine test_cutoff_convergence()
        type(ewald2p_recip_options_type) :: options_low
        type(ewald2p_recip_options_type) :: options_mid
        type(ewald2p_recip_options_type) :: options_hi
        type(ewald2p_recip_plan_type) :: plan_low
        type(ewald2p_recip_plan_type) :: plan_mid
        type(ewald2p_recip_plan_type) :: plan_hi
        type(ewald2p_recip_state_type) :: state_low
        type(ewald2p_recip_state_type) :: state_mid
        type(ewald2p_recip_state_type) :: state_hi
        real(dp) :: src_pos(3, 4)
        real(dp) :: src_q(4)
        real(dp) :: r(3)
        real(dp) :: e_low(3)
        real(dp) :: e_mid(3)
        real(dp) :: e_hi(3)
        real(dp) :: phi_low
        real(dp) :: phi_mid
        real(dp) :: phi_hi
        real(dp) :: err_low
        real(dp) :: err_mid
        integer(ip) :: ierr
        character(len=256) :: message

        call fill_options(options_low, 1.6_dp, 1.9_dp, 1.2_dp, 2_ip, 2_ip)
        call fill_options(options_mid, 1.6_dp, 1.9_dp, 1.2_dp, 4_ip, 4_ip)
        call fill_options(options_hi, 1.6_dp, 1.9_dp, 1.2_dp, 8_ip, 8_ip)

        src_pos(:, 1) = [0.10_dp, 0.20_dp, -0.60_dp]
        src_pos(:, 2) = [0.55_dp, 1.10_dp, -0.10_dp]
        src_pos(:, 3) = [1.25_dp, 0.75_dp, 0.35_dp]
        src_pos(:, 4) = [-0.30_dp, 1.55_dp, 0.70_dp]
        src_q = [1.2_dp, -0.5_dp, -0.4_dp, -0.3_dp]
        r = [0.62_dp, 0.48_dp, 0.18_dp]

        call build_recip_plan(plan_low, options_low, ierr, message)
        call require_success(ierr, message, "build_recip_plan low cutoff")
        call update_recip_state(plan_low, state_low, src_pos, src_q, ierr, message)
        call require_success(ierr, message, "update_recip_state low cutoff")
        call eval_recip_point(plan_low, state_low, r, e_low, phi_low, ierr, message)
        call require_success(ierr, message, "eval low cutoff")

        call build_recip_plan(plan_mid, options_mid, ierr, message)
        call require_success(ierr, message, "build_recip_plan mid cutoff")
        call update_recip_state(plan_mid, state_mid, src_pos, src_q, ierr, message)
        call require_success(ierr, message, "update_recip_state mid cutoff")
        call eval_recip_point(plan_mid, state_mid, r, e_mid, phi_mid, ierr, message)
        call require_success(ierr, message, "eval mid cutoff")

        call build_recip_plan(plan_hi, options_hi, ierr, message)
        call require_success(ierr, message, "build_recip_plan high cutoff")
        call update_recip_state(plan_hi, state_hi, src_pos, src_q, ierr, message)
        call require_success(ierr, message, "update_recip_state high cutoff")
        call eval_recip_point(plan_hi, state_hi, r, e_hi, phi_hi, ierr, message)
        call require_success(ierr, message, "eval high cutoff")

        err_low = combined_error(e_low, phi_low, e_hi, phi_hi)
        err_mid = combined_error(e_mid, phi_mid, e_hi, phi_hi)
        call assert_true(err_mid < err_low, "larger cutoff should improve agreement with a higher-cutoff reference")

        call destroy_recip_state(state_low)
        call destroy_recip_state(state_mid)
        call destroy_recip_state(state_hi)
        call destroy_recip_plan(plan_low)
        call destroy_recip_plan(plan_mid)
        call destroy_recip_plan(plan_hi)
    end subroutine test_cutoff_convergence

    subroutine fill_options(options, lx, ly, xi, h_max, l_max)
        type(ewald2p_recip_options_type), intent(out) :: options
        real(dp), intent(in) :: lx
        real(dp), intent(in) :: ly
        real(dp), intent(in) :: xi
        integer(ip), intent(in) :: h_max
        integer(ip), intent(in) :: l_max

        options%lx = lx
        options%ly = ly
        options%xi = xi
        options%h_max = h_max
        options%l_max = l_max
        options%neutrality_tol = 1.0e-12_dp
        options%return_potential = .true.
        options%enable_zero_mode = .true.
    end subroutine fill_options

    subroutine require_success(ierr, message, label)
        integer(ip), intent(in) :: ierr
        character(len=*), intent(in) :: message
        character(len=*), intent(in) :: label

        if (ierr /= 0_ip) then
            call fail(label // ": " // trim(message))
        end if
    end subroutine require_success

    subroutine assert_true(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label

        if (.not. condition) call fail(label)
    end subroutine assert_true

    subroutine assert_close_scalar(actual, expected, tol, label)
        real(dp), intent(in) :: actual
        real(dp), intent(in) :: expected
        real(dp), intent(in) :: tol
        character(len=*), intent(in) :: label

        real(dp) :: scale

        scale = max(1.0_dp, abs(actual))
        scale = max(scale, abs(expected))

        if (abs(actual - expected) > tol * scale) then
            call fail(label)
        end if
    end subroutine assert_close_scalar

    subroutine assert_close_vector(actual, expected, tol, label)
        real(dp), intent(in) :: actual(3)
        real(dp), intent(in) :: expected(3)
        real(dp), intent(in) :: tol
        character(len=*), intent(in) :: label

        real(dp) :: scale

        scale = max(1.0_dp, maxval(abs(actual)))
        scale = max(scale, maxval(abs(expected)))

        if (maxval(abs(actual - expected)) > tol * scale) then
            call fail(label)
        end if
    end subroutine assert_close_vector

    subroutine assert_finite_scalar(value, label)
        real(dp), intent(in) :: value
        character(len=*), intent(in) :: label

        if (.not. ieee_is_finite(value)) call fail(label)
    end subroutine assert_finite_scalar

    subroutine assert_finite_vector(value, label)
        real(dp), intent(in) :: value(3)
        character(len=*), intent(in) :: label

        if (.not. all(ieee_is_finite(value))) call fail(label)
    end subroutine assert_finite_vector

    subroutine fail(message)
        character(len=*), intent(in) :: message

        write (*, '(a)') trim(message)
        error stop 1
    end subroutine fail

    subroutine reference_eval_no_self(plan, state, r_in, e, phi)
        type(ewald2p_recip_plan_type), intent(in) :: plan
        type(ewald2p_recip_state_type), intent(in) :: state
        real(dp), intent(in) :: r_in(3)
        real(dp), intent(out) :: e(3)
        real(dp), intent(out) :: phi

        integer(ip) :: ik
        integer(ip) :: isrc
        real(dp) :: r(3)
        real(dp) :: prefactor
        real(dp) :: phase
        real(dp) :: cos_phase
        real(dp) :: sin_phase
        real(dp) :: dz
        real(dp) :: kernel_sum
        real(dp) :: kernel_diff
        real(dp) :: phi0
        real(dp) :: ez0

        call wrap_xy_local(plan%lx, plan%ly, r_in, r)

        prefactor = pi_local() / (plan%lx * plan%ly)
        e = 0.0_dp
        phi = 0.0_dp

        do ik = 1, plan%nk
            do isrc = 1, state%nsrc
                dz = r(3) - state%src_pos(3, isrc)
                call nonzero_kernel_local(plan%xi, plan%kabs(ik), dz, kernel_sum, kernel_diff)

                phase = plan%kx(ik) * (r(1) - state%src_pos(1, isrc)) + &
                        plan%ky(ik) * (r(2) - state%src_pos(2, isrc))
                cos_phase = cos(phase)
                sin_phase = sin(phase)

                phi = phi + prefactor * state%src_q(isrc) * kernel_sum * cos_phase / plan%kabs(ik)
                e(1) = e(1) + prefactor * state%src_q(isrc) * (plan%kx(ik) / plan%kabs(ik)) * kernel_sum * sin_phase
                e(2) = e(2) + prefactor * state%src_q(isrc) * (plan%ky(ik) / plan%kabs(ik)) * kernel_sum * sin_phase
                e(3) = e(3) - prefactor * state%src_q(isrc) * kernel_diff * cos_phase
            end do
        end do

        call reference_zero_mode(plan%lx, plan%ly, plan%xi, state%src_pos(3, :), state%src_q, r(3), phi0, ez0)
        phi = phi + phi0
        e(3) = e(3) + ez0
    end subroutine reference_eval_no_self

    subroutine reference_zero_mode(lx, ly, xi, src_z, src_q, z, phi0, ez0)
        real(dp), intent(in) :: lx
        real(dp), intent(in) :: ly
        real(dp), intent(in) :: xi
        real(dp), intent(in) :: src_z(:)
        real(dp), intent(in) :: src_q(:)
        real(dp), intent(in) :: z
        real(dp), intent(out) :: phi0
        real(dp), intent(out) :: ez0

        integer :: isrc
        real(dp) :: sqrt_pi
        real(dp) :: dz
        real(dp) :: xi_dz

        sqrt_pi = sqrt(pi_local())
        phi0 = 0.0_dp
        ez0 = 0.0_dp

        do isrc = 1, size(src_q)
            dz = z - src_z(isrc)
            xi_dz = xi * dz
            phi0 = phi0 - (2.0_dp * sqrt_pi / (lx * ly)) * src_q(isrc) * &
                   ((exp(-(xi_dz * xi_dz)) / xi) + sqrt_pi * dz * erf(xi_dz))
            ez0 = ez0 + (2.0_dp * pi_local() / (lx * ly)) * src_q(isrc) * erf(xi_dz)
        end do
    end subroutine reference_zero_mode

    pure subroutine nonzero_kernel_local(xi, kabs, dz, kernel_sum, kernel_diff)
        real(dp), intent(in) :: xi
        real(dp), intent(in) :: kabs
        real(dp), intent(in) :: dz
        real(dp), intent(out) :: kernel_sum
        real(dp), intent(out) :: kernel_diff

        real(dp) :: a
        real(dp) :: term_plus
        real(dp) :: term_minus

        a = kabs / (2.0_dp * xi)
        term_plus = exp(kabs * dz) * erfc(a + xi * dz)
        term_minus = exp(-kabs * dz) * erfc(a - xi * dz)

        kernel_sum = term_plus + term_minus
        kernel_diff = term_plus - term_minus
    end subroutine nonzero_kernel_local

    pure subroutine wrap_xy_local(lx, ly, r_in, r_out)
        real(dp), intent(in) :: lx
        real(dp), intent(in) :: ly
        real(dp), intent(in) :: r_in(3)
        real(dp), intent(out) :: r_out(3)

        r_out = r_in
        r_out(1) = modulo(r_in(1), lx)
        r_out(2) = modulo(r_in(2), ly)

        if (r_out(1) < 0.0_dp) r_out(1) = r_out(1) + lx
        if (r_out(2) < 0.0_dp) r_out(2) = r_out(2) + ly
    end subroutine wrap_xy_local

    pure real(dp) function self_scale_local(xi)
        real(dp), intent(in) :: xi

        self_scale_local = -2.0_dp * xi / sqrt(pi_local())
    end function self_scale_local

    pure real(dp) function pi_local()
        pi_local = acos(-1.0_dp)
    end function pi_local

    pure real(dp) function combined_error(e_a, phi_a, e_b, phi_b)
        real(dp), intent(in) :: e_a(3)
        real(dp), intent(in) :: phi_a
        real(dp), intent(in) :: e_b(3)
        real(dp), intent(in) :: phi_b

        combined_error = maxval(abs(e_a - e_b)) + abs(phi_a - phi_b)
    end function combined_error
end program test_ewald2p_reciprocal
