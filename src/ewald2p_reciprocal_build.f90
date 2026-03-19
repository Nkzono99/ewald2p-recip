module ewald2p_reciprocal_build
    use ewald2p_recip_kspace, only: enumerate_nonzero_modes
    use ewald2p_recip_types, only: dp, ip, ewald2p_recip_options_type, ewald2p_recip_plan_type
    use ewald2p_recip_utils, only: ERR_INVALID_OPTION, ERR_ZERO_MODE_UNSUPPORTED, clear_error, set_error
    implicit none
    private

    public :: build_recip_plan
    public :: destroy_recip_plan

contains

    subroutine build_recip_plan(plan, options, ierr, message)
        type(ewald2p_recip_plan_type), intent(out) :: plan
        type(ewald2p_recip_options_type), intent(in) :: options
        integer(ip), intent(out) :: ierr
        character(len=*), intent(out), optional :: message

        call clear_error(ierr, message)
        call destroy_recip_plan(plan)

        if (options%lx <= 0.0_dp) then
            call set_error(ierr, ERR_INVALID_OPTION, "Invalid option: lx must be positive.", message)
            return
        end if

        if (options%ly <= 0.0_dp) then
            call set_error(ierr, ERR_INVALID_OPTION, "Invalid option: ly must be positive.", message)
            return
        end if

        if (options%xi <= 0.0_dp) then
            call set_error(ierr, ERR_INVALID_OPTION, "Invalid option: xi must be positive.", message)
            return
        end if

        if (options%h_max < 0_ip) then
            call set_error(ierr, ERR_INVALID_OPTION, "Invalid option: h_max must be nonnegative.", message)
            return
        end if

        if (options%l_max < 0_ip) then
            call set_error(ierr, ERR_INVALID_OPTION, "Invalid option: l_max must be nonnegative.", message)
            return
        end if

        if (options%neutrality_tol < 0.0_dp) then
            call set_error(ierr, ERR_INVALID_OPTION, "Invalid option: neutrality_tol must be nonnegative.", message)
            return
        end if

        if (.not. options%enable_zero_mode) then
            call set_error(ierr, ERR_ZERO_MODE_UNSUPPORTED, "enable_zero_mode=.false. is unsupported in v1.", message)
            return
        end if

        plan%lx = options%lx
        plan%ly = options%ly
        plan%xi = options%xi
        plan%h_max = options%h_max
        plan%l_max = options%l_max
        plan%neutrality_tol = options%neutrality_tol

        call enumerate_nonzero_modes(plan%lx, plan%ly, plan%h_max, plan%l_max, plan%nk, &
                                     plan%h_list, plan%l_list, plan%kx, plan%ky, plan%kabs, ierr, message)
        if (ierr /= 0_ip) then
            call destroy_recip_plan(plan)
            return
        end if

        plan%is_built = .true.
    end subroutine build_recip_plan

    subroutine destroy_recip_plan(plan)
        type(ewald2p_recip_plan_type), intent(inout) :: plan

        if (allocated(plan%kx)) deallocate (plan%kx)
        if (allocated(plan%ky)) deallocate (plan%ky)
        if (allocated(plan%kabs)) deallocate (plan%kabs)
        if (allocated(plan%h_list)) deallocate (plan%h_list)
        if (allocated(plan%l_list)) deallocate (plan%l_list)

        plan%is_built = .false.
        plan%lx = 0.0_dp
        plan%ly = 0.0_dp
        plan%xi = 0.0_dp
        plan%h_max = 0_ip
        plan%l_max = 0_ip
        plan%neutrality_tol = 0.0_dp
        plan%nk = 0_ip
    end subroutine destroy_recip_plan
end module ewald2p_reciprocal_build
