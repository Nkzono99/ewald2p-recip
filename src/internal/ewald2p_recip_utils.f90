module ewald2p_recip_utils
    use ewald2p_recip_types, only: dp, ip
    implicit none
    private

    public :: ERR_SUCCESS
    public :: ERR_INVALID_OPTION
    public :: ERR_PLAN_NOT_BUILT
    public :: ERR_STATE_NOT_READY
    public :: ERR_SIZE_MISMATCH
    public :: ERR_NON_NEUTRAL_SYSTEM
    public :: ERR_ZERO_MODE_UNSUPPORTED
    public :: ERR_MEMORY_ALLOCATION
    public :: ERR_INTERNAL_NUMERICAL
    public :: clear_error
    public :: set_error
    public :: charge_scale
    public :: total_charge_value
    public :: coincidence_tol
    public :: pi_dp

    integer(ip), parameter :: ERR_SUCCESS = 0_ip
    integer(ip), parameter :: ERR_INVALID_OPTION = 1_ip
    integer(ip), parameter :: ERR_PLAN_NOT_BUILT = 2_ip
    integer(ip), parameter :: ERR_STATE_NOT_READY = 3_ip
    integer(ip), parameter :: ERR_SIZE_MISMATCH = 4_ip
    integer(ip), parameter :: ERR_NON_NEUTRAL_SYSTEM = 5_ip
    integer(ip), parameter :: ERR_ZERO_MODE_UNSUPPORTED = 6_ip
    integer(ip), parameter :: ERR_MEMORY_ALLOCATION = 7_ip
    integer(ip), parameter :: ERR_INTERNAL_NUMERICAL = 8_ip

contains

    subroutine clear_error(ierr, message)
        integer(ip), intent(out) :: ierr
        character(len=*), intent(out), optional :: message

        ierr = ERR_SUCCESS
        if (present(message)) message = ""
    end subroutine clear_error

    subroutine set_error(ierr, code, text, message)
        integer(ip), intent(out) :: ierr
        integer(ip), intent(in) :: code
        character(len=*), intent(in) :: text
        character(len=*), intent(out), optional :: message

        ierr = code
        if (present(message)) message = text
    end subroutine set_error

    pure real(dp) function charge_scale(src_q)
        real(dp), intent(in) :: src_q(:)

        charge_scale = max(1.0_dp, sum(abs(src_q)))
    end function charge_scale

    pure real(dp) function total_charge_value(src_q)
        real(dp), intent(in) :: src_q(:)

        total_charge_value = sum(src_q)
    end function total_charge_value

    pure real(dp) function coincidence_tol(lx, ly, z_value)
        real(dp), intent(in) :: lx
        real(dp), intent(in) :: ly
        real(dp), intent(in) :: z_value
        real(dp) :: scale

        scale = max(1.0_dp, abs(z_value))
        scale = max(scale, lx)
        scale = max(scale, ly)

        coincidence_tol = 128.0_dp * epsilon(1.0_dp) * scale
    end function coincidence_tol

    pure real(dp) function pi_dp()
        pi_dp = acos(-1.0_dp)
    end function pi_dp
end module ewald2p_recip_utils
