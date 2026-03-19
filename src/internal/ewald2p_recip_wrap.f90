module ewald2p_recip_wrap
    use ewald2p_recip_types, only: dp
    implicit none
    private

    public :: wrap_periodic
    public :: wrap_xy_point

contains

    pure real(dp) function wrap_periodic(value, period)
        real(dp), intent(in) :: value
        real(dp), intent(in) :: period

        if (period <= 0.0_dp) then
            wrap_periodic = value
            return
        end if

        wrap_periodic = modulo(value, period)

        if (wrap_periodic >= period) wrap_periodic = wrap_periodic - period
        if (wrap_periodic < 0.0_dp) wrap_periodic = wrap_periodic + period
    end function wrap_periodic

    pure subroutine wrap_xy_point(lx, ly, r_in, r_out)
        real(dp), intent(in) :: lx
        real(dp), intent(in) :: ly
        real(dp), intent(in) :: r_in(3)
        real(dp), intent(out) :: r_out(3)

        r_out = r_in
        r_out(1) = wrap_periodic(r_in(1), lx)
        r_out(2) = wrap_periodic(r_in(2), ly)
    end subroutine wrap_xy_point
end module ewald2p_recip_wrap
