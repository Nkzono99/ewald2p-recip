module ewald2p_recip_kspace
    use ewald2p_recip_types, only: dp, ip
    use ewald2p_recip_utils, only: ERR_MEMORY_ALLOCATION, clear_error, pi_dp, set_error
    implicit none
    private

    public :: enumerate_nonzero_modes
    public :: nonzero_mode_kernel

contains

    subroutine enumerate_nonzero_modes(lx, ly, h_max, l_max, nk, h_list, l_list, kx, ky, kabs, ierr, message)
        real(dp), intent(in) :: lx
        real(dp), intent(in) :: ly
        integer(ip), intent(in) :: h_max
        integer(ip), intent(in) :: l_max
        integer(ip), intent(out) :: nk
        integer(ip), allocatable, intent(out) :: h_list(:)
        integer(ip), allocatable, intent(out) :: l_list(:)
        real(dp), allocatable, intent(out) :: kx(:)
        real(dp), allocatable, intent(out) :: ky(:)
        real(dp), allocatable, intent(out) :: kabs(:)
        integer(ip), intent(out) :: ierr
        character(len=*), intent(out), optional :: message

        integer(ip) :: expected_nk
        integer(ip) :: h
        integer(ip) :: l
        integer(ip) :: idx
        integer :: stat
        real(dp) :: two_pi

        call clear_error(ierr, message)

        expected_nk = (2_ip * h_max + 1_ip) * (2_ip * l_max + 1_ip) - 1_ip
        nk = expected_nk

        allocate (h_list(nk), l_list(nk), kx(nk), ky(nk), kabs(nk), stat=stat)
        if (stat /= 0) then
            call set_error(ierr, ERR_MEMORY_ALLOCATION, "Failed to allocate k-space tables.", message)
            return
        end if

        two_pi = 2.0_dp * pi_dp()
        idx = 0_ip

        do l = -l_max, l_max
            do h = -h_max, h_max
                if (h == 0_ip .and. l == 0_ip) cycle

                idx = idx + 1_ip
                h_list(idx) = h
                l_list(idx) = l
                kx(idx) = two_pi * real(h, dp) / lx
                ky(idx) = two_pi * real(l, dp) / ly
                kabs(idx) = sqrt(kx(idx) * kx(idx) + ky(idx) * ky(idx))
            end do
        end do
    end subroutine enumerate_nonzero_modes

    pure subroutine nonzero_mode_kernel(xi, kabs, dz, kernel_sum, kernel_diff)
        real(dp), intent(in) :: xi
        real(dp), intent(in) :: kabs
        real(dp), intent(in) :: dz
        real(dp), intent(out) :: kernel_sum
        real(dp), intent(out) :: kernel_diff

        real(dp) :: a
        real(dp) :: xi_dz
        real(dp) :: term_plus
        real(dp) :: term_minus

        a = kabs / (2.0_dp * xi)
        xi_dz = xi * dz

        term_plus = exp(kabs * dz) * erfc(a + xi_dz)
        term_minus = exp(-kabs * dz) * erfc(a - xi_dz)

        kernel_sum = term_plus + term_minus
        kernel_diff = term_plus - term_minus
    end subroutine nonzero_mode_kernel
end module ewald2p_recip_kspace
