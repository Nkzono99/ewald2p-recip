module ewald2p_recip_zero_mode
    use ewald2p_recip_types, only: dp
    use ewald2p_recip_utils, only: pi_dp
    implicit none
    private

    public :: eval_zero_mode_contribution
    public :: self_potential_scale

contains

    pure subroutine eval_zero_mode_contribution(lx, ly, xi, src_z, src_q, z, phi0, ez0)
        real(dp), intent(in) :: lx
        real(dp), intent(in) :: ly
        real(dp), intent(in) :: xi
        real(dp), intent(in) :: src_z(:)
        real(dp), intent(in) :: src_q(:)
        real(dp), intent(in) :: z
        real(dp), intent(out) :: phi0
        real(dp), intent(out) :: ez0

        integer :: isrc
        real(dp) :: area
        real(dp) :: dz
        real(dp) :: xi_dz
        real(dp) :: gaussian
        real(dp) :: erf_term
        real(dp) :: sqrt_pi
        real(dp) :: phi_prefactor
        real(dp) :: ez_prefactor

        area = lx * ly
        sqrt_pi = sqrt(pi_dp())
        phi_prefactor = -2.0_dp * sqrt_pi / area
        ez_prefactor = 2.0_dp * pi_dp() / area

        phi0 = 0.0_dp
        ez0 = 0.0_dp

        do isrc = 1, size(src_q)
            dz = z - src_z(isrc)
            xi_dz = xi * dz
            gaussian = exp(-(xi_dz * xi_dz))
            erf_term = erf(xi_dz)

            phi0 = phi0 + phi_prefactor * src_q(isrc) * ((gaussian / xi) + sqrt_pi * dz * erf_term)
            ez0 = ez0 + ez_prefactor * src_q(isrc) * erf_term
        end do
    end subroutine eval_zero_mode_contribution

    pure real(dp) function self_potential_scale(xi)
        real(dp), intent(in) :: xi

        self_potential_scale = -2.0_dp * xi / sqrt(pi_dp())
    end function self_potential_scale
end module ewald2p_recip_zero_mode
