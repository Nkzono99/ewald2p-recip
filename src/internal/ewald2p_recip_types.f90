module ewald2p_recip_types
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none
    private

    public :: dp, ip
    public :: ewald2p_recip_options_type
    public :: ewald2p_recip_plan_type
    public :: ewald2p_recip_state_type

    integer, parameter :: dp = real64
    integer, parameter :: ip = int32

    type :: ewald2p_recip_options_type
        real(dp) :: lx = 0.0_dp
        real(dp) :: ly = 0.0_dp
        real(dp) :: xi = 0.0_dp
        integer(ip) :: h_max = 0_ip
        integer(ip) :: l_max = 0_ip
        real(dp) :: neutrality_tol = 1.0e-12_dp
        logical :: return_potential = .false.
        logical :: enable_zero_mode = .true.
    end type ewald2p_recip_options_type

    type :: ewald2p_recip_plan_type
        logical :: is_built = .false.

        real(dp) :: lx = 0.0_dp
        real(dp) :: ly = 0.0_dp
        real(dp) :: xi = 0.0_dp
        integer(ip) :: h_max = 0_ip
        integer(ip) :: l_max = 0_ip
        real(dp) :: neutrality_tol = 0.0_dp

        integer(ip) :: nk = 0_ip
        real(dp), allocatable :: kx(:)
        real(dp), allocatable :: ky(:)
        real(dp), allocatable :: kabs(:)

        integer(ip), allocatable :: h_list(:)
        integer(ip), allocatable :: l_list(:)
    end type ewald2p_recip_plan_type

    type :: ewald2p_recip_state_type
        logical :: is_ready = .false.

        integer(ip) :: nsrc = 0_ip

        real(dp), allocatable :: src_pos(:, :)
        real(dp), allocatable :: src_q(:)

        ! Cached plane-wave moments. The reference evaluator retains the full
        ! source set because the 2P z-kernel is not separable in target z.
        complex(dp), allocatable :: s_plus(:)
        complex(dp), allocatable :: s_minus(:)
    end type ewald2p_recip_state_type
end module ewald2p_recip_types
