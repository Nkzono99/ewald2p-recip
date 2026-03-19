module ewald2p_reciprocal
    use ewald2p_recip_types, only: dp, ip, ewald2p_recip_options_type, ewald2p_recip_plan_type, &
                                   ewald2p_recip_state_type
    use ewald2p_reciprocal_build, only: build_recip_plan, destroy_recip_plan
    use ewald2p_reciprocal_eval, only: eval_recip_point, eval_recip_points
    use ewald2p_reciprocal_state, only: check_neutrality, destroy_recip_state, update_recip_state
    implicit none
    private

    public :: dp
    public :: ip
    public :: ewald2p_recip_options_type
    public :: ewald2p_recip_plan_type
    public :: ewald2p_recip_state_type
    public :: build_recip_plan
    public :: destroy_recip_plan
    public :: update_recip_state
    public :: destroy_recip_state
    public :: eval_recip_point
    public :: eval_recip_points
    public :: check_neutrality
end module ewald2p_reciprocal
