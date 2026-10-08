! Copyright (c) 2026, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!
!> Implements `duprat_t`.
module duprat
  use field, only : field_t
  use num_types, only : rp
  use json_module, only : json_file
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use wall_model, only : wall_model_t
  use wall_sampler, only : wall_sampler_t
  use wall_sampler_fctry, only : wall_sampler_factory
  use user_intf, only : user_t
  use registry, only : neko_registry
  use json_utils, only : json_get_or_lookup, json_get_or_lookup_or_default, &
       json_get_or_default, json_get_subdict_or_empty
  use duprat_cpu, only : duprat_compute_cpu
  use field_math, only : field_invcol3
  use vector, only : vector_t
  use math, only : masked_gather_copy_0
  use scratch_registry, only : neko_scratch_registry
  use logger, only : LOG_SIZE, neko_log
  use utils, only : neko_error

  implicit none
  private

  !> Wall model of Duprat et al. (2011), which accounts for the streamwise
  !! pressure gradient through the extended velocity scale
  !! \f$ u_{\tau p} = \sqrt{u_\tau^2 + u_p^2} \f$, with the pressure velocity
  !! \f$ u_p = |(\nu / \rho) \partial p / \partial s|^{1/3} \f$, and a
  !! pressure-gradient-dependent mixing-length eddy viscosity. The velocity
  !! profile is obtained by integrating the simplified thin-boundary-layer
  !! equation across the wall layer.
  !! Reference: https://doi.org/10.1063/1.3529358
  !! @note Currently only a constant, prescribed pressure gradient is
  !! supported, and only on the CPU backend.
  type, public, extends(wall_model_t) :: duprat_t
     !> The von Karman coefficient.
     real(kind=rp) :: kappa = 0.41_rp
     !> The exponent of the pressure-gradient term in the eddy viscosity.
     real(kind=rp) :: beta = 0.78_rp
     !> The damping constant of the eddy viscosity.
     real(kind=rp) :: A = 17.0_rp
     !> The prescribed wall-tangential pressure gradient, positive when
     !! adverse to the local flow direction.
     real(kind=rp) :: dpds_constant = 0.0_rp
     !> The kinematic viscosity.
     type(vector_t) :: nu
     !> The fluid density at the boundary.
     type(vector_t) :: rho_w
     !> The wall-tangential pressure gradient at the wall nodes.
     type(vector_t) :: dpds
     !> Velocity sampled away from the wall.
     type(vector_t) :: u_s, v_s, w_s
   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => duprat_init
     !> Partial constructor from JSON, meant to work as the first stage of
     !! initialization before the `finalize` call.
     procedure, pass(this) :: partial_init => duprat_partial_init
     !> Finalize the construction using the mask and facet arrays of the bc.
     procedure, pass(this) :: finalize => duprat_finalize
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          duprat_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => duprat_free
     !> Compute the kinematic viscosity at the wall.
     procedure, pass(this) :: compute_nu => duprat_compute_nu
     !> Compute the wall shear stress.
     procedure, pass(this) :: compute => duprat_compute
  end type duprat_t

contains
  !> Constructor from JSON.
  !! @param scheme_name The name of the scheme for which the wall model is used.
  !! @param coef SEM coefficients.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param json A dictionary with parameters.
  subroutine duprat_init(this, scheme_name, coef, msk, facet, json)
    class(duprat_t), intent(inout) :: this
    character(len=*), intent(in) :: scheme_name
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    type(json_file), intent(inout) :: json
    real(kind=rp) :: kappa, beta, A, dpds_constant
    class(wall_sampler_t), allocatable :: sampler

    call duprat_read_json(json, kappa, beta, A, dpds_constant)

    call wall_sampler_factory(sampler, json)
    call this%init_from_components(scheme_name, coef, msk, facet, sampler, &
         kappa, beta, A, dpds_constant)
  end subroutine duprat_init

  !> Constructor from JSON.
  !! @param coef SEM coefficients.
  !! @param scheme_name The name of the scheme for which the wall model is used.
  !! @param json A dictionary with parameters.
  subroutine duprat_partial_init(this, coef, scheme_name, json)
    class(duprat_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    character(len=*), intent(in) :: scheme_name
    type(json_file), intent(inout) :: json
    character(len=LOG_SIZE) :: log_buf

    call duprat_check_backend()

    call this%partial_init_base(coef, scheme_name, json)
    call duprat_read_json(json, this%kappa, this%beta, this%A, &
         this%dpds_constant)

    call neko_log%section('Wall model')
    write(log_buf, '(A)') 'Model : Duprat'
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'kappa : ', this%kappa
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'beta : ', this%beta
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'A : ', this%A
    call neko_log%message(log_buf)
    write(log_buf, '(A, E15.7)') 'dp/ds (constant) : ', this%dpds_constant
    call neko_log%message(log_buf)
    call neko_log%end_section()

  end subroutine duprat_partial_init

  !> Read the model parameters from JSON.
  !! @param json A dictionary with parameters.
  !! @param kappa The von Karman coefficient.
  !! @param beta The exponent of the pressure-gradient term.
  !! @param A The damping constant.
  !! @param dpds_constant The prescribed wall-tangential pressure gradient.
  subroutine duprat_read_json(json, kappa, beta, A, dpds_constant)
    type(json_file), intent(inout) :: json
    real(kind=rp), intent(out) :: kappa, beta, A, dpds_constant
    type(json_file) :: pg_json
    character(len=:), allocatable :: pg_type

    call json_get_or_lookup(json, "kappa", kappa)
    call json_get_or_lookup_or_default(json, "beta", beta, 0.78_rp)
    call json_get_or_lookup_or_default(json, "A", A, 17.0_rp)

    ! The pressure gradient settings. Without them, the model runs with a
    ! zero pressure gradient, i.e. in its equilibrium limit.
    call json_get_subdict_or_empty(json, "pressure_gradient", pg_json)
    call json_get_or_default(pg_json, "type", pg_type, "constant")
    select case (trim(pg_type))
    case ("constant")
       call json_get_or_lookup_or_default(pg_json, "value", dpds_constant, &
            0.0_rp)
    case default
       call neko_error("Unknown pressure_gradient type '" // trim(pg_type) &
            // "' for the duprat wall model. Supported: constant")
    end select
    call pg_json%destroy()

  end subroutine duprat_read_json

  !> Finalize the construction using the mask and facet arrays of the bc.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param bc_name The name of the boundary condition.
  !! @param user The user interface.
  subroutine duprat_finalize(this, msk, facet, bc_name, user)
    class(duprat_t), intent(inout) :: this
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    character(len=*), optional, intent(in) :: bc_name
    type(user_t), target, optional, intent(in) :: user

    call this%finalize_base(msk, facet, bc_name, user)
    call duprat_init_vectors(this)
  end subroutine duprat_finalize

  !> Constructor from components.
  !! @param scheme_name The name of the scheme for which the wall model is used.
  !! @param coef SEM coefficients.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param sampler The sampling strategy. Ownership is transferred.
  !! @param kappa The von Karman coefficient.
  !! @param beta The exponent of the pressure-gradient term.
  !! @param A The damping constant.
  !! @param dpds_constant The prescribed wall-tangential pressure gradient,
  !! positive when adverse to the local flow direction.
  subroutine duprat_init_from_components(this, scheme_name, coef, msk, &
       facet, sampler, kappa, beta, A, dpds_constant)
    class(duprat_t), intent(inout) :: this
    character(len=*), intent(in) :: scheme_name
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    class(wall_sampler_t), allocatable, intent(inout) :: sampler
    real(kind=rp), intent(in) :: kappa, beta, A, dpds_constant

    call duprat_check_backend()

    call this%free()
    call this%init_base(scheme_name, coef, msk, facet, sampler)

    this%kappa = kappa
    this%beta = beta
    this%A = A
    this%dpds_constant = dpds_constant

    call duprat_init_vectors(this)
  end subroutine duprat_init_from_components

  !> Allocate the work vectors on the wall nodes.
  subroutine duprat_init_vectors(this)
    class(duprat_t), intent(inout) :: this

    call this%nu%init(this%n_nodes)
    call this%rho_w%init(this%n_nodes)
    call this%dpds%init(this%n_nodes)
    this%dpds%x = this%dpds_constant
    call this%validate_single_sample()
    call this%u_s%init(this%n_nodes)
    call this%v_s%init(this%n_nodes)
    call this%w_s%init(this%n_nodes)
  end subroutine duprat_init_vectors

  !> Stop with an error if a device backend is used, since the Duprat
  !! wall model is so far only implemented for the CPU backend.
  subroutine duprat_check_backend()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error("The duprat wall model is only implemented " // &
            "for the CPU backend")
    end if
  end subroutine duprat_check_backend

  !> Compute the kinematic viscosity vector.
  subroutine duprat_compute_nu(this)
    class(duprat_t), intent(inout) :: this
    type(field_t), pointer :: temp
    integer :: idx

    call neko_scratch_registry%request_field(temp, idx, .false.)
    call field_invcol3(temp, this%mu, this%rho)

    call masked_gather_copy_0(this%nu%x, temp%x, this%msk, temp%size(), &
         this%nu%size())
    call masked_gather_copy_0(this%rho_w%x, this%rho%x, this%msk, &
         this%rho%size(), this%rho_w%size())

    call neko_scratch_registry%relinquish_field(idx)
  end subroutine duprat_compute_nu

  !> Destructor for the duprat_t class.
  subroutine duprat_free(this)
    class(duprat_t), intent(inout) :: this

    call this%nu%free()
    call this%rho_w%free()
    call this%dpds%free()
    call this%u_s%free()
    call this%v_s%free()
    call this%w_s%free()
    call this%free_base()

  end subroutine duprat_free

  !> Compute the wall shear stress.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine duprat_compute(this, t, tstep)
    class(duprat_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), pointer :: u
    type(field_t), pointer :: v
    type(field_t), pointer :: w

    call this%compute_nu()

    u => neko_registry%get_field("u")
    v => neko_registry%get_field("v")
    w => neko_registry%get_field("w")

    call this%sampler%sample(u, this%u_s)
    call this%sampler%sample(v, this%v_s)
    call this%sampler%sample(w, this%w_s)

    call duprat_compute_cpu(this%u_s%x, this%v_s%x, this%w_s%x, &
         this%n_x%x, this%n_y%x, this%n_z%x, &
         this%nu%x, this%rho_w%x, this%sampler%h%x, this%dpds%x, &
         this%tau_x%x, this%tau_y%x, this%tau_z%x, &
         this%n_nodes, this%kappa, this%beta, this%A, tstep)

    nullify(u, v, w)

  end subroutine duprat_compute
end module duprat
