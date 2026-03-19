
# AGENTS.md

やりとりはすべて日本語で行うこと.

## Purpose

This repository contains a **Fortran implementation** of the reciprocal-space part of the **two-periodic Ewald sum** for **point charges** with the Coulomb kernel \(1/r\).

The code is intended to be used as an **independent internal library**.
It must not depend directly on BEM mesh types, simulator configuration types, or FMM-specific internal plan/state types from other subsystems.

The main use case is:

- real-space part: evaluated by a separate solver, e.g. Cartesian / Taylor FMM
- reciprocal-space part: evaluated by this library
- outer adapter: combines both contributions

---

## Language and coding rules

- Primary language: **Fortran**
- Prefer **Fortran 2008 compatible** style unless there is a strong reason otherwise
- Use `implicit none` everywhere
- Use explicit `public` / `private`
- Keep module interfaces stable and minimal
- Avoid hidden global state
- Avoid side effects outside plan/state objects
- Use clear and boring names over clever names

Do not introduce Python, C, C++, or build-system-specific logic unless explicitly requested.

---

## Scope boundaries

### This library owns

- reciprocal-space evaluation for **2-periodic slab Ewald**
- nonzero Fourier mode handling
- analytic `k=0` contribution
- neutrality checking
- self-interaction policy for reciprocal-side potential evaluation
- source-dependent reciprocal state
- source-independent reciprocal plan

### This library does not own

- real-space Ewald evaluation
- FMM implementation
- BEM panel integration
- mesh element definitions
- simulator configuration
- Coulomb constant scaling in external adapters unless explicitly specified
- non-neutral regularization beyond explicit future extensions

If a requested change crosses these boundaries, keep the reciprocal module independent and introduce an adapter rather than coupling directly to external solver internals.

---

## Physical model assumptions

Current supported model:

- kernel: Coulomb \(1/r\)
- source model: point charges
- boundary condition: periodic in `x` and `y`, open in `z`
- only **neutral systems** are supported
- standard neutral 2P Ewald `k=0` term is included
- self interaction is excluded

Do not silently add support for:

- non-neutral systems
- charged walls
- uniform neutralizing background
- softened kernels
- panel form factors
- FFT/NFFT acceleration

These belong to future extensions and must remain opt-in.

---

## Required design philosophy

Keep the implementation split into three layers.

### 1. Public API layer

Small wrapper modules exposing procedures like:

- `build_recip_plan`
- `destroy_recip_plan`
- `update_recip_state`
- `destroy_recip_state`
- `eval_recip_point`
- `eval_recip_points`
- `check_neutrality`

This layer should be thin.

### 2. Internal algorithm layer

Contains:

- k-space mode enumeration
- source coefficient accumulation
- reciprocal evaluation kernels
- zero-mode evaluation
- coordinate wrapping helpers
- validation utilities

### 3. External adapter layer

Not part of this library.

External code may convert:

- BEM element centroids -> `src_pos`
- element charges -> `src_q`

This library must not know where those arrays came from.

---

## File structure guidance

Recommended structure:

```text
src/
  physics/
    field_solver/
      ewald/
        ewald2p_reciprocal.f90
        ewald2p_reciprocal_build.f90
        ewald2p_reciprocal_state.f90
        ewald2p_reciprocal_eval.f90
        internal/
          ewald2p_recip_types.f90
          ewald2p_recip_kspace.f90
          ewald2p_recip_zero_mode.f90
          ewald2p_recip_wrap.f90
          ewald2p_recip_utils.f90
````

Keep source-independent and source-dependent logic separate.

---

## Data ownership rules

### Plan objects

Plan objects contain only source-independent information, such as:

* box lengths
* Ewald parameter
* mode index tables
* wave numbers
* precomputed constants

A plan must remain reusable across charge/state updates.

### State objects

State objects contain source-dependent information, such as:

* source positions if retained
* source charges if retained
* reciprocal mode coefficients

A state becomes invalid if the source set changes and must then be rebuilt.

---

## Neutrality policy

The library currently supports only:

```text
neutrality_mode = require_neutral
```

Meaning:

* compute total charge
* compare against tolerance
* reject non-neutral systems

Never silently continue on a non-neutral input.

Never auto-add background charge or charged walls.

If future support for non-neutral systems is added, it must be introduced through an explicit new option and documented as a different physical model.

---

## Zero-mode policy

The library currently supports only:

```text
zero_mode_policy = standard_2p_neutral
```

Meaning:

* the `k=0` contribution is part of the physical model
* it must be evaluated explicitly
* it must not be silently disabled

Do not add a runtime option that skips zero mode unless the specification is also extended and the physical consequences are documented.

---

## Self-interaction policy

The library currently supports only:

```text
self_policy = exclude
```

Meaning:

* self interaction is excluded
* for potential evaluation at source locations, apply the standard reciprocal-side self correction when appropriate
* do not invent vector self corrections for electric field unless explicitly derived and specified

Keep this policy explicit in code and tests.

---

## Numerical implementation policy

### First implementation target

Use **direct mode summation** for reciprocal-space evaluation.

This means:

* no FFT backend in v1
* no NFFT/NUFFT backend in v1
* correctness before acceleration

### Why

The direct mode-sum implementation is easier to verify and is a better reference backend for future fast implementations.

If a fast backend is introduced later, it must reproduce the direct backend within documented tolerances.

---

## Error handling

Public procedures should return:

* integer error code
* optional message string

Do not stop the program from deep internal routines unless the repository already has a consistent fatal-error convention and the caller explicitly expects it.

Prefer returning errors over hidden failure.

Typical failures include:

* invalid box lengths
* invalid Ewald parameter
* invalid mode cutoff
* unbuilt plan
* uninitialized state
* size mismatch
* non-neutral system
* unsupported option combination

---

## Testing expectations

Every feature addition should include tests when feasible.

Minimum tests to preserve:

* neutrality rejection
* periodicity under `x -> x + Lx`
* periodicity under `y -> y + Ly`
* reciprocal convergence with increasing mode cutoffs
* consistency between pointwise and batched evaluation
* correct handling of self policy
* correct inclusion of zero mode

When adding optimization, compare against the direct mode-sum reference implementation.

---

## Performance rules

* First make it correct
* Then make it measurable
* Then optimize

Avoid premature optimization that obscures the physics or numerical meaning.

Safe optimizations include:

* precomputing mode tables
* precomputing repeated constants
* batched evaluation
* OpenMP over targets or modes where race conditions are controlled

Unsafe optimizations include:

* undocumented approximations
* silently dropping terms
* changing the physical model to gain speed

---

## Documentation rules

When changing formulas, also update:

* `SPEC.md`
* public comments near the API
* any derivation notes if present

Do not let implementation drift away from the written specification.

If the code and `SPEC.md` disagree, align them immediately.

---

## Style rules for contributors

When editing or generating code:

* preserve existing naming conventions if reasonable
* keep procedures short and single-purpose
* separate validation, math kernels, and orchestration logic
* add comments for non-obvious formulas, not for trivial assignments
* prefer readable loops over overly compressed logic
* avoid introducing hidden assumptions
* ignore `*.i90` files, as they are automatically generated backup files.

When adding a new option, document:

* default
* allowed values
* physical meaning
* numerical consequences

---

## Change management guidance

### Acceptable changes without redesign

* adding helper routines
* improving comments
* refactoring internal loops
* adding tests
* adding OpenMP pragmas in independent loops
* improving batched evaluation

### Changes that require spec update first

* non-neutral support
* softened kernel support
* FFT/NFFT backends
* panel form factors
* energy or virial outputs
* altered self policy
* altered zero-mode handling

If a change modifies the physical meaning of the solver, update `SPEC.md` before or together with the implementation.

---

## What to avoid

Do not:

* couple directly to BEM element types
* couple directly to simulator config types
* mix reciprocal and real-space logic in one large procedure
* hide zero-mode decisions
* silently ignore non-neutral inputs
* silently switch kernels
* mix unit scaling responsibilities without documenting it

---

## Preferred development order

When building from scratch, the recommended order is:

1. define types
2. implement plan builder
3. implement neutrality check
4. implement state update with mode coefficients
5. implement direct reciprocal point evaluation
6. implement batched evaluation
7. add zero-mode routine
8. add self handling for potential evaluation
9. add tests
10. add OpenMP only after validation

---

## Final principle

This library should remain a **clean, reusable, Fortran reciprocal-space engine** for the neutral 2-periodic Coulomb Ewald problem with point sources.

Keep it physically explicit, numerically verifiable, and weakly coupled to the rest of the codebase.
