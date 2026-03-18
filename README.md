# ewald2p-recip

Fortran library for the reciprocal-space part of the neutral 2-periodic Ewald sum for Coulomb point charges.

## Overview

`ewald2p-recip` is a Fortran library that evaluates the **reciprocal-space contribution** of the **neutral 2-periodic Ewald sum** for **point charges** with the Coulomb kernel `1/r`.

The intended use is to combine this library with a separate **real-space solver**, such as:

- Cartesian / Taylor FMM
- direct real-space Ewald sum
- other short-range solvers

An outer adapter is expected to combine:

- real-space contribution
- reciprocal-space contribution from this library
- zero-mode contribution
- self correction, when applicable

## Scope of v0.1

This repository currently targets the following model:

- Coulomb kernel `1/r`
- point sources only
- periodic in `x` and `y`
- open in `z`
- neutral systems only
- standard neutral 2P Ewald `k=0` term included
- reciprocal-space evaluated by direct mode summation

## Non-goals for v0.1

The following are intentionally out of scope for the initial version:

- non-neutral systems
- charged-wall regularization
- uniform neutralizing background
- softened kernels
- panel form factors for BEM elements
- FFT / PME / NFFT / NUFFT acceleration
- real-space Ewald evaluation
- energy / virial APIs

## Design philosophy

This library is designed to remain **independent** from:

- BEM mesh types
- simulator configuration types
- FMM internal types

The public API is array-based and works with:

- `src_pos(3, nsrc)`
- `src_q(nsrc)`
- `target_pos(3, ntgt)`

This makes the library suitable as a reusable internal module behind a solver adapter.

## Planned public API

The intended public interface is:

```fortran
call build_recip_plan(plan, options)
call update_recip_state(plan, state, src_pos, src_q)
call eval_recip_point(plan, state, r, e, phi)
call eval_recip_points(plan, state, target_pos, e, phi)
````

## Repository layout

```text
.
├─ README.md
├─ SPEC.md
├─ AGENTS.md
├─ fpm.toml
├─ src/
│  └─ ewald2p_recip.f90
└─ test/
```

As the implementation grows, the source tree is expected to be split into smaller modules such as:

```text
src/
  ewald2p_recip.f90
  ewald2p_recip_build.f90
  ewald2p_recip_state.f90
  ewald2p_recip_eval.f90
  internal/
    ewald2p_recip_types.f90
    ewald2p_recip_kspace.f90
    ewald2p_recip_zero_mode.f90
    ewald2p_recip_wrap.f90
    ewald2p_recip_utils.f90
```

## Build

This project uses [fpm](https://fpm.fortran-lang.org/).

```bash
fpm build
```

Run tests with:

```bash
fpm test
```

## Status

This repository is currently specification-first.

Current focus:

1. establish the reciprocal-space specification
2. implement a direct mode-sum reference backend
3. validate neutrality, zero-mode, and self-interaction handling
4. add tests before optimization

## Notes

* This library evaluates the reciprocal-space part only.
* The real-space part must be supplied separately.
* Non-neutral systems are rejected by design in the initial implementation.

## License

MIT