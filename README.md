# xDDPhaseSpace.jl

A Julia package for calculating 3-body phase space with analytic continuation into the complex plane. xDDPhaseSpace.jl calculates 3-body phase space integrals with complex energy support. It handles analytic continuation by carefully selecting integration paths to avoid singularities.

## Covariant model description

- Study of the doubly charmed tetraquark $T_{cc}^+$​ - [INSPIRE:1915358](https://inspirehep.net/literature/1915358)
- Effective-range expansion of the $T_{cc}^+$ - [INSPIRE:2048996](https://inspirehep.net/literature/2048996)

## Architecture

### Abstract xDD

Is an abstract type for all xDD models.
To extend the package, one needs to implement the following methods:
  - `decay_matrix_element_squared` - Decay matrix elements
  - `integrand_mapped_thr` - Mapped integrand for threshold
  - `mapdalitzmethod` - Dalitz plot mapping
  - `ρ_thr` - Threshold phase space density
  - `masses` - Particle mass accessor


### Available models

#### Symmetrized VP S-wave

Decay amplitude structure:
```
A = A12 + A13
A12 = R(12) P-wave, +3 in S-wave   
A13 = R(13) P-wave, +2 in S-wave   
```

Examples:
- Tcc → D*+ D0 (symmetrized)
- χc1 → D*0 D̄0 (c.c.)

Two main decay models are implemented:
1. `πDD` - Hadronic decay with pion emission
2. `γDD` - Radiative decay with photon emission  


### Complex Plane Handling

The `integrationpaths.jl` defines the core integration components:
- Integration contours in complex plane
- Path deformation around poles and branch points
- Dalitz plot mapping implementations

## Installation

```julia
using Pkg
Pkg.add("xDDPhaseSpace")
```

## Usage

### Example: πDD and γDD decay channels with D*+ resonances

Both `πDD` and `γDD` models follow similar API calls, but call different decay amplitudes.
The radiative transition needs M1 transition matrix elements for both topologies in addition to masses and the resonant parametrizations.

```julia
using xDDPhaseSpace

# Define particle masses (GeV)
const mπ⁺ = 0.13957039
const mD⁰ = 1.86483
const mDˣ⁺ = 2.0102558
const ΓDˣ⁺ = 0.00083  # D*+ decay width

# πDD: Pion decay to D meson pairs
decay_channel_pi = πDD(
    (m1 = mπ⁺, m2 = mD⁰, m3 = mD⁰), 
    BW(m = mDˣ⁺, Γ = ΓDˣ⁺), 
    BW(m = mDˣ⁺, Γ = ΓDˣ⁺)
)

# γDD: Photon decay to D meson pairs with electromagnetic couplings
μ₀ = -3.77  # Transition matrix element < D0 | mu | D*0 >
decay_channel_gamma = γDD(
    (m1 = 0.0, m2 = mD⁰, m3 = mD⁰),  # m1 = 0.0 for photon
    BW(m = mDˣ⁺, Γ = ΓDˣ⁺), 
    BW(m = mDˣ⁺, Γ = ΓDˣ⁺),
    μ₀, -μ₀
)

# Both models work with the same phase space infrastructure
ρ_value_pi = ρ_thr(decay_channel_pi, WithThrE(1.0, mD⁰ + mDˣ⁺))
ρ_value_gamma = ρ_thr(decay_channel_gamma, WithThrE(1.0, mD⁰ + mDˣ⁺))

# Complex energies are supported for both
ρ_value_pi_c = ρ_thr(decay_channel_pi, WithThrE(1.0 + 0.1im, mD⁰ + mDˣ⁺))
ρ_value_gamma_c = ρ_thr(decay_channel_gamma, WithThrE(1.0 + 0.1im, mD⁰ + mDˣ⁺))

# Get complex branch points
branch_pts_pi = branch_points(decay_channel_pi)
branch_pts_gamma = branch_points(decay_channel_gamma)
```

### Interpolated threshold density

For channels that are expensive to evaluate repeatedly near threshold, you can
tabulate `ρ_thr` on a fixed mass grid and then reuse the interpolated
representation. The package keeps the exact tabulated values below the chosen
cutoff and matches smoothly to a rescaled `ρ_tb` above it.

```julia
using xDDPhaseSpace

const mγ = 0.0
const mD0 = 1.86483
const mDstar0 = 2.00685
const ΓDstar0 = 55.2e-6
const μ0 = -3.77

BW0 = BW(m = mDstar0, Γ = ΓDstar0)

ch_gamma = γDD(
    (m1 = mγ, m2 = mD0, m3 = mD0),
    BW0,
    BW0,
    μ0, -μ0,
)

# Build an interpolated version of ρ_thr up to the selected cutoff.
ch_gamma_interpolated = interpolated(ch_gamma, 3.8725; estep_GeV = 1.3e-3)

thr = mDstar0 + mD0
m = thr + 5 * 1.3e-3

# In the interpolation region, ρ_thr is read from the tabulated grid.
ρ_direct = ρ_thr(ch_gamma, m)
ρ_interp = ρ_thr(ch_gamma_interpolated, m)

# The dispersive construction uses the interpolated density as input.
disp = dispersive(ch_gamma_interpolated, m)

ρ_direct
ρ_interp
imag(disp)
```

At grid points used to build the interpolant, `ρ_thr(ch_gamma_interpolated, m)`
matches the tabulated `ρ_thr(ch_gamma, m)`. For real `m`, the numerical
imaginary part of `dispersive(ch_gamma_interpolated, m)` reproduces
`ρ_thr(ch_gamma_interpolated, m)`.
