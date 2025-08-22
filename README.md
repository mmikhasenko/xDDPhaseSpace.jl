# xDDPhaseSpace.jl

A Julia package for calculating 3-body phase space with analytic continuation into the complex plane.

## Overview

xDDPhaseSpace.jl calculates 3-body phase space integrals with complex energy support. It handles analytic continuation by carefully selecting integration paths to avoid singularities.

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
- a1 → ρπ S (symmetrized)

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

### Example: πDD decay channel with D*+ resonances

```julia
using xDDPhaseSpace

# Define particle masses (GeV)
const mπ⁺ = 0.13957039
const mD⁰ = 1.86483
const mDˣ⁺ = 2.0102558
const ΓDˣ⁺ = 0.00083  # D*+ decay width

# Create decay channel with πDD and two D*+ Breit-Wigner resonances
decay_channel = πDD(
    (m1 = mπ⁺, m2 = mD⁰, m3 = mD⁰), 
    BW(m = mDˣ⁺, Γ = ΓDˣ⁺), 
    BW(m = mDˣ⁺, Γ = ΓDˣ⁺)
)

# Calculate phase space density 1 MeV above threshold
ρ_value = ρ_thr(decay_channel, WithThrE(1.0, mD⁰ + mDˣ⁺))

# Get complex branch points
branch_pts = branch_points(decay_channel)
```
