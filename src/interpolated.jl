"""
    interpolated{T<:AbstractxDD,F1<:Real,F2<:Real,V}

Wrapper around an `AbstractxDD` channel that tabulates `ρ_thr` from the
three-body threshold up to a chosen `cutoff` and uses linear interpolation in
that region.

Above the cutoff, the threshold density is continued with the two-body proxy
`ρ_tb`, rescaled so that the interpolated and asymptotic descriptions match at
`cutoff`.
"""
struct interpolated{T<:AbstractxDD,F1<:Real,F2<:Real,V}
    channel::T
    cutoff::F1
    cutoffratio::F2
    itr::V
end

"""
    interpolated(d::AbstractxDD, cutoff::Real; estep_GeV=0.01e-3)

Construct an interpolated representation of the threshold phase-space density
for channel `d`.

The function samples `ρ_thr(d, m)` on a uniform mass grid from the threshold
`sum(masses(d))` to slightly above `cutoff`, builds a linear interpolant, and
stores a matching factor

```julia
real(ρ_thr(d, cutoff)) / ρ_tb(d, cutoff)
```

so that for `m >= cutoff` the interpolated channel follows `ρ_tb(d, m)` with a
continuous normalization.

# Arguments
- `d`: underlying decay channel.
- `cutoff`: upper end of the interpolation window in GeV.

# Keywords
- `estep_GeV`: spacing of the interpolation grid in GeV.
"""
function interpolated(d::AbstractxDD, cutoff::Real; estep_GeV=0.01e-3)
    mth = masses(d) |> sum
    mv = mth:estep_GeV:(cutoff+2*estep_GeV)
    calv = ρ_thr.(Ref(d), mv)
    itr = interpolate((mv,), calv, Gridded(Linear()))
    cutoffratio = real(ρ_thr(d, cutoff)) / ρ_tb(d, cutoff)
    interpolated(d, cutoff, cutoffratio, itr)
end

"""
    ρ_thr(d::interpolated, m::Real)

Evaluate the threshold phase-space density for an interpolated channel.

- For `m < sum(masses(d.channel))`, the result is `0.0`.
- For threshold `<= m < d.cutoff`, the value is taken from the stored linear
  interpolant.
- For `m >= d.cutoff`, the function returns `ρ_tb(d.channel, m) *
  d.cutoffratio`, which preserves the normalization at the cutoff.
"""
function ρ_thr(d::interpolated, m::Real)
    mth = masses(d.channel) |> sum
    m < mth && return 0.0
    m < d.cutoff ?
    d.itr(m) :
    ρ_tb(d.channel, m) * d.cutoffratio
end
ρ_thr(d::interpolated, m::Complex) = ρ_thr(d.channel, m)



"""
    dispersive(d::interpolated, m)

Calculate the dispersive part associated with an interpolated channel.

This evaluates the dispersion integral built from `ρ_thr(d, m')` and returns
the first-sheet analytic continuation as a function of the external mass `m`.
For real `m`, the method below evaluates the integral slightly above the real
axis, `m + iϵ`, so that the expected imaginary part is recovered numerically.

# Arguments
- `d::interpolated`: The interpolated phase space density.
- `m`: The mass at which to calculate the dispersive part.
"""
function dispersive(d::interpolated, m)
    s = m^2
    mth = masses(d.channel) |> sum
    function integrand(s′)
        m′ = sqrt(s′)
        ρ_thr(d, m′) / s′ / (s′ - s)
    end
    s / π * quadgk(integrand, mth^2, Inf)[1]
end
# 
dispersive(d::interpolated, m::Real) = dispersive(d, m + iϵ)
