module xDDPhaseSpace

using Cuba
using Parameters

import Base.^
^(ms::NamedTuple{(:m1, :m2, :m3),T} where {T}, n::Int) = Tuple(ms) .^ n  # type piracy -- not good

export WithThrE
export e2m
include("withthr.jl")

export AbstractxDD
export mapdalitzmethod
export integrand_mapped_thr
export ρ_thr, ρ_tb
export masses
include("abstractxdd.jl")

export BW, BW_norm, BW_Swave, ZeroBW
export Jᴵ, Jᴵᴵ
include("isobar.jl")

export AbstractDalitzMapping
export HookSqrtDalitzMapping
export LinearDalitzMapping
export mapdalitz
include("integrationpaths.jl")

export branch_points
export πDD
include("covariant_exptessions.jl")
include("mainmodel-pi.jl")
export γDD
include("mainmodel-gamma.jl")

export obj2nt
export constructchannel
include("deserialization.jl")

end # xDDPhaseSpace
