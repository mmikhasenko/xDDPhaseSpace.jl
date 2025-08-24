# γDD
struct γDD{T1, T2, T3} <: AbstractxDD
    ms::NamedTuple{(:m1, :m2, :m3), T1}
    R12::T2
    R13::T3
    μ12::Float64 # Transition matrix element < D | mu | D* >
    μ13::Float64 # Transition matrix element < D | mu | D* >
end

const h² = 20.13e-3

function γDD_𝔐²_nonana3(v, (J₁₂ᴵ, J₁₂ᴵᴵ), (J₁₃ᴵ, J₁₃ᴵᴵ), μ12, μ13)
    _p1_p2 = p1_p2(v)
    _p1_p3 = p1_p3(v)
    _G = G(v)
    #
    𝔐² = μ12^2 * (_p1_p2^2 + _G) * J₁₂ᴵ * J₁₂ᴵᴵ -
         μ12 * μ13 * (_p1_p2 * _p1_p3 - _G) * J₁₃ᴵ * J₁₂ᴵᴵ
    # 	
    return h² * 𝔐² / 3
end

function γDD_𝔐²_nonana2(v, (J₁₂ᴵ, J₁₂ᴵᴵ), (J₁₃ᴵ, J₁₃ᴵᴵ), μ12, μ13)
    _p1_p2 = p1_p2(v)
    _p1_p3 = p1_p3(v)
    _G = G(v)
    # 
    𝔐² = μ13^2 * (_p1_p3^2 + _G) * J₁₃ᴵ * J₁₃ᴵᴵ -
         μ12 * μ13 * (_p1_p2 * _p1_p3 - _G) * J₁₂ᴵ * J₁₃ᴵᴵ
    # 
    return h² * 𝔐² / 3
end

decay_matrix_element_squared(d::γDD, s, σ3, σ2) = covertapply(
    (v, J₁₂, J₁₃) -> γDD_𝔐²_nonana3(v, J₁₂, J₁₃, d.μ12, d.μ13) +
                     γDD_𝔐²_nonana2(v, J₁₂, J₁₃, d.μ12, d.μ13),
    d, s, σ3, σ2)

branch_points(d::γDD) = (
    d.ms[3] + sqrt(pole_position(d.R12)),
    d.ms[2] + sqrt(pole_position(d.R13)))

function ρ_tb(d::γDD, e::Real)
    M, m = d.R13.m, d.ms.m2
    sqrts = e2m(e)
    sqrts < M + m ? 0.0 :
    sqrt(λ(e2m(e)^2, M^2, m^2)) / e2m(e)^2
end

#                                _|  _|  _|                        _|      _|                      
#    _|_|_|    _|_|    _|  _|_|      _|      _|_|_|_|    _|_|_|  _|_|_|_|        _|_|    _|_|_|    
#  _|_|      _|_|_|_|  _|_|      _|  _|  _|      _|    _|    _|    _|      _|  _|    _|  _|    _|  
#      _|_|  _|        _|        _|  _|  _|    _|      _|    _|    _|      _|  _|    _|  _|    _|  
#  _|_|_|      _|_|_|  _|        _|  _|  _|  _|_|_|_|    _|_|_|      _|_|  _|    _|_|    _|    _|  

#
obj2nt(ch::γDD) =
    (type = string(typeof(ch)),
        ms = ch.ms,
        R12 = obj2nt(ch.R12), R13 = obj2nt(ch.R13),
        μ12 = ch.μ12, μ13 = ch.μ13)
#
