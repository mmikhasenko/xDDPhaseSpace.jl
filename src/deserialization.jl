obj2nt(nt::NamedTuple) = nt

# deserilization
function constructchannel(ser::NamedTuple{(:R12, :R13, :ms, :type), T} where {T})
    @assert startswith(ser.type, "πDD")
    type = eval(Meta.parse(ser.type))
    ms = NamedTuple{(:m1, :m2, :m3)}(ser.ms)
    R12 = constructlineshape(ser.R12)
    R13 = constructlineshape(ser.R13)
    type(ms, R12, R13)
end

function constructchannel(ser::NamedTuple{(:R12, :R13, :μ12, :μ13, :ms, :type), T} where {T})
    @assert startswith(ser.type, "γDD")
    type = eval(Meta.parse(ser.type))
    ms = NamedTuple{(:m1, :m2, :m3)}(ser.ms)
    R12 = constructlineshape(ser.R12)
    R13 = constructlineshape(ser.R13)
    type(ms, R12, R13, ser.μ12, ser.μ13)
end
constructlineshape(R::NamedTuple{(:m, :Γ), T} where {T}) = R
# 
function constructlineshape(R::NamedTuple)
    type = eval(Meta.parse(R.type))
    keysnotype = filter(x -> x != :type, keys(R))
    args = NamedTuple{keysnotype}(R)
    type(args...)
end
