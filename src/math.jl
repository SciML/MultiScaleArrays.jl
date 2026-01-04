# Abbreviation to help keep code short!
const AMSA = AbstractMultiScaleArray

function Base.map!(f::F, m::AMSA, A0::AbstractArray, As::AbstractArray...) where {F}
    return broadcast!(f, m, A0, As...)
end
Base.map!(f::F, m::AMSA, A0, As...) where {F} = broadcast!(f, m, A0, As...)

struct AMSAStyle <: Broadcast.AbstractArrayStyle{Any} end
Broadcast.BroadcastStyle(::AMSAStyle, ::Broadcast.DefaultArrayStyle{0}) = AMSAStyle()
Broadcast.BroadcastStyle(::Broadcast.DefaultArrayStyle{0}, ::AMSAStyle) = AMSAStyle()
function Broadcast.BroadcastStyle(::AMSAStyle, ::Broadcast.DefaultArrayStyle{N}) where {N}
    return Broadcast.DefaultArrayStyle{N}()
end
Broadcast.BroadcastStyle(::Type{<:AMSA}) = AMSAStyle()

@inline function Base.copy(bc::Broadcast.Broadcasted{<:AMSAStyle})
    first_amsa = find_amsa(bc)::AMSA  # Type assertion: AMSAStyle guarantees an AMSA exists

    out = similar(first_amsa, Base.Broadcast._broadcast_getindex_eltype(bc))

    #=
    ElType = Base.Broadcast.combine_eltypes(bc.f, bc.args)
    if Base.isconcretetype(ElType)
        # We can trust it and defer to the simpler `copyto!`
        return copyto!(similar(bc, ElType), bc)
    end
    =#

    copyto!(out, bc)
    return out
end

@inline function Base.copyto!(dest::AMSA, bc::Broadcast.Broadcasted{<:AMSAStyle})
    N = length(dest.nodes)
    for i in 1:N
        copyto!(dest.nodes[i], unpack(bc, i))
    end
    copyto!(dest.values, unpack(bc, nothing))
    return dest
end

@inline function Base.copyto!(
        dest::AbstractMultiScaleArrayLeaf,
        bc::Broadcast.Broadcasted{<:AMSAStyle}
    )
    copyto!(dest.values, unpack(bc, nothing))
    return dest
end

# drop axes because it is easier to recompute
@inline function unpack(bc::Broadcast.Broadcasted, i)
    return Broadcast.Broadcasted(bc.f, unpack_args(i, bc.args))
end
unpack(x, ::Any) = x
unpack(x::AMSA, i) = x.nodes[i]
unpack(x::AMSA, ::Nothing) = x.values

@inline function unpack_args(i, args::Tuple)
    return (unpack(args[1], i), unpack_args(i, Base.tail(args))...)
end
unpack_args(i, args::Tuple{Any}) = (unpack(args[1], i),)
unpack_args(::Any, args::Tuple{}) = ()

nnodes(A) = 0
nnodes(A::AMSA) = length(A.nodes)
nnodes(bc::Broadcast.Broadcasted) = _nnodes(bc.args)
nnodes(A, Bs...) = common_number(nnodes(A), _nnodes(Bs))

@inline _nnodes(args::Tuple) = common_number(nnodes(args[1]), _nnodes(Base.tail(args)))
_nnodes(args::Tuple{Any}) = nnodes(args[1])
_nnodes(args::Tuple{}) = 0

"""
`A = find_amsa(As)` returns the first AMSA among the arguments.

Returns the first `AbstractMultiScaleArray` found in the arguments, or `nothing` if none found.
"""
find_amsa(bc::Base.Broadcast.Broadcasted) = find_amsa(bc.args)
find_amsa(args::Tuple) = isempty(args) ? nothing : find_amsa(find_amsa(args[1]), Base.tail(args))
find_amsa(x) = nothing
find_amsa(x::AMSA) = x
find_amsa(a::AMSA, rest) = a
find_amsa(::Any, rest) = find_amsa(rest)
find_amsa(::Nothing, rest) = find_amsa(rest)

any_non_amsa(bc::Base.Broadcast.Broadcasted) = any_non_amsa(bc.args)
any_non_amsa(args::Tuple) = any_non_amsa(any_non_amsa(args[1]), Base.tail(args))
any_non_amsa(x::AMSA) = false
any_non_amsa(x::Number) = false
any_non_amsa(x::Any) = true
any_non_amsa(x::AbstractArray) = true
any_non_amsa(x::Bool, rest) = isempty(rest) ? x : x || any_non_amsa(rest)

## utils
function common_number(a, b)
    return a == 0 ? b :
        (
            b == 0 ? a :
            (
                a == b ? a :
                throw(DimensionMismatch("number of nodes must be equal"))
            )
        )
end

## Linear Algebra

function LinearAlgebra.ldiv!(A::LinearAlgebra.LU, b::AMSA)
    x = Array(b)
    ldiv!(A, x)
    return b .= x
end
function LinearAlgebra.ldiv!(A::LinearAlgebra.QR, b::AMSA)
    x = Array(b)
    ldiv!(A, x)
    return b .= x
end
function LinearAlgebra.ldiv!(A::LinearAlgebra.SVD, b::AMSA)
    x = Array(b)
    ldiv!(A, x)
    return b .= x
end
