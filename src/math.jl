# Abbreviation to help keep code short!
const AMSA = AbstractMultiScaleArray

Base.map!(f::F, m::AMSA, A0::AbstractArray, As::AbstractArray...) where {F} =
        broadcast!(f, m, A0, As...)
Base.map!(f::F, m::AMSA, A0, As...) where {F} =
        broadcast!(f, m, A0, As...)

struct AMSAStyle{Style <: Broadcast.BroadcastStyle} <: Broadcast.AbstractArrayStyle{Any} end
AMSAStyle(::S) where {S} = AMSAStyle{S}()
AMSAStyle(::S, ::Val{N}) where {S,N} = AMSAStyle(S(Val(N)))
AMSAStyle(::Val{N}) where N = AMSAStyle{Broadcast.DefaultArrayStyle{N}}()

# promotion rules
function Broadcast.BroadcastStyle(::AMSAStyle{AStyle}, ::AMSAStyle{BStyle}) where {AStyle, BStyle}
    AMSAStyle(Broadcast.BroadcastStyle(AStyle(), BStyle()))
end

#=
combine_styles(args::Tuple{})         = Broadcast.DefaultArrayStyle{0}()
combine_styles(args::Tuple{Any})      = Broadcast.result_style(Broadcast.BroadcastStyle(args[1]))
combine_styles(args::Tuple{Any, Any}) = Broadcast.result_style(Broadcast.BroadcastStyle(args[1]), Broadcast.BroadcastStyle(args[2]))
@inline combine_styles(args::Tuple)   = Broadcast.result_style(Broadcast.BroadcastStyle(args[1]), combine_styles(Base.tail(args)))

function Broadcast.BroadcastStyle(::Type{AMSA{T}}) where {T}
    Style = combine_styles((T.parameters...,))
    AMSAStyle(Style)
end
=#

@inline function Base.copy(bc::Broadcast.Broadcasted{AMSAStyle{Style}}) where Style
    nnodes(bc)
    @inline function f(i)
        copy(unpack(bc, i))
    end
    CommonType = get_common_type(bc)
    construct(CommonType, map(f,N), f(nothing))
end

@inline function Base.copyto!(dest::AMSA, bc::Broadcast.Broadcasted)
    N = length(dest.nodes)
    for i in 1:N
        copyto!(dest.nodes[i], unpack(bc, i))
    end
    copyto!(dest.values,unpack(bc, nothing))
end

@inline function Base.copyto!(dest::AbstractMultiScaleArrayLeaf, bc::Broadcast.Broadcasted)
    copyto!(dest.values,unpack(bc,nothing))
end

# drop axes because it is easier to recompute
@inline unpack(bc::Broadcast.Broadcasted{Style}, i) where Style = Broadcast.Broadcasted{Style}(bc.f, unpack_args(i, bc.args))
@inline unpack(bc::Broadcast.Broadcasted{AMSAStyle{Style}}, i) where Style = Broadcast.Broadcasted{Style}(bc.f, unpack_args(i, bc.args))
unpack(x,::Any) = x
unpack(x::AMSA, i) = x.nodes[i]
unpack(x::AMSA, ::Nothing) = x.values

@inline unpack_args(i, args::Tuple) = (unpack(args[1], i), unpack_args(i, Base.tail(args))...)
unpack_args(i, args::Tuple{Any}) = (unpack(args[1], i),)
unpack_args(::Any, args::Tuple{}) = ()

nnodes(A) = 0
nnodes(A::AMSA) = length(A.nodes)
nnodes(bc::Broadcast.Broadcasted) = _nnodes(bc.args)
nnodes(A, Bs...) = common_number(nnodes(A), _nnodes(Bs))

@inline _nnodes(args::Tuple) = common_number(nnodes(args[1]), _nnodes(Base.tail(args)))
_nnodes(args::Tuple{Any}) = nnodes(args[1])
_nnodes(args::Tuple{}) = 0

get_common_type(A) = Nothing
get_common_type(A::AMSA) = typeof(A)
get_common_type(bc::Broadcast.Broadcasted) = _nnodes(bc.args)
get_common_type(A, Bs...) = common_type(get_common_type(A), _get_common_type(Bs))

@inline _get_common_type(args::Tuple) = get_common_type(get_common_type(args[1]), _get_common_type(Base.tail(args)))
_get_common_type(args::Tuple{Any}) = get_common_type(args[1])
_get_common_type(args::Tuple{}) = Nothing

## utils
common_number(a, b) =
    a == 0 ? b :
    (b == 0 ? a :
     (a == b ? a :
      throw(DimensionMismatch("number of nodes must be equal"))))

common_type(a::T, b::T) where T = T
common_type(a::T, b::Nothing) where T = T
common_type(a::Nothing, b::T) where T = T
common_type(a::Nothing, b::Nothing) where T = Nothing
