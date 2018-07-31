# Abbreviation to help keep code short!
const AMSA = AbstractMultiScaleArray

Base.map!(f::F, m::AMSA, A0::AbstractArray, As::AbstractArray...) where {F} =
        broadcast!(f, m, A0, As...)
Base.map!(f::F, m::AMSA, A0, As...) where {F} =
        broadcast!(f, m, A0, As...)

const AMSAStyle = Broadcast.ArrayStyle{AMSA}
Base.BroadcastStyle(::Type{<:AMSA}) = Broadcast.ArrayStyle{AMSA}()
Base.BroadcastStyle(::Broadcast.ArrayStyle{AMSA},::Broadcast.DefaultArrayStyle{1}) = Broadcast.DefaultArrayStyle{1}()
Base.BroadcastStyle(::Broadcast.DefaultArrayStyle{1},::Broadcast.ArrayStyle{AMSA}) = Broadcast.DefaultArrayStyle{1}()
Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{AMSA}},::Type{ElType}) where ElType = similar(bc)

function Base.copy(bc::Broadcast.Broadcasted{AMSAStyle})
    ret = Broadcast.flatten(bc)
    __broadcast(ret.f,ret.args...)
end

function Base.copyto!(dest::AMSA, bc::Broadcast.Broadcasted{AMSAStyle})
    ret = Broadcast.flatten(bc)
    __broadcast!(ret.f,dest,ret.args...)
end

function __broadcast(f, A::AMSA, Bs...)
    broadcast!(f, similar(A), A, Bs...)
end

function __broadcast!(f, A::AbstractMultiScaleArrayLeaf, Bs::Union{Number,AbstractMultiScaleArrayLeaf}...)
    broadcast!(f, A.values, (typeof(B)<:AbstractMultiScaleArrayLeaf ? B.values : B for B in Bs)...)
    A
end

function __broadcast!(f, A::AMSA, Bs::Union{Number,AMSA}...)
    for i in eachindex(A.nodes)
            broadcast!(f, A.nodes[i], (typeof(B)<:AMSA ? B.nodes[i] : B for B in Bs)...)
    end
    broadcast!(f, A.values, (typeof(B)<:AMSA ? B.values : B for B in Bs)...)
    A
end

+(m::AbstractMultiScaleArray, y::Number) = m .+ y
+(y::Number, m::AbstractMultiScaleArray) = m .+ y

-(m::AbstractMultiScaleArray, y::Number) = m .- y
-(y::Number, m::AbstractMultiScaleArray) = y .- m

*(m::AbstractMultiScaleArray, y::Number) = m .* y
*(y::Number, m::AbstractMultiScaleArray) = m .* y

/(m::AbstractMultiScaleArray, y::Number) = m ./ y
/(y::Number, m::AbstractMultiScaleArray) = y ./ m
