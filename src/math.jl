# Abbreviation to help keep code short!
const AMSA = AbstractMultiScaleArray

Base.map!(f::F, m::AMSA, A0::AbstractArray, As::AbstractArray...) where {F} =
        broadcast!(f, m, A0, As...)
Base.map!(f::F, m::AMSA, A0, As...) where {F} =
        broadcast!(f, m, A0, As...)

Base.Broadcast.promote_containertype(::Type{T}, ::Type{T}) where {T<:AMSA} = T
Base.Broadcast.promote_containertype(::Type{T}, ::Type{S}) where {T<:AMSA, S<:AbstractArray} = T
Base.Broadcast.promote_containertype(::Type{S}, ::Type{T}) where {T<:AMSA, S<:AbstractArray} = T
Base.Broadcast.promote_containertype(::Type{T}, ::Type{<:Any}) where {T<:AMSA} = T
Base.Broadcast.promote_containertype(::Type{<:Any}, ::Type{T}) where {T<:AMSA} = T
Base.Broadcast.promote_containertype(::Type{Array}, ::Type{T}) where {T<:AMSA} = T
Base.Broadcast.promote_containertype(::Type{T}, ::Type{Array}) where {T<:AMSA} = T
Base.Broadcast._containertype(::Type{T}) where {T<:AMSA} = T
Base.Broadcast.broadcast_indices(::Type{<:AMSA}, A) = indices(A)

@inline function Base.Broadcast.broadcast_c(f, ::Type{S}, A, Bs...) where S<:AMSA
    T = Base.Broadcast._broadcast_eltype(f, A, Bs...)
    shape = Base.Broadcast.broadcast_indices(A, Bs...)
    broadcast!(f, similar(A), A, Bs...)
end

@inline function Base.Broadcast.broadcast_c(f, ::Type{S}, A::AMSA, Bs::AMSA...) where S<:AMSA
    new_A = similar(A)
    broadcast!(f,new_A,A,Bs...)
    new_A
end

@inline function Base.Broadcast.broadcast_c!(f, ::Type{S}, ::Type, A::AbstractMultiScaleArrayLeaf, Bs::AbstractMultiScaleArrayLeaf...) where S<:AbstractMultiScaleArrayLeaf
    broadcast!(f, A.values, (B.values for B in Bs)...)
    A
end

@inline function Base.Broadcast.broadcast_c!(f, ::Type{S}, ::Type, A::AMSA, Bs::AMSA...) where S<:AMSA
    for i in eachindex(A.nodes)
            broadcast!(f, A.nodes[i], (B.nodes[i] for B in Bs)...)
    end
    broadcast!(f, A.values, (B.values for B in Bs)...)
    A
end

+(m::AbstractMultiScaleArray, y::AbstractMultiScaleArray) = m .+ y
+(m::AbstractMultiScaleArray, y::Number) = m .+ y
+(y::Number, m::AbstractMultiScaleArray) = m .+ y

-(m::AbstractMultiScaleArray, y::AbstractMultiScaleArray) = m .- y
-(m::AbstractMultiScaleArray, y::Number) = m .- y
-(y::Number, m::AbstractMultiScaleArray) = y .- m

*(m::AbstractMultiScaleArray, y::Number) = m .* y
*(y::Number, m::AbstractMultiScaleArray) = m .* y

/(m::AbstractMultiScaleArray, y::Number) = m ./ y
/(y::Number, m::AbstractMultiScaleArray) = y ./ m
