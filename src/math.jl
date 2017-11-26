# Abbreviation to help keep code short!
const AMSA = AbstractMultiScaleArray

Base.map!{F}(f::F, m::AMSA, A0::AbstractArray, As::AbstractArray...) =
        broadcast!(f, m, A0, As...)
Base.map!{F}(f::F, m::AMSA, A0, As...) =
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

@inline function Base.Broadcast.broadcast_c{S<:AMSA}(f, ::Type{S}, A, Bs...)
    T = Base.Broadcast._broadcast_eltype(f, A, Bs...)
    shape = Base.Broadcast.broadcast_indices(A, Bs...)
    broadcast!(f, similar(A), A, Bs...)
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
