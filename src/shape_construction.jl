length(m::AbstractMultiScaleArrayLeaf) = length(m.x)
length(m::AbstractMultiScaleArray) = m.end_idxs[end]
num_daughters(m::AbstractMultiScaleArrayLeaf) = 0
num_daughters(m::AbstractMultiScaleArray) = size(m.x)[1]
ndims(m::AbstractMultiScaleArray) = 1
size(m::AbstractMultiScaleArray, i::Int) = i == 1 ? length(m) : 0
size(m::AbstractMultiScaleArray) = (length(m),)

parameterless_type(T::Type) = Base.typename(T).wrapper
parameterless_type(x) = parameterless_type(typeof(x))

similar(m::AbstractMultiScaleArray) =
    construct(typeof(m), deepcopy(m.x), deepcopy(m.y))
similar(m::AbstractMultiScaleArrayLeaf) =
    construct(typeof(m),similar(m.x))
similar(m::AbstractMultiScaleArrayLeaf,T::Type) =
    construct(parameterless_type(m),similar(m.x,T))
similar(m::AbstractMultiScaleArray,T::Type) =
    construct(parameterless_type(m), [similar(v,T) for v in m.x], similar(m.y,T))

construct(::Type{T},x) where {T<:AbstractMultiScaleArrayLeaf} = T(x)

function __construct(x::Vector{<:AbstractMultiScaleArray})
    end_idxs = Vector{Int}(length(x))
    off = 0
    @inbounds for i in 1:length(x)
        end_idxs[i] = (off += length(x[i]))
    end
    end_idxs
end

construct(::Type{T}, x::Vector{<:AbstractMultiScaleArray}) where {T<:AbstractMultiScaleArray} =
    T(x, Float64[], __construct(x))

function construct(::Type{T}, x::Vector{<:AbstractMultiScaleArray},
                   y::Vector{Float64}) where {T<:AbstractMultiScaleArray}
    ylen = length(y)
    end_idxs = Vector{Int}(length(x) + ifelse(ylen == 0, 0, 1))
    off = 0
    @inbounds for i in 1:length(x)
        end_idxs[i] = (off += length(x[i]))
    end
    ylen == 0 || (end_idxs[end] = off + ylen)
    T(x, y, end_idxs)
end

vcat(m1::AbstractMultiScaleArray,m2::AbstractMultiScaleArray) =
    error("AbstractMultiScaleArrays cannot be concatenated")

hcat(m1::AbstractMultiScaleArray,m2::AbstractMultiScaleArray) =
    error("AbstractMultiScaleArrays cannot be concatenated")

==(m1::AbstractMultiScaleArray,m2::AbstractMultiScaleArray) = (m1 === m2)

function recursivecopy!(b::AbstractMultiScaleArrayLeaf, a::AbstractMultiScaleArrayLeaf)
    @inbounds copy!(b,a)
end

function recursivecopy!(b::AbstractMultiScaleArray, a::AbstractMultiScaleArray)
    @inbounds for i in eachindex(a.x)
        recursivecopy!(b.x[i], a.x[i])
    end
    @inbounds for i in eachindex(a.y)
        recursivecopy!(b.y[i], a.y[i])
    end
end
