length(m::AbstractMultiScaleArrayLeaf) = length(m.nodes)
length(m::AbstractMultiScaleArray) = m.end_idxs[end]
num_daughters(m::AbstractMultiScaleArrayLeaf) = 0
num_daughters(m::AbstractMultiScaleArray) = size(m.nodes)[1]
ndims(m::AbstractMultiScaleArray) = 1
size(m::AbstractMultiScaleArray, i::Int) = i == 1 ? length(m) : 0
size(m::AbstractMultiScaleArray) = (length(m),)

parameterless_type(T::Type) = Base.typename(T).wrapper
parameterless_type(x) = parameterless_type(typeof(x))

similar(m::AbstractMultiScaleArray) =
    construct(typeof(m), deepcopy(m.nodes), deepcopy(m.values))
similar(m::AbstractMultiScaleArrayLeaf) =
    construct(typeof(m), similar(m.nodes))
similar(m::AbstractMultiScaleArrayLeaf, T::Type) =
    construct(parameterless_type(m), similar(m.nodes,T))
similar(m::AbstractMultiScaleArray,T::Type) =
    construct(parameterless_type(m), [similar(v, T) for v in m.nodes], similar(m.values, T))

construct(::Type{T}, x) where {T<:AbstractMultiScaleArrayLeaf} = T(x)

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

function construct(::Type{T}, nodes::Vector{<:AbstractMultiScaleArray},
                   values::Vector{Float64}) where {T<:AbstractMultiScaleArray}
    vallen = length(values)
    end_idxs = Vector{Int}(length(nodes) + ifelse(vallen == 0, 0, 1))
    off = 0
    @inbounds for i in 1:length(nodes)
        end_idxs[i] = (off += length(nodes[i]))
    end
    vallen == 0 || (end_idxs[end] = off + vallen)
    T(nodes, values, end_idxs)
end

vcat(m1::AbstractMultiScaleArray, m2::AbstractMultiScaleArray) =
    error("AbstractMultiScaleArrays cannot be concatenated")

hcat(m1::AbstractMultiScaleArray, m2::AbstractMultiScaleArray) =
    error("AbstractMultiScaleArrays cannot be concatenated")

==(m1::AbstractMultiScaleArray, m2::AbstractMultiScaleArray) = (m1 === m2)

function recursivecopy!(b::AbstractMultiScaleArrayLeaf, a::AbstractMultiScaleArrayLeaf)
    @inbounds copy!(b,a)
end

function recursivecopy!(b::AbstractMultiScaleArray, a::AbstractMultiScaleArray)
    @inbounds for i in eachindex(a.nodes)
        recursivecopy!(b.nodes[i], a.nodes[i])
    end
    @inbounds for i in eachindex(a.values)
        recursivecopy!(b.values[i], a.values[i])
    end
end
