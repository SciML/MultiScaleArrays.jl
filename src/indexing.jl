nodeselect(ns, i, I...) = ns[i][I...]
nodechild(ns, i, j) = ns[i].nodes[j]

function bisect_search(a, i)
    first(searchsorted(a, i))
end

Base.IndexStyle(::Type{<:AbstractMultiScaleArray}) = IndexLinear()

@inline function getindex(m::AbstractMultiScaleArray, i::Int)
    idx = bisect_search(m.end_idxs, i)
    idx > 1 && (i -= m.end_idxs[idx - 1]) # also works with values
    (isempty(m.values) || idx < length(m.end_idxs)) ? nodeselect(m.nodes, idx, i) :
    m.values[i]
end

@inline function setindex!(m::AbstractMultiScaleArray, nodes, i::Int)
    idx = bisect_search(m.end_idxs, i) # +1 for 1-based indexing
    idx > 1 && (i -= m.end_idxs[idx - 1])
    if isempty(m.values) || idx < length(m.end_idxs)
        m.nodes[idx][i] = nodes
    else
        m.values[i] = nodes
    end
end

@inline getindex(m::AbstractMultiScaleArrayLeaf, i::Int) = m.values[i]
@inline getindex(m::AbstractMultiScaleArrayLeaf, i::Int...) = m.values[i[1]]

@inline function getindex(m::AbstractMultiScaleArray, i, I...)
    if isempty(m.values) || i < length(m.end_idxs)
        length(I) == 1 ? nodechild(m.nodes, i, I[1]) : nodeselect(m.nodes, i, I...)
    else
        m.values[I...]
    end
end

@inline getindex(m::AbstractMultiScaleArrayLeaf, i...) = m.values[i[1]]
@inline getindex(m::AbstractMultiScaleArrayLeaf, i::CartesianIndex{1}) = m.values[i[1]]

@inline setindex!(m::AbstractMultiScaleArrayLeaf, values, i::Int) = (m.values[i] = values)
@inline function setindex!(m::AbstractMultiScaleArrayLeaf, values, i::Int...)
    (m.values[i[1]] = values)
end

@inline function setindex!(m::AbstractMultiScaleArray, values, i::CartesianIndex{1})
    (m[i[1]] = values)
end

@inline function setindex!(m::AbstractMultiScaleArray, values, i, I::Int...)
    if isempty(m.values) || i < length(m.end_idxs)
        mx = m.nodes[i]
        length(I) == 1 ? (mx[I[1]] = values) : (mx[I...] = values)
    else
        m.values[I...] = values
    end
end

@inline getindex(m::AbstractMultiScaleArray, ::Colon) = [m[i] for i in 1:length(m)]
@inline getindex(m::AbstractMultiScaleArrayLeaf, ::Colon) = m.values
@inline getindex(m::AbstractMultiScaleArray, i::CartesianIndex{1}) = m[i[1]] # (i, )

eachindex(m::AbstractMultiScaleArray) = 1:length(m)
endof(m::AbstractMultiScaleArray) = length(m)

Base.eltype(::Type{<:AbstractMultiScaleArray{B}}) where {B <: Number} = B
# For now, default to Float64 (like v0.5 version did) if not specified
function Base.eltype(::Type{T}) where {T <: AbstractMultiScaleArray}
    eltype(fieldtype(T, :values)) === Any && return Float64
    error("Invalid AbstractMultiScaleArray, type of values not specified")
end
Base.eltype(::T) where {T <: AbstractMultiScaleArray} = eltype(T)

#broadcast_getindex(m::AbstractMultiScaleArrayLeaf, i::Int)    =  (println("here");m[i])
#broadcast_getindex(m::AbstractMultiScaleArray, i::Int)    =  (println("here");m[i])
#broadcast_getindex(m::AbstractMultiScaleArray, i::Int...) = m[i]

getindices(m::AbstractMultiScaleArrayHead) = 1:length(m)

function getindices(m::AbstractMultiScaleArrayHead, i::Int)
    ((i > 1) ? (m.end_idxs[i - 1] + 1) : 1):m.end_idxs[i]
end

function getindices(m::AbstractMultiScaleArrayHead, i, I::Int...)
    getindices(m.nodes[i], ((i > 1) ? (m.end_idxs[i - 1] + 1) : 1), m.end_idxs[i], I...)
end

function getindices(m::AbstractMultiScaleArray, bot_idx, top_idx, i::Int)
    (bot_idx + ((i > 1) ? m.end_idxs[i - 1] : 0)):(top_idx + m.end_idxs[i] - length(m))
end

function getindices(m::AbstractMultiScaleArray, bot_idx, top_idx, i, I::Int...)
    getindices(m.nodes[i],
        bot_idx + ((i > 1) ? m.end_idxs[i - 1] : 0),
        top_idx + m.end_idxs[i] - length(m),
        I...)
end

function getindices(m::AbstractMultiScaleArrayLeaf, bot_idx, top_idx, i::Int)
    i > length(m) && error("Final index is larger than length of leaf")
    top_idx -= length(m) - i
    top_idx:top_idx
end

function getindices(m::AbstractMultiScaleArrayLeaf, bot_idx, top_idx, I::Int...)
    error("Indexes past the bottom.")
end
