function bisect_search(a, i)
    L=1; R=length(a)
    while true
        L>R && error("Left > Right in Bisection. Failed")
        m = floor(Int, .5*(L+R))
        if a[m] == i
            return m
        elseif R==L
            return L
        elseif a[m] < i
            L = m+1
        elseif a[m] > i
            R = m
        end
    end
    #length(a) - sum(map((x)->x>=i, a))
end

Base.IndexStyle(::Type{<:AbstractMultiScaleArray}) = IndexLinear()

@inline function getindex(m::AbstractMultiScaleArray, i::Int)
    idx = bisect_search(m.end_idxs, i)
    idx > 1 && (i -= m.end_idxs[idx-1]) # also works with y
    (isempty(m.y) || idx < length(m.end_idxs)) ? m.x[idx][i] : m.y[i]
end

@inline function setindex!(m::AbstractMultiScaleArray, x, i::Int)
    idx = bisect_search(m.end_idxs, i) # +1 for 1-based indexing
    idx > 1 && (i -= m.end_idxs[idx-1])
    if isempty(m.y) || idx < length(m.end_idxs)
        m.x[idx][i] = x
    else
        m.y[i] = x
    end
end

@inline getindex(m::AbstractMultiScaleArrayLeaf, i::Int) = m.x[i]
@inline getindex(m::AbstractMultiScaleArrayLeaf, i::Int...) = m.x[i[1]]

@inline function getindex(m::AbstractMultiScaleArray, i, I...)
    if isempty(m.y) || i < length(m.end_idxs)
        if length(I)==1
            return m.x[i].x[I[1]]
        else
            return m.x[i][I...]
        end
    else
        return m.y[I...]
    end
end

@inline getindex(m::AbstractMultiScaleArrayLeaf, i...) = m.x[i[1]]
@inline getindex(m::AbstractMultiScaleArrayLeaf, i::CartesianIndex{1}) = m.x[i[1]]

@inline setindex!(m::AbstractMultiScaleArrayLeaf, x, i::Int) = (m.x[i] = x)
@inline setindex!(m::AbstractMultiScaleArray, x, i::CartesianIndex{1}) = (m[i[1]] = x)
@inline setindex!(m::AbstractMultiScaleArrayLeaf, x, i::Int...) = (m.x[i[1]] = x)

@inline function setindex!(m::AbstractMultiScaleArray, x, i, I::Int...)
    if isempty(m.y) || i < length(m.end_idxs)
        if length(I)==1
            return m.x[i].x[I[1]] = x
        else
            return m.x[i][I...] = x
        end
    else
        m.y[I...] = x
    end
end

@inline getindex(m::AbstractMultiScaleArray, ::Colon) = [m[i] for i in 1:length(m)]
@inline getindex(m::AbstractMultiScaleArrayLeaf, ::Colon) = m.x
@inline getindex(m::AbstractMultiScaleArray, i::CartesianIndex{1}) = m[i[1]] # (i, )

eachindex(m::AbstractMultiScaleArray) = 1:length(m)
endof(m::AbstractMultiScaleArray) = length(m)

Base.eltype(S::AbstractMultiScaleArray{B}) where {B} = B
#broadcast_getindex(m::AbstractMultiScaleArrayLeaf, i::Int)    =  (println("here");m[i])
#broadcast_getindex(m::AbstractMultiScaleArray, i::Int)    =  (println("here");m[i])
#broadcast_getindex(m::AbstractMultiScaleArray, i::Int...) = m[i]

getindices(m::AbstractMultiScaleArrayHead) = 1:length(m)

getindices(m::AbstractMultiScaleArrayHead, i::Int) =
    ((i > 1) ? (m.end_idxs[i-1] + 1) : 1) : m.end_idxs[i]

getindices(m::AbstractMultiScaleArrayHead, i, I::Int...) =
    getindices(m.x[i], ((i > 1) ? (m.end_idxs[i-1] + 1) : 1), m.end_idxs[i], I...)

getindices(m::AbstractMultiScaleArray, bot_idx, top_idx, i::Int) =
    (bot_idx + ((i > 1) ? m.end_idxs[i-1] : 0)) : (top_idx + m.end_idxs[i] - length(m))

getindices(m::AbstractMultiScaleArray, bot_idx, top_idx, i, I::Int...) =
    getindices(m.x[i],
               bot_idx + ((i > 1) ? m.end_idxs[i-1] : 0),
               top_idx + m.end_idxs[i] - length(m),
               I...)

function getindices(m::AbstractMultiScaleArrayLeaf, bot_idx, top_idx, i::Int)
    i > length(m) && error("Final index is larger than length of leaf")
    top_idx -= length(m) - i
    top_idx:top_idx
end

getindices(m::AbstractMultiScaleArrayLeaf, bot_idx, top_idx, I::Int...) =
    error("Indexes past the bottom.")
