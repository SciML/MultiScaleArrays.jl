function __add_daughter!(m::AbstractMultiScaleArray, x::AbstractMultiScaleArray)
    push!(m.end_idxs, m.end_idxs[end] + length(x) + length(m.y))
    push!(m.x, x)
    nothing
end

function __update_lengths(m::AbstractMultiScaleArray, beg, len)
    for j = beg:num_daughters(m)
        m.end_idxs[j] += len
    end
    isempty(m.y) || (m.end_idxs[end] += len)
    nothing
end

function __add_daughter!(m::AbstractMultiScaleArray, x::AbstractMultiScaleArray, i::Int)
    __add_daughter!(m.x[i], x)
    __update_lengths(m, i, length(x))
end

function __add_daughter!(m::AbstractMultiScaleArray, x::AbstractMultiScaleArray, i, I::Int...)
    __add_daughter!(m.x[i], x, I...)
    __update_lengths(m, i, length(x))
end

add_daughter!(m::AbstractMultiScaleArrayHead, x::AbstractMultiScaleArray) =
    __add_daughter!(m, x)

add_daughter!(m::AbstractMultiScaleArrayHead, x::AbstractMultiScaleArray, i::Int) =
    __add_daughter!(m, x, i)

add_daughter!(m::AbstractMultiScaleArrayHead, x::AbstractMultiScaleArray, i, I::Int...) =
    __add_daughter!(m, x, i, I...)

function __remove_daughter!(m::AbstractMultiScaleArray, i::Int)
    del_length = length(m.x[i])
    deleteat!(m.x, i)
    deleteat!(m.end_idxs, i)
    for j = i:num_daughters(m)
        m.end_idxs[j] -= del_length
    end
    isempty(m.y) || (m.end_idxs[end] -= del_length)
    del_length
end

function __remove_daughter!(m::AbstractMultiScaleArrayLeaf, i::Int)
    deleteat!(m.x, i)
    1
end

function remove_daughter!(m::AbstractMultiScaleArrayHead, i, I::Int...)
    del_length = __remove_daughter!(m.x[i], I...)
    for j = i:num_daughters(m)
        m.end_idxs[j] -= del_length
    end
    isempty(m.y) || (m.end_idxs[end] -= del_length)
    if size(m.x[i].x) == (0, )
        deleteat!(m.x, i)
        deleteat!(m.end_idxs, i)
    end
    nothing
end

function __remove_daughter!(m::AbstractMultiScaleArray, i, I::Int...)
    del_length = __remove_daughter!(m.x[i], I...)
    for j = i:num_daughters(m)
        m.end_idxs[j] -= del_length
    end
    isempty(m.y) || (m.end_idxs[end] -= del_length)
    isempty(m.x[i]) && deleteat!(m.x, i)
    del_length
end

remove_daughter!(m::AbstractMultiScaleArrayHead, i::Int) =
    (__remove_daughter!(m, i); nothing)
