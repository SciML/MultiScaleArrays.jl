function __add_daughter!(m::AbstractMultiScaleArray, nodes::AbstractMultiScaleArray)
    push!(m.end_idxs, m.end_idxs[end] + length(nodes) + length(m.values))
    push!(m.nodes, nodes)
    nothing
end

function __update_lengths(m::AbstractMultiScaleArray, beg, len)
    for j = beg:num_daughters(m)
        m.end_idxs[j] += len
    end
    isempty(m.values) || (m.end_idxs[end] += len)
    nothing
end

function __add_daughter!(m::AbstractMultiScaleArray, nodes::AbstractMultiScaleArray, i::Int)
    __add_daughter!(m.nodes[i], nodes)
    __update_lengths(m, i, length(nodes))
end

function __add_daughter!(m::AbstractMultiScaleArray, nodes::AbstractMultiScaleArray, i, I::Int...)
    __add_daughter!(m.nodes[i], nodes, I...)
    __update_lengths(m, i, length(nodes))
end

add_daughter!(m::AbstractMultiScaleArrayHead, nodes::AbstractMultiScaleArray) =
    __add_daughter!(m, nodes)

add_daughter!(m::AbstractMultiScaleArrayHead, nodes::AbstractMultiScaleArray, i::Int) =
    __add_daughter!(m, nodes, i)

add_daughter!(m::AbstractMultiScaleArrayHead, nodes::AbstractMultiScaleArray, i, I::Int...) =
    __add_daughter!(m, nodes, i, I...)

function __remove_daughter!(m::AbstractMultiScaleArray, i::Int)
    del_length = length(m.nodes[i])
    deleteat!(m.nodes, i)
    deleteat!(m.end_idxs, i)
    for j = i:num_daughters(m)
        m.end_idxs[j] -= del_length
    end
    isempty(m.values) || (m.end_idxs[end] -= del_length)
    del_length
end

function __remove_daughter!(m::AbstractMultiScaleArrayLeaf, i::Int)
    deleteat!(m.nodes, i)
    1
end

function remove_daughter!(m::AbstractMultiScaleArrayHead, i, I::Int...)
    del_length = __remove_daughter!(m.nodes[i], I...)
    for j = i:num_daughters(m)
        m.end_idxs[j] -= del_length
    end
    isempty(m.values) || (m.end_idxs[end] -= del_length)
    if size(m.nodes[i].nodes) == (0, )
        deleteat!(m.nodes, i)
        deleteat!(m.end_idxs, i)
    end
    nothing
end

function __remove_daughter!(m::AbstractMultiScaleArray, i, I::Int...)
    del_length = __remove_daughter!(m.nodes[i], I...)
    for j = i:num_daughters(m)
        m.end_idxs[j] -= del_length
    end
    isempty(m.values) || (m.end_idxs[end] -= del_length)
    isempty(m.nodes[i]) && deleteat!(m.nodes, i)
    del_length
end

remove_daughter!(m::AbstractMultiScaleArrayHead, i::Int) =
    (__remove_daughter!(m, i); nothing)
