function __add_node!(m::AbstractMultiScaleArray, node::AbstractMultiScaleArray)
    push!(m.end_idxs, m.end_idxs[end] + length(node) + length(m.values))
    push!(m.nodes, node)
    return nothing
end

function __update_lengths(m::AbstractMultiScaleArray, beg, len)
    for j in beg:num_nodes(m)
        m.end_idxs[j] += len
    end
    isempty(m.values) || (m.end_idxs[end] += len)
    return nothing
end

function __add_node!(m::AbstractMultiScaleArray, node::AbstractMultiScaleArray, i::Int)
    __add_node!(m.nodes[i], node)
    return __update_lengths(m, i, length(node))
end

function __add_node!(m::AbstractMultiScaleArray, node::AbstractMultiScaleArray, i, I::Int...)
    __add_node!(m.nodes[i], node, I...)
    return __update_lengths(m, i, length(node))
end

add_node!(m::AbstractMultiScaleArrayHead, node::AbstractMultiScaleArray) = __add_node!(m, node)

add_node!(m::AbstractMultiScaleArrayHead, node::AbstractMultiScaleArray, i::Int) = __add_node!(m,
                                                                                               node,
                                                                                               i)

add_node!(m::AbstractMultiScaleArrayHead, node::AbstractMultiScaleArray, i, I::Int...) = __add_node!(m,
                                                                                                     node,
                                                                                                     i,
                                                                                                     I...)

function __remove_node!(m::AbstractMultiScaleArray, i::Int)
    del_length = length(m.nodes[i])
    deleteat!(m.nodes, i)
    deleteat!(m.end_idxs, i)
    for j in i:num_nodes(m)
        m.end_idxs[j] -= del_length
    end
    isempty(m.values) || (m.end_idxs[end] -= del_length)
    return del_length
end

function __remove_node!(m::AbstractMultiScaleArrayLeaf, i::Int)
    deleteat!(m.nodes, i)
    return 1
end

function remove_node!(m::AbstractMultiScaleArrayHead, i, I::Int...)
    del_length = __remove_node!(m.nodes[i], I...)
    for j in i:num_nodes(m)
        m.end_idxs[j] -= del_length
    end
    isempty(m.values) || (m.end_idxs[end] -= del_length)
    if size(m.nodes[i].nodes) == (0,)
        deleteat!(m.nodes, i)
        deleteat!(m.end_idxs, i)
    end
    return nothing
end

function __remove_node!(m::AbstractMultiScaleArray, i, I::Int...)
    del_length = __remove_node!(m.nodes[i], I...)
    for j in i:num_nodes(m)
        m.end_idxs[j] -= del_length
    end
    isempty(m.values) || (m.end_idxs[end] -= del_length)
    isempty(m.nodes[i]) && deleteat!(m.nodes, i)
    return del_length
end

remove_node!(m::AbstractMultiScaleArrayHead, i::Int) = (__remove_node!(m, i); nothing)
