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

function __add_node!(
        m::AbstractMultiScaleArray, node::AbstractMultiScaleArray, i,
        I::Int...
    )
    __add_node!(m.nodes[i], node, I...)
    return __update_lengths(m, i, length(node))
end

"""
    add_node!(m::AbstractMultiScaleArrayHead, node::AbstractMultiScaleArray, I...)

Insert `node` at the location identified by `I` and update cached linear indices.

# Arguments

- `m`: An [`AbstractMultiScaleArrayHead`](@ref) to mutate.
- `node`: The multiscale array node to append or insert.
- `I...`: Optional indices that select the parent node. With no indices, `node` is appended to
  `m`; with one or more indices, it is inserted into the selected descendant.

# Examples

```julia
add_node!(model, Cell([0.0]))
add_node!(model, Cell([0.0]), 2)
```
"""
function add_node!(m::AbstractMultiScaleArrayHead, node::AbstractMultiScaleArray)
    return __add_node!(m, node)
end

function add_node!(m::AbstractMultiScaleArrayHead, node::AbstractMultiScaleArray, i::Int)
    return __add_node!(m, node, i)
end

function add_node!(
        m::AbstractMultiScaleArrayHead, node::AbstractMultiScaleArray, i,
        I::Int...
    )
    return __add_node!(m, node, i, I...)
end

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

"""
    remove_node!(m::AbstractMultiScaleArrayHead, I...)

Remove the node at the location identified by `I` and update cached linear indices.

# Arguments

- `m`: An [`AbstractMultiScaleArrayHead`](@ref) to mutate.
- `I...`: Indices selecting the node to remove. A single index removes a direct child; additional
  indices descend through the hierarchy.

# Examples

```julia
remove_node!(model, 2)
remove_node!(model, 2, 1)
```
"""
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
