"""
    level_iter(S::AbstractMultiScaleArray, n::Int)

Iterate over nodes `n` levels below the top of `S`.

# Arguments

- `S`: A multiscale array node whose descendants are traversed.
- `n`: Number of levels below `S`. `n = 1` iterates over direct children.

# Examples

```julia
for cell in level_iter(model, 2)
    println(length(cell))
end
```
"""
function level_iter(S, n::Int)
    return n == 1 ? S.nodes : RecursiveArrayTools.chain((level_iter(node, n - 1) for node in S.nodes)...)
end

"""
    LevelIterIdx(S::AbstractMultiScaleArray, n::Int)

Iterator over nodes `n` levels below `S`, yielding each node and its linear index range.

# Arguments

- `S`: A multiscale array node whose descendants are traversed.
- `n`: Number of levels below `S` to visit.

# Examples

```julia
for node, first_index, last_index in LevelIterIdx(model, 2)
    @view model[first_index:last_index]
end
```
"""
struct LevelIterIdx{T}
    iter::T
end
LevelIterIdx(S::AbstractMultiScaleArray, n::Int) = LevelIterIdx(level_iter(S, n))

function Base.iterate(l::LevelIterIdx)
    x = iterate(l.iter)
    x == nothing && return nothing
    val, new_state = x
    end_idx = 1 + length(val) - 1
    return ((val, 1, end_idx), (new_state, end_idx + 1))
end

function Base.iterate(l::LevelIterIdx, state)
    x = iterate(l.iter, state[1])
    x == nothing && return nothing
    val, new_state = x
    end_idx = state[2] + length(val) - 1
    return ((val, state[2], end_idx), (new_state, end_idx + 1))
end

"""
    LevelIter(n::Int, S::AbstractMultiScaleArray...)

Zip `level_iter(s, n)` across one or more multiscale arrays.

# Arguments

- `n`: Number of levels below each input to visit.
- `S...`: One or more multiscale arrays with matching hierarchy structure.

# Examples

```julia
for (state_node, derivative_node) in LevelIter(2, state, derivative)
    derivative_node .= state_node
end
```
"""
LevelIter(n::Int, S::AbstractMultiScaleArray...) = zip((level_iter(s, n) for s in S)...)
