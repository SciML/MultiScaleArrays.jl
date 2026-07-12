"""
    level_iter(S::AbstractMultiScaleArray, n::Int)

Iterate over nodes `n` levels below the top of `S`.
"""
function level_iter(S, n::Int)
    return n == 1 ? S.nodes : chain((level_iter(node, n - 1) for node in S.nodes)...)
end

"""
    LevelIterIdx(S::AbstractMultiScaleArray, n::Int)

Iterator over nodes `n` levels below `S`, yielding each node and its linear index range.
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
"""
LevelIter(n::Int, S::AbstractMultiScaleArray...) = zip((level_iter(s, n) for s in S)...)
