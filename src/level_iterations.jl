function level_iter(S, n::Int)
    n == 1 ? S.nodes : chain((level_iter(node, n - 1) for node in S.nodes)...)
end

struct LevelIterIdx{T}
    iter::T
end
LevelIterIdx(S::AbstractMultiScaleArray, n::Int) = LevelIterIdx(level_iter(S, n))

function Base.iterate(l::LevelIterIdx)
    x = iterate(l.iter)
    x == nothing && return nothing
    val, new_state = x
    end_idx = 1 + length(val) - 1
    ((val, 1, end_idx), (new_state, end_idx + 1))
end

function Base.iterate(l::LevelIterIdx, state)
    x = iterate(l.iter, state[1])
    x == nothing && return nothing
    val, new_state = x
    end_idx = state[2] + length(val) - 1
    ((val, state[2], end_idx), (new_state, end_idx + 1))
end

LevelIter(n::Int, S::AbstractMultiScaleArray...) = zip((level_iter(s, n) for s in S)...)
