level_iter(S, n::Int) = n == 1 ? S.nodes : chain((level_iter(node, n-1) for node in S.nodes)...)

struct LevelIterIdx{T}
    iter::T
end
LevelIterIdx(S::AbstractMultiScaleArray, n::Int) = LevelIterIdx(level_iter(S, n))

start(l::LevelIterIdx) = (start(l.iter), 1)
function next(l::LevelIterIdx, state)
    val, new_state = next(l.iter, state[1])
    end_idx = state[2] + length(val) - 1
    ((val, state[2], end_idx), (new_state, end_idx + 1))
end
done(l::LevelIterIdx, state) = done(l.iter, state[1])

function Base.iterate(l::LevelIterIdx)
    x = iterate(l.iter)
    x == nothing && return nothing
    val, new_state = x
    end_idx = 1 + length(val) - 1
    ((val, 1, end_idx), (new_state, end_idx + 1))
end

function Base.iterate(l::LevelIterIdx,state)
    x = iterate(l.iter,state[1])
    x == nothing && return nothing
    val, new_state = x
    end_idx = state[2] + length(val) - 1
    ((val, state[2], end_idx), (new_state, end_idx + 1))
end
