function level_iter(S,n::Int)
  if n == 1
    return S.x
  else
    return chain((level_iter(x,n-1) for x in S.x)...)
  end
end

type LevelIterIdx{T}
  iter::T
end
LevelIterIdx(S::AbstractMultiScaleModel,n::Int) = LevelIterIdx(level_iter(S,n))

function Base.start(l::LevelIterIdx)
  (start(l.iter),1)
end
function Base.next(l::LevelIterIdx,state)
  val,new_state = next(l.iter,state[1])
  end_idx = state[2]+length(val)-1
  ((val,state[2],end_idx),(new_state,end_idx+1))
end
Base.done(l::LevelIterIdx,state) = done(l.iter,state[1])
