function bisect_search(a,i)
  L=1; R=length(a)
  while true
    L>R && error("Left > Right in Bisection. Failed")
    m = floor(Int,.5*(L+R))
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
  #length(a) - sum(map((x)->x>=i,a))
end

linearindexing(m::AbstractMultiScaleModel) = Base.LinearFast()

function getindex(m::AbstractMultiScaleModel,i::Int)
  idx = bisect_search(m.end_idxs,i)
  if idx > 1
    i = i-m.end_idxs[idx-1]
  end
  m.x[idx][i]
end

function setindex!(m::AbstractMultiScaleModel,x,i::Int)
  idx = bisect_search(m.end_idxs,i) # +1 for 1-based indexing
  if idx > 1
    i = i-m.end_idxs[idx-1]
  end
  m.x[idx][i] = x
end

function getindex(m::MultiScaleModelLeaf,i::Int)
  m.x[i]
end

function getindex(m::MultiScaleModelLeaf,i::Int...)
  m.x[i[1]]
end

function getindex(m::AbstractMultiScaleModel,i...)
  @show i
  m.x[i[1]][i[2:end]...]
end

function getindex(m::MultiScaleModelLeaf,i...)
  m.x[i[1]]
end

function getindex(m::MultiScaleModelLeaf,i::CartesianIndex{1})
  m.x[i[1]]
end

function setindex!(m::MultiScaleModelLeaf,x,i::Int)
  m.x[i] = x
end

function setindex!(m::MultiScaleModelLeaf,x,i::Int...)
  m.x[i[1]] = x
end

function setindex!(m::AbstractMultiScaleModel,x,i::Int...)
  m.x[i[1]][i[2:end]...] = x
end

function getindex(m::AbstractMultiScaleModel,::Colon)
  [m[i] for i in 1:length(m)]
end

function getindex(m::AbstractMultiScaleModel,i::CartesianIndex{1}) # (i,)
  m[i[1]]
end

eachindex(m::AbstractMultiScaleModel) = 1:length(m)
endof(m::AbstractMultiScaleModel) = length(m)

Base.start(::AbstractMultiScaleModel) = 1
Base.next(S::AbstractMultiScaleModel, state) = (S[state], state+1)
Base.done(S::AbstractMultiScaleModel, state) = state > length(S);

#broadcast_getindex(m::MultiScaleModelLeaf,i::Int)    =  (println("here");m[i])
#broadcast_getindex(m::AbstractMultiScaleModel,i::Int)    =  (println("here");m[i])
#broadcast_getindex(m::AbstractMultiScaleModel,i::Int...) = m[i]
