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
    i = i-m.end_idxs[idx-1] # also works with y
  end
  if isempty(m.y) || idx < length(m.end_idxs)
    return m.x[idx][i]
  else
    return m.y[i]
  end
end

function setindex!(m::AbstractMultiScaleModel,x,i::Int)
  idx = bisect_search(m.end_idxs,i) # +1 for 1-based indexing
  if idx > 1
    i = i-m.end_idxs[idx-1]
  end
  if isempty(m.y) || idx < length(m.end_idxs)
    m.x[idx][i] = x
  else
    m.y[i] = x
  end
end

function getindex(m::MultiScaleModelLeaf,i::Int)
  m.x[i]
end

function getindex(m::MultiScaleModelLeaf,i::Int...)
  m.x[i[1]]
end

function getindex(m::AbstractMultiScaleModel,i...)
  if isempty(m.y) || i[1] < length(m.end_idxs)
    m.x[i[1]][i[2:end]...]
  else
    m.y[i[2:end]...]
  end
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
  if isempty(m.y) || i[1] < length(m.end_idxs)
    m.x[i[1]][i[2:end]...] = x
  else
    m.y[i[2:end]...] = x
  end
end

function getindex(m::AbstractMultiScaleModel,::Colon)
  [m[i] for i in 1:length(m)]
end

function getindex(m::AbstractMultiScaleModel,i::CartesianIndex{1}) # (i,)
  m[i[1]]
end

eachindex(m::AbstractMultiScaleModel) = 1:length(m)
endof(m::AbstractMultiScaleModel) = length(m)

Base.eltype{B}(S::AbstractMultiScaleModel{B}) = B
#broadcast_getindex(m::MultiScaleModelLeaf,i::Int)    =  (println("here");m[i])
#broadcast_getindex(m::AbstractMultiScaleModel,i::Int)    =  (println("here");m[i])
#broadcast_getindex(m::AbstractMultiScaleModel,i::Int...) = m[i]
