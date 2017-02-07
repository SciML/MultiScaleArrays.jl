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

linearindexing{T<:AbstractMultiScaleArray}(::Type{T}) = Base.LinearFast()

@inline function getindex(m::AbstractMultiScaleArray,i::Int)
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

@inline function setindex!(m::AbstractMultiScaleArray,x,i::Int)
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

@inline function getindex(m::AbstractMultiScaleArrayLeaf,i::Int)
  m.x[i]
end

@inline function getindex(m::AbstractMultiScaleArrayLeaf,i::Int...)
  m.x[i[1]]
end

@inline function getindex(m::AbstractMultiScaleArray,i...)
  if isempty(m.y) || i[1] < length(m.end_idxs)
    m.x[i[1]][i[2:end]...]
  else
    m.y[i[2:end]...]
  end
end

@inline function getindex(m::AbstractMultiScaleArrayLeaf,i...)
  m.x[i[1]]
end

@inline function getindex(m::AbstractMultiScaleArrayLeaf,i::CartesianIndex{1})
  m.x[i[1]]
end

@inline function setindex!(m::AbstractMultiScaleArrayLeaf,x,i::Int)
  m.x[i] = x
end

@inline function setindex!(m::AbstractMultiScaleArrayLeaf,x,i::Int...)
  m.x[i[1]] = x
end

@inline function setindex!(m::AbstractMultiScaleArray,x,i::Int...)
  if isempty(m.y) || i[1] < length(m.end_idxs)
    m.x[i[1]][i[2:end]...] = x
  else
    m.y[i[2:end]...] = x
  end
end

@inline function getindex(m::AbstractMultiScaleArray,::Colon)
  [m[i] for i in 1:length(m)]
end

@inline function getindex(m::AbstractMultiScaleArrayLeaf,::Colon)
  m.x
end

@inline function getindex(m::AbstractMultiScaleArray,i::CartesianIndex{1}) # (i,)
  m[i[1]]
end

eachindex(m::AbstractMultiScaleArray) = 1:length(m)
endof(m::AbstractMultiScaleArray) = length(m)

Base.eltype{B}(S::AbstractMultiScaleArray{B}) = B
#broadcast_getindex(m::AbstractMultiScaleArrayLeaf,i::Int)    =  (println("here");m[i])
#broadcast_getindex(m::AbstractMultiScaleArray,i::Int)    =  (println("here");m[i])
#broadcast_getindex(m::AbstractMultiScaleArray,i::Int...) = m[i]
