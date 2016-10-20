function bisect_search(a,i)
  length(a) - sum(map((x)->x>=i,a))

  ## My bisection sucks... please fix it.

  #=
  ### Bisection serach for the index closet index s.t. a[idx]<i
  if a[end]<=i
    return length(a)
  end
  L = 0; R = length(a)-1; M = 0
  while true
    M = (L+R)÷2
    println("$L $M $R")
    if R-L <= 1
      break
    end
    if a[M] > i
      R = M
    elseif a[M] < i
      L = M
    end
  end
  if R == 1
    return a[1]>i
  end
  if L == R
    return L
  end
  println(a[M])
  if a[M] == R
    idx = R
  else
    idx = L
  end
  idx
  =#
end

linearindexing(m::AbstractMultiScaleModel) = Base.LinearFast()

# In general, we simply re-index the parent indices by the provided ones
typealias SlowSubArray{T,N,P,I} SubArray{T,N,P,I,false}
@inline function getindex{T<:AbstractMultiScaleModel}(V::SlowSubArray{T,1}, I::Vararg{Real,1})
    @boundscheck checkbounds(V, I...)
    idxs = Base.reindex(V, V.indexes, Base.to_indexes(I...))
    @show I
    @show typeof(I)
    @show idxs[1]
    @show V.parent
    @show typeof(idxs[1])
    @show @which V.parent[idxs[1]]
    @inbounds r = V.parent[idxs[1]]
    r
end
@inline function getindex{T<:AbstractMultiScaleModel,N}(V::SlowSubArray{T,N}, I::Tuple{Int64})
    @boundscheck checkbounds(V, I...)
    idxs = Base.reindex(V, V.indexes, Base.to_indexes(I...))
    @inbounds r = V.parent[idxs[1]]
    r
end

function getindex(m::AbstractMultiScaleModel,i::Int)
  idx = bisect_search(m.end_idxs,i)+1 # +1 for 1-based indexing
  if idx > 1
    i = i-m.end_idxs[idx-1]
  end
  m.x[idx][i]
end

function setindex!(m::AbstractMultiScaleModel,x,i::Int)
  idx = bisect_search(m.end_idxs,i)+1 # +1 for 1-based indexing
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
  m.x[i[1]][i[2:end]...]
end

function getindex(m::MultiScaleModelLeaf,i...)
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

eachindex(m::AbstractMultiScaleModel) = 1:length(m)

#broadcast_getindex(m::MultiScaleModelLeaf,i::Int)    =  (println("here");m[i])
#broadcast_getindex(m::AbstractMultiScaleModel,i::Int)    =  (println("here");m[i])
#broadcast_getindex(m::AbstractMultiScaleModel,i::Int...) = m[i]
