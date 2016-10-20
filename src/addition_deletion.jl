function add_daughter!(m::MultiScaleModelHead,x::AbstractMultiScaleModel)
  push!(m.x,x)
  push!(m.end_idxs,m.end_idxs[end]+length(x))
  nothing
end

function __add_daughter!(m::AbstractMultiScaleModel,x::AbstractMultiScaleModel,i::Int...)
  __add_daughter!(m.x[i[1]],x,i[2:end])
  for j = i[1]:num_daughters(m)
    m.end_idxs[j]  += length(x)
  end
  nothing
end

function __add_daughter!(m::AbstractMultiScaleModel,x::AbstractMultiScaleModel)
  push!(m.x,x)
  push!(m.end_idxs,m.end_idxs[end]+length(x))
  nothing
end

function __add_daughter!(m::AbstractMultiScaleModel,x::AbstractMultiScaleModel,i::Tuple{Int64})
  __add_daughter!(m,x,i[1])
end

function __add_daughter!(m::AbstractMultiScaleModel,x::AbstractMultiScaleModel,i::Int64)
  __add_daughter!(m.x[i],x)
  for j = i:num_daughters(m)
    m.end_idxs[j]  += length(x)
  end
  nothing
end

function add_daughter!(m::MultiScaleModelHead,x::AbstractMultiScaleModel,i::Int)
  __add_daughter!(m.x[i],x)
  for j = i:num_daughters(m)
    m.end_idxs[j]  += length(x)
  end
  nothing
end

function add_daughter!(m::MultiScaleModelHead,x::AbstractMultiScaleModel,i::Int...)
  __add_daughter!(m.x[i[1]],x,i[2:end])
  for j = i[1]:num_daughters(m)
    m.end_idxs[j]  += length(x)
  end
  nothing
end

function __remove_daughter!(m::AbstractMultiScaleModel,i::Int)
  del_length = length(m.x[i])
  deleteat!(m.x,i); deleteat!(m.end_idxs,i)
  for j = i:num_daughters(m)
    m.end_idxs[j] -= del_length
  end
  del_length
end

function remove_daughter!(m::MultiScaleModelHead,i::Int)
  del_length = length(m.x[i])
  deleteat!(m.x,i); deleteat!(m.end_idxs,i)
  for j = i:num_daughters(m)
    m.end_idxs[j] -= del_length
  end
  nothing
end

function __remove_daughter!(m::MultiScaleModelLeaf,i::Int)
  deleteat!(m.x,i)
  1
end

function __remove_daughter!(m::AbstractMultiScaleModel,i::Tuple{Int64})
  __remove_daughter!(m,i[1])
end

function remove_daughter!(m::MultiScaleModelHead,i::Int...)
  del_length = __remove_daughter!(m.x[i[1]],i[2:end])
  for j = i[1]:num_daughters(m)
    m.end_idxs[j] -= del_length
  end
  if size(m.x[i[1]].x) == (0,)
    deleteat!(m.x,i[1])
    deleteat!(m.end_idxs,i[1])
  end
  nothing
end

function __remove_daughter!(m::AbstractMultiScaleModel,i::Int...)
  del_length = __remove_daughter!(m.x[i[1]],i[2:end])
  for j = i[1]:num_daughters(m)
    m.end_idxs[j] -= del_length
  end
  if length(m.x[i[1]]) == 0
    deleteat!(m.x,i[1])
  end
  del_length
end
