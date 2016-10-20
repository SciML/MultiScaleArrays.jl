module MultiScaleModels

import Base: length, push!, deleteat!,getindex, setindex!, eachindex
abstract AbstractMultiScaleModel
abstract MultiScaleModelLeaf <: AbstractMultiScaleModel
abstract MultiScaleModelHead <: AbstractMultiScaleModel
length(m::MultiScaleModelLeaf) = length(m.x)
length(m::AbstractMultiScaleModel) = m.end_idxs[end]
num_daughters(m::MultiScaleModelLeaf) = 0
num_daughters(m::AbstractMultiScaleModel) = size(m.x)[1]

function construct{T<:AbstractMultiScaleModel,T2<:AbstractMultiScaleModel}(::Type{T},x::Vector{T2})
  end_idxs = Vector{Int}(length(x))
  end_idxs[1] = length(x[1])
  for i in 2:length(x)
    end_idxs[i] = end_idxs[i-1] + length(x[i])
  end
  m = T(x,end_idxs)
end

include("addition_deletion.jl")
include("indexing.jl")

# Types
export AbstractMultiScaleModel, MultiScaleModelLeaf, MultiScaleModelHead

# Constructors
export construct

# Addition Deletion
export add_daughter!, remove_daughter!

# Indexing
export getindex, setindex!, eachindex, length, num_daughters
end # module
