module MultiScaleModels

import Base: length, push!, deleteat!,getindex, setindex!, eachindex,
       ndims, size, print_matrix, similar, broadcast_getindex, hcat, vcat, ==,
       linearindexing, .*, .+, *, +,/,./,-,.-,show
import RecursiveArrayTools: recursivecopy!
abstract AbstractMultiScaleModel{T} <: AbstractArray{T,1}
abstract MultiScaleModelLeaf{T} <: AbstractMultiScaleModel{T}
abstract MultiScaleModelHead{T} <: AbstractMultiScaleModel{T}

Base.show(io::IO, x::AbstractMultiScaleModel) = invoke(show, Tuple{IO, Any}, io, x)
Base.show(io::IO, ::MIME"text/plain", x::AbstractMultiScaleModel) = show(io, x)

include("shape_construction.jl")
include("addition_deletion.jl")
include("indexing.jl")
include("math.jl")

# Types
export AbstractMultiScaleModel, MultiScaleModelLeaf, MultiScaleModelHead

# Constructors
export construct, similar, deepcopy, recursivecopy!

# Addition Deletion
export add_daughter!, remove_daughter!

# Indexing
export getindex, setindex!, eachindex, length, num_daughters,
       ndims, size, broadcast_getindex, hcat, vcat, linearindexing

# Math and Logic
export ==, .*, .+, *, +,/,./,-,.-

# Misc
export print_matrix

end # module
