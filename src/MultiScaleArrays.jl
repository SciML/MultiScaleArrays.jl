module MultiScaleArrays

import Base: length, push!, deleteat!,getindex, setindex!, eachindex,
       ndims, size, print_matrix, similar, broadcast_getindex, hcat, vcat, ==,
       linearindexing, .*, .+, *, +,/,./,-,.-,show
import RecursiveArrayTools: recursivecopy!
using Iterators
abstract AbstractMultiScaleArray{B} <: AbstractArray{B,1}
abstract MultiScaleArrayLeaf{B} <: AbstractMultiScaleArray{B}
abstract MultiScaleArrayHead{B} <: AbstractMultiScaleArray{B}

Base.show(io::IO, x::AbstractMultiScaleArray) = invoke(show, Tuple{IO, Any}, io, x)
Base.show(io::IO, ::MIME"text/plain", x::AbstractMultiScaleArray) = show(io, x)

include("shape_construction.jl")
include("addition_deletion.jl")
include("indexing.jl")
include("math.jl")
include("level_iterations.jl")

# Types
export AbstractMultiScaleArray, MultiScaleArrayLeaf, MultiScaleArrayHead

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

export level_iter,LevelIterIdx, level_iter_idx

end # module
