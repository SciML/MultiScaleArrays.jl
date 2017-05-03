__precompile__()

module MultiScaleArrays

import Base: length, push!, deleteat!,getindex, setindex!, eachindex,
       ndims, size, print_matrix, similar, broadcast_getindex, hcat, vcat, ==,
       linearindexing, .*, .+, *, +,/,./,-,.-,show, vec, reshape

using Compat
import RecursiveArrayTools: recursivecopy!
using Iterators
@compat abstract type AbstractMultiScaleArray{B} <: AbstractArray{B,1} end
@compat abstract type AbstractMultiScaleArrayLeaf{B} <: AbstractMultiScaleArray{B} end
@compat abstract type AbstractMultiScaleArrayHead{B} <: AbstractMultiScaleArray{B} end

using DiffEqBase

Base.show(io::IO, x::AbstractMultiScaleArray) = invoke(show, Tuple{IO, Any}, io, x)
Base.show(io::IO, ::MIME"text/plain", x::AbstractMultiScaleArray) = show(io, x)

include("utils.jl")
include("shape_construction.jl")
include("addition_deletion.jl")
include("indexing.jl")
include("math.jl")
include("level_iterations.jl")
include("diffeq.jl")

# Types
export AbstractMultiScaleArray, AbstractMultiScaleArrayLeaf, AbstractMultiScaleArrayHead

# Constructors
export construct, similar, deepcopy, recursivecopy!

# Addition Deletion
export add_daughter!, remove_daughter!

# Indexing
export getindex, setindex!, eachindex, length, num_daughters, getindices,
       ndims, size, broadcast_getindex, hcat, vcat, linearindexing

# Math and Logic
export ==, .*, .+, *, +,/,./,-,.-

# Misc
export print_matrix

export LevelIterIdx, level_iter

end # module
