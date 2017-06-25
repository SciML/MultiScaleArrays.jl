__precompile__()

module MultiScaleArrays

import Base: length, push!, deleteat!, getindex, setindex!, eachindex, ndims, size,
       print_matrix, similar, broadcast_getindex, hcat, vcat, linearindexing,
       ==, *, +, /, -, show, vec, reshape

import RecursiveArrayTools: recursivecopy!

import RecursiveArrayTools: chain

abstract type AbstractMultiScaleArray{B}     <: AbstractVector{B} end
abstract type AbstractMultiScaleArrayLeaf{B} <: AbstractMultiScaleArray{B} end
abstract type AbstractMultiScaleArrayHead{B} <: AbstractMultiScaleArray{B} end

using DiffEqBase
using StringLiterals

Base.show(io::IO, x::AbstractMultiScaleArray) = invoke(show, Tuple{IO, Any}, io, x)
Base.show(io::IO, ::MIME"text/plain", x::AbstractMultiScaleArray) = show(io, x)

include("shape_construction.jl")
include("addition_deletion.jl")
include("indexing.jl")
include("math.jl")
include("level_iterations.jl")
include("diffeq.jl")

# Types
export AbstractMultiScaleArray, AbstractMultiScaleArrayLeaf, AbstractMultiScaleArrayHead

# Constructors
export construct, recursivecopy!

# Addition Deletion
export add_daughter!, remove_daughter!

# Indexing
export num_daughters, getindices

# Misc
export LevelIterIdx, level_iter

end # module
