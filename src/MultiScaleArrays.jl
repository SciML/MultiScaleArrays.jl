__precompile__()

module MultiScaleArrays

import Base: length, push!, deleteat!, getindex, setindex!, eachindex, ndims,
       size, print_matrix, similar, hcat, vcat,
       ==, *, +, /, -, show, vec, reshape

import RecursiveArrayTools: recursivecopy,recursivecopy!

import RecursiveArrayTools: chain

abstract type AbstractMultiScaleArray{B}     <: AbstractVector{B} end
abstract type AbstractMultiScaleArrayLeaf{B} <: AbstractMultiScaleArray{B} end
abstract type AbstractMultiScaleArrayHead{B} <: AbstractMultiScaleArray{B} end

using DiffEqBase, Statistics, LinearAlgebra, FiniteDiff
import OrdinaryDiffEq, StochasticDiffEq, ForwardDiff, SparseDiffTools

include("shape_construction.jl")
include("addition_deletion.jl")
include("indexing.jl")
include("math.jl")
include("level_iterations.jl")
include("diffeq.jl")
include("show_MultiScaleArray.jl")

# Types
export AbstractMultiScaleArray, AbstractMultiScaleArrayLeaf,
       AbstractMultiScaleArrayHead

# Constructors
export construct, recursivecopy!

# Addition Deletion
export add_node!, remove_node!

# Indexing
export num_nodes, getindices

# Misc
export LevelIterIdx, LevelIter, level_iter

end # module
