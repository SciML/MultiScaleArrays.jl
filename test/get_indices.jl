using MultiScaleArrays
using OrdinaryDiffEq, DiffEqBase, Test, StochasticDiffEq

#=
struct Cell{B} <: AbstractMultiScaleArrayLeaf{B}
    values::Vector{B}
end
struct Population{T<:AbstractMultiScaleArray,B<:Number} <: AbstractMultiScaleArray{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
end
struct Tissue{T<:AbstractMultiScaleArray,B<:Number} <: AbstractMultiScaleArray{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
end
struct Embryo{T<:AbstractMultiScaleArray,B<:Number} <: AbstractMultiScaleArrayHead{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
end
=#

cell1 = Cell([1.0; 2.0; 3.0])
cell2 = Cell([4.0; 5])
population = construct(Population, deepcopy([cell1, cell2])) # Make a Population from cells
cell3 = Cell([3.0; 2.0; 5.0])
cell4 = Cell([4.0; 6])
population2 = construct(Population, deepcopy([cell3, cell4]))
tissue1 = construct(Tissue, deepcopy([population, population2])) # Make a Tissue from Populations
tissue2 = construct(Tissue, deepcopy([population2, population]))
embryo = construct(Embryo, deepcopy([tissue1, tissue2])) # Make an embryo from Tissues

@test getindices(embryo) == 1:20
@test getindices(embryo, 1) == 1:10
@test getindices(embryo, 2) == 11:20
@test getindices(embryo, 1, 1) == 1:5
@test getindices(embryo, 1, 2) == 6:10
@test getindices(embryo, 1, 1, 1) == 1:3
@test getindices(embryo, 1, 1, 2) == 4:5
@test getindices(embryo, 1, 2, 1) == 6:8
@test getindices(embryo, 1, 2, 2) == 9:10
@test getindices(embryo, 1, 2, 2, 1) == 9:9
@test getindices(embryo, 1, 2, 2, 2) == 10:10
@test getindices(embryo, 2, 1) == 11:15
@test getindices(embryo, 2, 2) == 16:20
@test getindices(embryo, 2, 1, 1) == 11:13
@test getindices(embryo, 2, 1, 2) == 14:15
@test getindices(embryo, 2, 2, 1) == 16:18
@test getindices(embryo, 2, 2, 2) == 19:20
@test getindices(embryo, 2, 2, 2, 1) == 19:19
@test getindices(embryo, 2, 2, 2, 2) == 20:20
