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
population = construct(Population, deepcopy([cell1, cell2]), [11.0, 12.0, 13.0]) # Make a Population from cells
cell3 = Cell([3.0; 2.0; 5.0])
cell4 = Cell([4.0; 6])
population2 = construct(Population, deepcopy([cell3, cell4]), [8.0, 9.0])

@test population[1] == 1.0
@test population[2] == 2.0
@test population[3] == 3.0
@test population[4] == 4.0
@test population[5] == 5.0
@test population[6] == 11.0
@test population[7] == 12.0
@test population[8] == 13.0

tissue1 = construct(Tissue, deepcopy([population, population2])) # Make a Tissue from Populations
tissue2 = construct(Tissue, deepcopy([population2, population]))
embryo = construct(Embryo, deepcopy([tissue1, tissue2]), [31.0, 32.0]) # Make an embryo from Tissues

@test tissue1[1] == 1.0
@test tissue1[2] == 2.0
@test tissue1[3] == 3.0
@test tissue1[4] == 4.0
@test tissue1[5] == 5.0
@test tissue1[6] == 11.0
@test tissue1[7] == 12.0
@test tissue1[8] == 13.0

@test tissue1[9] == 3.0
@test tissue1[10] == 2.0
@test tissue1[11] == 5.0
@test tissue1[12] == 4.0
@test tissue1[13] == 6.0
@test tissue1[14] == 8.0
@test tissue1[15] == 9.0

@test embryo[1] == 1.0
@test embryo[2] == 2.0
@test embryo[3] == 3.0
@test embryo[4] == 4.0
@test embryo[5] == 5.0
@test embryo[6] == 11.0
@test embryo[7] == 12.0
@test embryo[8] == 13.0

@test embryo[9] == 3.0
@test embryo[10] == 2.0
@test embryo[11] == 5.0
@test embryo[12] == 4.0
@test embryo[13] == 6.0
@test embryo[14] == 8.0
@test embryo[15] == 9.0

@test embryo[16] == 3.0
@test embryo[17] == 2.0
@test embryo[18] == 5.0
@test embryo[19] == 4.0
@test embryo[20] == 6.0
@test embryo[21] == 8.0
@test embryo[22] == 9.0

@test embryo[23] == 1.0
@test embryo[24] == 2.0
@test embryo[25] == 3.0
@test embryo[26] == 4.0
@test embryo[27] == 5.0
@test embryo[28] == 11.0
@test embryo[29] == 12.0
@test embryo[30] == 13.0

@test embryo[31] == 31.0
@test embryo[32] == 32.0
