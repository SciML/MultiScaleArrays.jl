using OrdinaryDiffEq, MultiScaleArrays

# Create Structures that define a Cell and Population
struct Cell{B} <: AbstractMultiScaleArrayLeaf{B}
    values::Vector{B}
end

struct Population{T <: AbstractMultiScaleArray, B <: Number} <:
       AbstractMultiScaleArrayHead{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
end

nCells = 2
cellConstructor = [Cell([1.0, 2.0]) for i in 1:nCells]
cellPop = construct(Population, deepcopy(cellConstructor))

function StateModel(dpop, pop, p, t) #This is arbitrary
    for (cell, dcell) in LevelIter(1, pop, dpop)
        dcell = 1.0
    end
end

prob = ODEProblem(StateModel, cellPop, (0.0, 1.0))

integrator = init(prob, Tsit5())

add_node!(integrator, Cell([5.0, 20.0]))

@test integrator.u.nodes[end].values == [5.0, 20.0]
