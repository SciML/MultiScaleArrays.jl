using MultiScaleArrays, Test

struct Cell{B} <: AbstractMultiScaleArrayLeaf{B}
    values::Vector{B}
end
struct Population{T <: AbstractMultiScaleArray, B <: Number} <: AbstractMultiScaleArray{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
end
struct Tissue{T <: AbstractMultiScaleArray, B <: Number} <: AbstractMultiScaleArray{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
end
struct Embryo{T <: AbstractMultiScaleArray, B <: Number} <: AbstractMultiScaleArrayHead{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
end

cell1 = Cell([1.0, 2.0, 3.0])
cell2 = Cell([4.0, 5.0])
cell3 = Cell([3.0, 2.0, 5.0])
cell4 = Cell([4.0, 6.0])
p1 = construct(Population, deepcopy([cell1, cell2]))
p2 = construct(Population, deepcopy([cell3, cell4]))
tis1 = construct(Tissue, deepcopy([p1, p2]))
tis2 = construct(Tissue, deepcopy([p2, p1]))
em = construct(Embryo, deepcopy([tis1, tis2]))

@test level_iter(em, 1) == em.nodes

pops = collect(level_iter(em, 2))
@test length(pops) == 4
@test pops[1] == em.nodes[1].nodes[1]
@test pops[2] == em.nodes[1].nodes[2]
@test pops[3] == em.nodes[2].nodes[1]
@test pops[4] == em.nodes[2].nodes[2]

cells = collect(level_iter(em, 3))
@test length(cells) == 8
@test cells[1] == em.nodes[1].nodes[1].nodes[1]
@test cells[end] == em.nodes[2].nodes[2].nodes[2]
