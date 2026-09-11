using MultiScaleArrays, BenchmarkTools
using StableRNGs

const SUITE = BenchmarkGroup()
const rng = StableRNG(123)

# Minimal multi-scale hierarchy: Cell -> Population -> Tissue
struct Cell{B} <: AbstractMultiScaleArrayLeaf{B}
    values::Vector{B}
end
struct Population{T <: AbstractMultiScaleArray, B <: Number} <:
    AbstractMultiScaleArray{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
end
struct Tissue{T <: AbstractMultiScaleArray, B <: Number} <:
    AbstractMultiScaleArrayHead{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
end

function make_tissue(npops, ncells)
    cells = [Cell(rand(rng, 4)) for _ in 1:ncells]
    pops = [construct(Population, deepcopy(cells)) for _ in 1:npops]
    return construct(Tissue, pops)
end

tissue = make_tissue(10, 20)
tissue_big = make_tissue(50, 40)

# =============================================================================
# Construction
# =============================================================================

SUITE["construct"] = BenchmarkGroup()

cells_small = [Cell(rand(rng, 4)) for _ in 1:20]
pops_small = [construct(Population, deepcopy(cells_small)) for _ in 1:10]

SUITE["construct"]["population"] = @benchmarkable construct(
    Population, $(deepcopy(cells_small))
)
SUITE["construct"]["tissue"] = @benchmarkable construct(
    Tissue, $(deepcopy(pops_small))
)

# =============================================================================
# Traversal
# =============================================================================

SUITE["iterate"] = BenchmarkGroup()

SUITE["iterate"]["level_iter"] = @benchmarkable collect(level_iter($tissue, 1))
SUITE["iterate"]["level_iter_big"] = @benchmarkable collect(
    level_iter($tissue_big, 2)
)
SUITE["iterate"]["eachindex"] = @benchmarkable sum(x -> x, $tissue_big)

# =============================================================================
# Node mutation
# =============================================================================

SUITE["nodes"] = BenchmarkGroup()

SUITE["nodes"]["add_node!"] = @benchmarkable add_node!(t, p) setup = (
    t = make_tissue(10, 20);
    p = construct(Population, [Cell(rand($rng, 4)) for _ in 1:5])
)
SUITE["nodes"]["num_nodes"] = @benchmarkable num_nodes($tissue_big)
