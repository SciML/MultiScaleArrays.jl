using MultiScaleArrays
using OrdinaryDiffEq, DiffEqBase, Test, StochasticDiffEq, Statistics

struct Cell2{B} <: AbstractMultiScaleArrayLeaf{B}
    values::Vector{B}
end
struct Population2{T <: AbstractMultiScaleArray, B <: Number} <:
       AbstractMultiScaleArrayHead{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
end

cell1 = Cell2([1.0; 2.0; 3.0])
cell2 = Cell2([4.0; 5])
cell3 = Cell2([3.0; 2.0; 5.0])
cell4 = Cell2([4.0; 6])
pop = construct(Population2, deepcopy([cell1, cell2, cell3, cell4]))

function cell_ode4(dcell, cell, p, t)
    m = mean(cell)
    for i in eachindex(cell)
        dcell[i] = -m * cell[i]
    end
end

function f4(dembryo, embryo, p, t)
    for (cell, dcell) in LevelIter(1, embryo, dembryo)
        cell_ode4(dcell, cell, p, t)
    end
end

tstop = [0.5]

condition = function (u, t, integrator)
    t ∈ tstop
end

affect! = function (integrator)
    add_node!(integrator, integrator.u.nodes[1])
end

growing_cb = DiscreteCallback(condition, affect!)

prob = ODEProblem(f4, deepcopy(pop), (0.0, 1.0))

add_node!(pop, pop.nodes[1])

sol = solve(prob, Tsit5(), callback = growing_cb, tstops = tstop)

sol = solve(prob, Rosenbrock23(chunk_size = 1), callback = growing_cb, tstops = tstop)

@test length(sol[end]) == 13

affect_del! = function (integrator)
    remove_node!(integrator, 1)
end

shrinking_cb = DiscreteCallback(condition, affect_del!)

prob = ODEProblem(f4, deepcopy(pop), (0.0, 1.0))
sol = solve(prob, Tsit5(), callback = shrinking_cb, tstops = tstop)

prob = ODEProblem(f4, deepcopy(pop), (0.0, 1.0))
sol = solve(prob, Rosenbrock23(chunk_size = 1), callback = shrinking_cb, tstops = tstop)
@test length(sol[end]) == 10

println("Do the SDE Part")

function g4(du, u, p, t)
    for i in eachindex(u)
        du[i] = 0.1u[i]
    end
end

prob = SDEProblem(f4, g4, deepcopy(pop), (0.0, 1.0))

sol = solve(prob, SOSRI(), callback = growing_cb, tstops = tstop)

@test length(sol[end]) == 16

sol = solve(prob, SOSRI(), callback = shrinking_cb, tstops = tstop)

@test length(sol[end]) == 10
