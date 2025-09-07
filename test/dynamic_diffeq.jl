using MultiScaleArrays
using OrdinaryDiffEq, DiffEqBase, Test, StochasticDiffEq, Statistics

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

cell1 = Cell([1.0; 2.0; 3.0])
cell2 = Cell([4.0; 5])
population = construct(Population, deepcopy([cell1, cell2])) # Make a Population from cells
cell3 = Cell([3.0; 2.0; 5.0])
cell4 = Cell([4.0; 6])
population2 = construct(Population, deepcopy([cell3, cell4]))
tissue1 = construct(Tissue, deepcopy([population, population2])) # Make a Tissue from Populations
tissue2 = construct(Tissue, deepcopy([population2, population]))
_embryo = construct(Embryo, deepcopy([tissue1, tissue2])) # Make an embryo from Tissues
embryo = deepcopy(_embryo)

cell_ode = function (dcell, cell, p, t)
    m = mean(cell)
    for i in eachindex(cell)
        dcell[i] = -m * cell[i]
    end
end

f = function (dembryo, embryo, p, t)
    for (cell, dcell) in LevelIter(3, embryo, dembryo)
        cell_ode(dcell, cell, p, t)
    end
end

tstop = [0.5]

condition = function (u, t, integrator)
    t ∈ tstop
end

affect! = function (integrator)
    add_node!(integrator, integrator.u[1, 1, 1], 1, 1)
end

growing_cb = DiscreteCallback(condition, affect!)

println("Do the ODE Part")

prob = ODEProblem(f, embryo, (0.0, 1.0))
test_embryo = deepcopy(embryo)

sol = solve(prob, Tsit5(), callback = growing_cb, tstops = tstop)
sol = solve(prob, Rosenbrock23(autodiff = false), tstops = tstop)
sol = solve(prob, Rosenbrock23(autodiff = false), callback = growing_cb, tstops = tstop)
sol = solve(prob, Rosenbrock23(chunk_size = 1), callback = growing_cb, tstops = tstop)

affect_del! = function (integrator)
    remove_node!(integrator, 1, 1, 1)
end

shrinking_cb = DiscreteCallback(condition, affect_del!)

sol = solve(prob, Tsit5(), callback = shrinking_cb, tstops = tstop)

sol = solve(prob, Rosenbrock23(autodiff = false), callback = shrinking_cb, tstops = tstop)

sol = solve(prob, Rosenbrock23(chunk_size = 1), callback = shrinking_cb, tstops = tstop)

@test length(sol[end]) == 17

println("Do the SDE Part")

g = function (du, u, p, t)
    for i in eachindex(u)
        du[i] = 0.1u[i]
    end
end
prob = SDEProblem(f, g, embryo, (0.0, 1.0))

@show SRIW1

sol = solve(prob, SRIW1(), callback = growing_cb, tstops = tstop)

@show SRA1

sol = solve(prob, SRA1(), callback = growing_cb, tstops = tstop)

@show RKMil

sol = solve(prob, RKMil(), callback = growing_cb, dt = 1 / 10, tstops = tstop)

@show EM

sol = solve(prob, EM(), dt = 1 / 20, callback = growing_cb, tstops = tstop)

@test length(sol[end]) == 23

@show SRIW1

sol = solve(prob, SRIW1(), callback = shrinking_cb, tstops = tstop)

@show SRA1

sol = solve(prob, SRA1(), callback = shrinking_cb, tstops = tstop)

@show RKMil

sol = solve(prob, RKMil(), dt = 1 / 10, callback = shrinking_cb, tstops = tstop)

@show EM

sol = solve(prob, EM(), dt = 1 / 10, callback = shrinking_cb, tstops = tstop)

@test length(sol[end]) == 17
