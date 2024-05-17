using OrdinaryDiffEq, Test
using MultiScaleArrays

struct CellGenotype{B <: Float64} <: AbstractMultiScaleArrayLeaf{B}
    values::Vector{B}
    Iam::Symbol
end

struct Host{T <: AbstractMultiScaleArray, B <: Float64} <: AbstractMultiScaleArrayHead{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
end

host = construct(Host,
    [
        CellGenotype([Float64(1)],
            :aSymbiontGenotype),
        CellGenotype([Float64(1)],
            :aSymbiontGenotype)])

evolve = function (du, u::Host, p, t)
    for i in 1:length(u)
        du[i] = 2 * u[i]
    end
end

similar(host)
similar(host.nodes[1])
similar(host.nodes[1], Float64)
similar(host, Float64)

## Solve
tspan = (0.0, 50)
prob = ODEProblem(evolve, host, tspan)
sol = solve(prob, Tsit5())

struct Host2{T <: AbstractMultiScaleArray, B <: Float64} <: AbstractMultiScaleArrayHead{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
    Iam::Symbol
end

host2 = construct(Host2,
    [
        CellGenotype([Float64(1)],
            :aSymbiontGenotype),
        CellGenotype([Float64(1)],
            :aSymbiontGenotype)], Float64[],
    :my_cell_type)

evolve = function (du, u::Host2, p, t)
    for i in 1:length(u)
        du[i] = 2 * u[i]
    end
end

similar(host2)
similar(host2.nodes[1])
similar(host2.nodes[1], Float64)
similar(host2, Float64)

## Solve
tspan = (0.0, 50)
prob = ODEProblem(evolve, host2, tspan)
sol = solve(prob, Tsit5())
