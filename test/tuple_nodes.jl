using Revise
using MultiScaleArrays, OrdinaryDiffEq, DiffEqBase, StochasticDiffEq
using Base.Test

struct PlantSettings{T} x::T end
struct OrganParams{T} y::T end

struct Organ{B<:Number,P} <: AbstractMultiScaleArrayLeaf{B}
    values::Vector{B}
    name::Symbol
    params::P
end

struct Plant{B,S,N<:Tuple{Vararg{<:Organ{<:Number}}}} <: AbstractMultiScaleArray{B}
    nodes::N
    values::Vector{B}
    end_idxs::Vector{Int}
    settings::S
end

struct Community{B,N<:Tuple{Vararg{<:Plant{<:Number}}}} <: AbstractMultiScaleArray{B}
    nodes::N
    values::Vector{B}
    end_idxs::Vector{Int}
end

mutable struct Scenario{B,N<:Tuple{Vararg{<:Community{<:Number}}}} <: AbstractMultiScaleArrayHead{B}
    nodes::N
    values::Vector{B}
    end_idxs::Vector{Int}
end

organ1 = Organ([1.1,2.1,3.1], :Shoot, OrganParams(:grows_up))
organ2 = Organ([4.1,5.1,6.1], :Root, OrganParams("grows down"))
organ3 = Organ([1.2,2.2,3.2], :Shoot, OrganParams(true))
organ4 = Organ([4.2,5.2,6.2], :Root, OrganParams(1//3))
plant1 = construct(Plant, (deepcopy(organ1), deepcopy(organ2)), [], PlantSettings(1))
plant2 = construct(Plant, (deepcopy(organ3), deepcopy(organ4)), [], PlantSettings(1.0))
community = construct(Community, (deepcopy(plant1), deepcopy(plant2), []))
scenario = construct(Scenario, (deepcopy(community),))

@inferred getindex(organ1, 1)
@inferred getindex(plant1, 3)
@inferred getindex(community, 4)
@inferred getindex(scenario, 8)

@test scenario[1] == 1.1
@test scenario[2] == 2.1
@test scenario[3] == 3.1
@test scenario[4] == 4.1
@test scenario[5] == 5.1
@test scenario[6] == 6.1
@test scenario[7] == 10.0
@test scenario[8] == 1.2
@test scenario[9] == 2.2
@test scenario[10] == 3.2
@test scenario[11] == 4.2
@test scenario[12] == 5.2
@test scenario[13] == 6.2
@test scenario[14] == 20.0

@test getindices(scenario, 1) == 1:14
@test getindices(scenario, 1, 1) == 1:7
@test getindices(scenario, 1, 2) == 8:14
@test getindices(scenario, 1, 1, 1) == 1:3
@test getindices(scenario, 1, 1, 2) == 4:6
@test getindices(scenario, 1, 1, 3) == 7:7
@test getindices(scenario, 1, 2, 1) == 8:10

organ_ode = function (dorgan,organ,p,t)
    m = mean(organ)
    for i in eachindex(organ)
        dorgan[i] = -m*organ[i]
    end
end
f = function (dscenario,scenario,p,t)
    for (organ, y, z) in LevelIterIdx(scenario, 2)
        organ_ode(@view(dscenario[y:z]),organ,p,t)
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

println("ODE with tuple nodes")

prob = ODEProblem(f, scenario, (0.0, 1.0))
test_scenario = deepcopy(scenario)

sol = solve(prob, Tsit5(), callback=growing_cb, tstops=tstop)

@test length(sol[end]) == 23
