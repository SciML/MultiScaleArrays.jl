using MultiScaleArrays
using OrdinaryDiffEq, DiffEqBase, Base.Test, StochasticDiffEq

#=
immutable Cell{B} <: AbstractMultiScaleArrayLeaf{B}
  x::Vector{B}
end
immutable Population{T<:AbstractMultiScaleArray,B<:Number} <: AbstractMultiScaleArray{B}
  x::Vector{T}
  y::Vector{B}
  end_idxs::Vector{Int}
end
immutable Tissue{T<:AbstractMultiScaleArray,B<:Number} <: AbstractMultiScaleArray{B}
  x::Vector{T}
  y::Vector{B}
  end_idxs::Vector{Int}
end
immutable Embryo{T<:AbstractMultiScaleArray,B<:Number} <: AbstractMultiScaleArrayHead{B}
  x::Vector{T}
  y::Vector{B}
  end_idxs::Vector{Int}
end
=#

cell1 = Cell([1.0;2.0;3.0])
cell2 = Cell([4.0;5])
population = construct(Population,deepcopy([cell1,cell2])) # Make a Population from cells
cell3 = Cell([3.0;2.0;5.0])
cell4 = Cell([4.0;6])
population2 = construct(Population,deepcopy([cell3,cell4]))
tissue1 = construct(Tissue,deepcopy([population,population2])) # Make a Tissue from Populations
tissue2 = construct(Tissue,deepcopy([population2,population]))
embryo = construct(Embryo,deepcopy([tissue1,tissue2])) # Make an embryo from Tissues

cell_ode = function (t,cell,dcell)
  m = mean(cell)
  for i in eachindex(cell)
    dcell[i] = -m*cell[i]
  end
end

f = function (t,embryo,dembryo)
  for (cell,y,z) in LevelIterIdx(embryo,2)
    cell_ode(t,cell,@view dembryo[y:z])
  end
end

tstop = [0.5]

condition = function (t,u,integrator)
  t ∈ tstop
end

affect! = function (integrator)
  add_daughter!(integrator,integrator.u[1,1,1],1,1)
end

growing_cb = DiscreteCallback(condition,affect!)

prob = ODEProblem(f,embryo,(0.0,1.0))
test_embryo = deepcopy(embryo)

sol = solve(prob,Tsit5(),callback=growing_cb,tstops=tstop)

sol = solve(prob,Rosenbrock23(),callback=growing_cb,tstops=tstop)

@test length(sol[end]) == 23

affect_del! = function (integrator)
  remove_daughter!(integrator,1,1,1)
end

shrinking_cb = DiscreteCallback(condition,affect_del!)

sol = solve(prob,Tsit5(),callback=shrinking_cb,tstops=tstop)

sol = solve(prob,Rosenbrock23(),callback=shrinking_cb,tstops=tstop)

@test length(sol[end]) == 17


g = function (t,u,du)
  for i in eachindex(u)
    du[i] = 0.1u[i]
  end
end
prob = SDEProblem(f,g,embryo,(0.0,1.0))

sol = solve(prob,SRIW1(),callback=growing_cb,tstops=tstop)

sol = solve(prob,SRI(),callback=growing_cb,tstops=tstop)

sol = solve(prob,SRA(),callback=growing_cb,tstops=tstop)

sol = solve(prob,SRA1(),callback=growing_cb,tstops=tstop)

sol = solve(prob,RKMil(),callback=growing_cb,dt=1/10,tstops=tstop)

sol = solve(prob,EM(),dt=1/10,callback=growing_cb,tstops=tstop)

@test length(sol[end]) == 23

sol = solve(prob,SRIW1(),callback=shrinking_cb,tstops=tstop)

sol = solve(prob,SRI(),callback=shrinking_cb,tstops=tstop)

sol = solve(prob,SRA(),callback=shrinking_cb,tstops=tstop)

sol = solve(prob,SRA1(),callback=shrinking_cb,tstops=tstop)

sol = solve(prob,RKMil(),dt=1/10,callback=shrinking_cb,tstops=tstop)

sol = solve(prob,EM(),dt=1/10,callback=shrinking_cb,tstops=tstop)

@test length(sol[end]) == 17
