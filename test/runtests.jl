using MultiScaleModels, DifferentialEquations
using Base.Test

immutable Cell{T} <: MultiScaleModelLeaf{T}
  x::Vector{T}
end
immutable Population{T<:AbstractMultiScaleModel} <: AbstractMultiScaleModel{Float64}
  x::Vector{T}
  end_idxs::Vector{Int}
end
immutable Tissue{T<:AbstractMultiScaleModel} <: AbstractMultiScaleModel{Float64}
  x::Vector{T}
  end_idxs::Vector{Int}
end
immutable Embryo{T<:AbstractMultiScaleModel} <: MultiScaleModelHead{Float64}
  x::Vector{T}
  end_idxs::Vector{Int}
end
a = collect(1:3:30)
@test MultiScaleModels.bisect_search(a,20) == 7
@test MultiScaleModels.bisect_search(a,13) == 4
for i = 1:28
  @test MultiScaleModels.bisect_search(a,i) == ((i+1)÷3)
end


cell1 = Cell([1.0;2.0;3.0])
cell2 = Cell([4.0;5])

sim_cell = similar(cell1)

@test length(cell1) == 3
sim_cell.x[:] = [1.;2;3]
@test sim_cell.x == cell1.x
@test !(sim_cell.x === cell1.x)

cell3 = Cell([3.0;2.0;5.0])
cell4 = Cell([4.0;6])

@test cell4[2] == 6.0

p = construct(Population,deepcopy([cell1,cell2]))

@test p[1] == 1
@test p[2] == 2
@test p[3] == 3
@test p[4] == 4
@test p[5] == 5

sim_p  = similar(p)
@test length(sim_p) == length(p)
@test length(sim_p.x[1]) == length(p.x[1])
@test !(sim_p.x[1]===p.x[1])

p2 = construct(Population,deepcopy([cell3,cell4]))

tis = construct(Tissue,deepcopy([p,p2]))

sim_tis = similar(tis)
@test length(sim_tis) == length(tis)

@test tis[9] == 4

@test tis.x[2].end_idxs[1] == length(cell1)
@test tis.x[2].end_idxs[2] == length(cell1)+length(cell2)
@test tis.x[2].end_idxs[2] == length(p)

tis2 = construct(Tissue,deepcopy([p2,p]))

em = construct(Embryo,deepcopy([tis,tis2]))

tis3 = construct(Tissue,deepcopy([p,p2]))

@test length(em) == 20

add_daughter!(em,tis3)

@test length(em) == 30
@test em.x[1] != tis #There was a deepcopy
@test em.x[3] == tis3 #No deepcopy

@test em[12] == 2

em[12] = 50.0
@test em.x[2].x[1].x[1].x[2] == 50
p3 = construct(Population,deepcopy([cell1,cell2]))
add_daughter!(em,p3,2)
@test length(em.x[2].x) == 3
@test length(em.x[2]) == 15
@test em.end_idxs[2] == 25

cell3 = Cell([20.0;11;13])
add_daughter!(em,cell3,2,3)

@test length(em.x[2].x) == length(em.x[2].end_idxs)

remove_daughter!(em,3)
@test length(em.end_idxs)==2
add_daughter!(em,tis3)
remove_daughter!(em,1)
@test em.end_idxs[1]==18

remove_daughter!(em,2,1)
@test em.end_idxs[2] == 23
remove_daughter!(em,2,1)
@test length(em.x) == 1 # Should drop the 0
remove_daughter!(em,1,1)
@test length(em) == 13

for i in eachindex(em)
  em[i] *= 2
end

@test em[12] == 22

### Test math

#broadcast_getindex(cell1,1)

g = (x,y) -> x*y
@test g.(cell1,2) == [2.;4;6]
# cell1 .= g.(cell1,2) How to broadcast right???

cell3 = cell1.+2
@test (p.+2)[1] - p[1] == 2
cell1./2

f = function (t,u,du)
  for i in eachindex(u)
    du[i] = 0.52*u[i]
  end
end

prob = ODEProblem(f,em)
sol = solve(prob)
