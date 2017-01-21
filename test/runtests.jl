using MultiScaleModels, OrdinaryDiffEq, DiffEqBase, StochasticDiffEq
using Base.Test

#=
macro define_hierarchy(BottomType,names)
 quote
   name = $(esc(names))[1]
   immutable $(esc(name)){$(esc(BottomType))} <: MultiScaleModelLeaf{$(esc(BottomType))}
     x::Vector{$(esc(BottomType))}
   end

   $((quote
     immutable $(esc(names.args[i])){$(esc(BottomType))<:AbstractMultiScaleModel} <: AbstractMultiScaleModel{Float64}
       x::Vector{T}
       end_idxs::Vector{Int}
     end
      end for i in 2:length(names.args))...)

 end
end
      =#
#=
@show macroexpand(:(@define_hierarchy(Float64,[:Cell,:Population,:Tissue,:Embryo])))

@define_hierarchy(Float64,[:Cell,:Population,:Tissue,:Embryo])
=#

### Setup a hierarchy

immutable Cell{B} <: MultiScaleModelLeaf{B}
  x::Vector{B}
end
immutable Population{T<:AbstractMultiScaleModel,B<:Number} <: AbstractMultiScaleModel{B}
  x::Vector{T}
  y::Vector{B}
  end_idxs::Vector{Int}
end
immutable Tissue{T<:AbstractMultiScaleModel,B<:Number} <: AbstractMultiScaleModel{B}
  x::Vector{T}
  y::Vector{B}
  end_idxs::Vector{Int}
end
immutable Embryo{T<:AbstractMultiScaleModel,B<:Number} <: MultiScaleModelHead{B}
  x::Vector{T}
  y::Vector{B}
  end_idxs::Vector{Int}
end

Cell(x::Tuple) = Cell(Vector{Float64}(x))
function Population(x::Tuple)
  p = Population(Vector{Cell}(x[1]))
  for i in 1:x[1]
    p[i] = Cell(x[2:end])
  end
end
function Tissue(x::Tuple)
  tis = Tissue(Vector{Population}(x[1]))
  for i in 1:x[1]
    tis[i] = Population(x[2:end])
  end
end
function Embryo(x::Tuple)
  tis = Embryo(Vector{Tissue}(x[1]))
  for i in 1:x[1]
    tis[i] = Tissue(x[2:end])
  end
end

#### End Setup

a = collect(3:3:30)
@test MultiScaleModels.bisect_search(a,20) == 7
@test MultiScaleModels.bisect_search(a,13) == 5
@test MultiScaleModels.bisect_search(a,12) == 4
for i = 1:30
  @test MultiScaleModels.bisect_search(a,i) == ((i-1)÷3)+1 #+1 for 1-based indexing
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
sim_p_arr  = similar(p,indices(p))

@test typeof(sim_p) <: Population
@test typeof(sim_p_arr) <: Vector{Float64}
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

em_save = deepcopy(em)

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

size(cell1)
Cell(size(cell1))

t = 2
p/zero(t)

size(p)

f = function (t,u,du)
  for i in eachindex(u)
    du[i] = 0.42*u[i]
  end
end
g = function (t,u,du)
  for i in eachindex(u)
    du[i] = 0.42*u[i]
  end
end

vem = @view [em,em][1:2]

prob = ODEProblem(f,em,(0.0,1500.0))
@time sol1 = solve(prob,Tsit5(),save_timeseries=false)

prob = ODEProblem(f,em[:],(0.0,1500.0))
@time sol2 = solve(prob,Tsit5(),save_timeseries=false)
@test sol1.t == sol2.t

prob = ODEProblem(f,em,(0.0,1500.0))
sol1 = solve(prob,Tsit5())

# Check stepping behavior matches array
srand(100)
prob = SDEProblem(f,g,em,(0.0,1000.0))
@time sol1 = solve(prob,SRIW1(),progress=false,abstol=1e-2,reltol=1e-2,save_timeseries=false)

srand(100)
prob = SDEProblem(f,g,em[:],(0.0,1000.0))
@time sol2 = solve(prob,SRIW1(),progress=false,abstol=1e-2,reltol=1e-2,save_timeseries=false)
@test sol1.t == sol2.t

function test_loop(a)
  for i in eachindex(a)
    a[i] = a[i] + 1
  end
end
@time test_loop(em)
@time test_loop(em[:])
@time em[5]
a = em[:]
@time a[1]

#sol = solve(prob,EM())

em2 = similar(em)
recursivecopy!(em2,em)
@test em[5] == em2[5]
@test em != em2

## Level iterators

em = em_save
level_iter(em,1) == em.x
for (i,p) in enumerate(level_iter(em,2))
  if i == 1
    @test p == em.x[1].x[1]
  elseif i == 2
    @test p == em.x[1].x[2]
  elseif i == 3
    @test p == em.x[2].x[1]
  elseif i == 4
    @test p == em.x[2].x[2]
  elseif i == 5
    @test p == em.x[3].x[1]
  elseif i == 6
    @test p == em.x[3].x[2]
  elseif i > 6
    error("you shouldn't be here")
  end
end

em_arr = em[:]

for (x,y,z) in LevelIterIdx(em,1)
  @test maximum(em_arr[y:z]- x[:]) ==0
end

for (x,y,z) in LevelIterIdx(em,2)
  @test maximum(em_arr[y:z]- x[:]) ==0
end

for (x,y,z) in LevelIterIdx(em,3)
  @test maximum(em_arr[y:z]- x[:]) ==0
end

### Non-Empty y

p = construct(Population,deepcopy([cell1,cell2]),[1.;2;3])

@test p[1] == 1
@test p[2] == 2
@test p[3] == 3
@test p[4] == 4
@test p[5] == 5
@test p[6] == 1
@test p[7] == 2
@test p[8] == 3

sim_p  = similar(p)
@test length(sim_p) == length(p)
@test length(sim_p.x[1]) == length(p.x[1])
@test !(sim_p.x[1]===p.x[1])

p2 = construct(Population,deepcopy([cell3,cell4]),[11.;12;13])

tis = construct(Tissue,deepcopy([p,p2]))

@test length(tis) == length(p) + length(p2)

@test tis[1] == 1
@test tis[2] == 2
@test tis[3] == 3
@test tis[4] == 4
@test tis[5] == 5
@test tis[6] == 1
@test tis[7] == 2
@test tis[8] == 3
@test tis[9] == p2[1]
@test tis[10] == p2[2]
@test tis[11] == p2[3]
@test tis[11] == p2[3]
@test tis[12] == p2[4]
@test tis[13] == p2[5]
@test tis[14] == p2[6]
@test tis[15] == p2[7]
@test tis[16] == p2[8]
