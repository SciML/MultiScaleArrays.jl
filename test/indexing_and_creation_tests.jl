using MultiScaleArrays, DiffEqBase, OrdinaryDiffEq, StochasticDiffEq, Test,
      Random

### Setup a hierarchy

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

#### End Setup

cell1 = Cell([1.0; 2.0; 3.0])
cell2 = Cell([4.0; 5])

cell1 .+ cell1

sim_cell = similar(cell1)

@test length(cell1) == 3
sim_cell.values[:] = [1.0; 2; 3]
@test sim_cell.values == cell1.values
@test !(sim_cell.values === cell1.values)

cell3 = Cell([3.0; 2.0; 5.0])
cell4 = Cell([4.0; 6])

@test cell4[2] == 6.0

p = construct(Population, deepcopy([cell1, cell2]))

@test p[1] == 1
@test p[2] == 2
@test p[3] == 3
@test p[4] == 4
@test p[5] == 5

sim_p = similar(p)
sim_p_arr = similar(p, axes(p))

@test sim_p isa Population
@test sim_p_arr isa Vector{Float64}
@test length(sim_p) == length(p)
@test length(sim_p.nodes[1]) == length(p.nodes[1])
@test !(sim_p.nodes[1] === p.nodes[1])

p2 = construct(Population, deepcopy([cell3, cell4]))

tis = construct(Tissue, deepcopy([p, p2]))

sim_tis = similar(tis)
@test length(sim_tis) == length(tis)

@test tis[9] == 4

@test tis.nodes[2].end_idxs[1] == length(cell1)
@test tis.nodes[2].end_idxs[2] == length(cell1) + length(cell2)
@test tis.nodes[2].end_idxs[2] == length(p)

tis2 = construct(Tissue, deepcopy([p2, p]))

em = construct(Embryo, deepcopy([tis, tis2]))

tis3 = construct(Tissue, deepcopy([p, p2]))

@test length(em) == 20

add_node!(em, tis3)

em_save = deepcopy(em)

print_human_readable(em)
print_human_readable(em; n_char_per_name = 2)
print_human_readable(em; n_char_per_name = 2, fields = [:values])
print_human_readable(em.nodes[1].nodes[1]; fields = [:values])

@test length(em) == 30
@test em.nodes[1] != tis #There was a deepcopy
@test em.nodes[3] == tis3 #No deepcopy

@test em[12] == 2

em[12] = 50.0
@test em.nodes[2].nodes[1].nodes[1].values[2] == 50
p3 = construct(Population, deepcopy([cell1, cell2]))
add_node!(em, p3, 2)
@test length(em.nodes[2].nodes) == 3
@test length(em.nodes[2]) == 15
@test em.end_idxs[2] == 25

cell3 = Cell([20.0; 11; 13])
add_node!(em, cell3, 2, 3)

@test length(em.nodes[2].nodes) == length(em.nodes[2].end_idxs)

remove_node!(em, 3)
@test length(em.end_idxs) == 2
add_node!(em, tis3)
remove_node!(em, 1)
@test em.end_idxs[1] == 18

remove_node!(em, 2, 1)
@test em.end_idxs[2] == 23
remove_node!(em, 2, 1)
@test length(em.nodes) == 1 # Should drop the 0
remove_node!(em, 1, 1)

@test length(em) == 13

for i in eachindex(em)
    em[i] *= 2
end

@test em[12] == 22

### Test math

g = (x, y) -> x * y
@test g.(cell1, 2) == [2.0; 4; 6]
# cell1 .= g.(cell1, 2) How to broadcast right???

cell3 = cell1 .+ 2

@test cell3 isa AbstractMultiScaleArray

cell3 = similar(cell1)
cell3 .= [1, 2, 3]
cell3 .= 2cell3

@test (p .+ 2)[1] - p[1] == 2
cell1 ./ 2

size(cell1)

t = 2
p / zero(t)
p ./ randn(length(p))

f = function (du, u, p, t)
    for i in eachindex(u)
        du[i] = 0.42 * u[i]
    end
end
g = function (du, u, p, t)
    for i in eachindex(u)
        du[i] = 0.42 * u[i]
    end
end

vem = @view [em, em][1:2]

prob = ODEProblem(f, em, (0.0, 1500.0))
@time sol1 = solve(prob, Tsit5(), save_everystep = false)

prob = ODEProblem(f, em[:], (0.0, 1500.0))
@time sol2 = solve(prob, Tsit5(), save_everystep = false)
@test sol1.t == sol2.t

prob = ODEProblem(f, em, (0.0, 1500.0))
sol1 = solve(prob, Tsit5())

# Check stepping behavior matches array
Random.seed!(100)
prob = SDEProblem(f, g, em, (0.0, 1000.0))
@time sol1 = solve(prob, SRIW1(), progress = false, abstol = 1e-2, reltol = 1e-2,
    save_everystep = false)

Random.seed!(100)
prob = SDEProblem(f, g, em[:], (0.0, 1000.0))
@time sol2 = solve(prob, SRIW1(), progress = false, abstol = 1e-2, reltol = 1e-2,
    save_everystep = false)
sol1.t == sol2.t

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

#sol = solve(prob, EM())

em2 = similar(em)
recursivecopy!(em2, em)
@test em[5] == em2[5]
@test em != em2

## Level iterators

em = em_save
level_iter(em, 1) == em.nodes
for (i, p) in enumerate(level_iter(em, 2))
    if i == 1
        @test p == em.nodes[1].nodes[1]
    elseif i == 2
        @test p == em.nodes[1].nodes[2]
    elseif i == 3
        @test p == em.nodes[2].nodes[1]
    elseif i == 4
        @test p == em.nodes[2].nodes[2]
    elseif i == 5
        @test p == em.nodes[3].nodes[1]
    elseif i == 6
        @test p == em.nodes[3].nodes[2]
    elseif i > 6
        error("you shouldn't be here")
    end
end

em_arr = em[:]

for (x, y, z) in LevelIterIdx(em, 1)
    @test maximum(em_arr[y:z] - x[:]) == 0
end

for (x, y, z) in LevelIterIdx(em, 2)
    @test maximum(em_arr[y:z] - x[:]) == 0
end

for (x, y, z) in LevelIterIdx(em, 3)
    @test maximum(em_arr[y:z] - x[:]) == 0
end

em2 = copy(em)

for (pop1, pop2) in LevelIter(2, em, em2)
    @test pop1 isa Population
    @test pop1[:] == pop2[:]
end

for (cell1, cell2) in LevelIter(3, em, em2)
    @test cell1 isa Cell
    @test cell1[:] == cell2[:]
end

### Non-Empty y

p = construct(Population, deepcopy([cell1, cell2]), [1.0; 2; 3])

@test p[1] == 1
@test p[2] == 2
@test p[3] == 3
@test p[4] == 4
@test p[5] == 5
@test p[6] == 1
@test p[7] == 2
@test p[8] == 3

sim_p = similar(p)
@test length(sim_p) == length(p)
@test length(sim_p.nodes[1]) == length(p.nodes[1])
@test !(sim_p.nodes[1] === p.nodes[1])

p2 = construct(Population, deepcopy([cell3, cell4]), [11.0; 12; 13])

tis = construct(Tissue, deepcopy([p, p2]))

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
