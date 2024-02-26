using MultiScaleArrays, OrdinaryDiffEq, DiffEqBase, StochasticDiffEq, SafeTestsets
using Test

@time @testset "Bisect Search Tests" include("bisect_search_tests.jl")
@time @testset "Indexing and Creation Tests" include("indexing_and_creation_tests.jl")
@time @testset "Values Indexing" include("values_indexing.jl")
@time @testset "Get Indices Tests" include("get_indices.jl")
@time @testset "Additional Fields Test" include("additional_fields_test.jl")
@time @testset "Dynamic DiffEq Tests" include("dynamic_diffeq.jl")
@time @testset "Single Layer DiffEq Tests" include("single_layer_diffeq.jl")
@time @safetestset "New Nodes Tests" include("new_nodes.jl")
