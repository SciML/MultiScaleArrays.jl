using Pkg
using SafeTestsets, Test

const GROUP = get(ENV, "GROUP", "All")

# The QA group (Aqua/JET) lives in an isolated environment under test/qa so its
# tooling compat bounds don't constrain the base test resolve. Activate it before
# running the QA tests. On Julia < 1.11 the [sources] table is ignored, so the
# in-repo root path is developed explicitly.
function activate_qa_env()
    Pkg.activate(joinpath(@__DIR__, "qa"))
    if VERSION < v"1.11.0-DEV.0"
        Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    end
    return Pkg.instantiate()
end

@time begin
    if GROUP == "All" || GROUP == "Core"
        @eval using MultiScaleArrays, OrdinaryDiffEq, DiffEqBase, StochasticDiffEq
        @time @safetestset "Tuple Nodes" begin
            include("tuple_nodes.jl")
        end
        @time @safetestset "Bisect Search Tests" begin
            include("bisect_search_tests.jl")
        end
        @time @safetestset "Indexing and Creation Tests" begin
            include("indexing_and_creation_tests.jl")
        end
        @time @safetestset "Values Indexing" begin
            include("values_indexing.jl")
        end
        @time @safetestset "Get Indices Tests" begin
            include("get_indices.jl")
        end
        @time @safetestset "Additional Fields Test" begin
            include("additional_fields_test.jl")
        end
        @time @safetestset "Dynamic DiffEq Tests" begin
            include("dynamic_diffeq.jl")
        end
        @time @safetestset "Single Layer DiffEq Tests" begin
            include("single_layer_diffeq.jl")
        end
        @time @safetestset "New Nodes Tests" begin
            include("new_nodes.jl")
        end
    end

    if GROUP == "QA"
        activate_qa_env()
        @time @safetestset "Quality Assurance" begin
            include("qa/qa.jl")
        end
    end
end
