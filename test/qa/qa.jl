using MultiScaleArrays, Aqua, JET
using Test

@testset "Aqua" begin
    Aqua.test_all(MultiScaleArrays)
end

@testset "JET" begin
    JET.test_package(MultiScaleArrays; target_defined_modules = true)
end
