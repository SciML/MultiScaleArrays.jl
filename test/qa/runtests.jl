using JET
using MultiScaleArrays
using Test

@testset "Quality Assurance" begin
    @testset "JET static analysis" begin
        report = JET.report_package(
            MultiScaleArrays;
            target_modules = (MultiScaleArrays,)
        )
        @test length(JET.get_reports(report)) == 0
    end
end
