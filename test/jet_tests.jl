using JET
using MultiScaleArrays

@testset "JET static analysis" begin
    # Run JET analysis on the package
    report = JET.report_package(MultiScaleArrays;
        target_modules = (MultiScaleArrays,))
    @test length(JET.get_reports(report)) == 0
end
