using MultiScaleArrays, Aqua
@testset "Aqua" begin
    Aqua.find_persistent_tasks_deps(MultiScaleArrays)
    Aqua.test_ambiguities(MultiScaleArrays, recursive = false)
    Aqua.test_deps_compat(MultiScaleArrays)
    Aqua.test_piracies(MultiScaleArrays)
    Aqua.test_project_extras(MultiScaleArrays)
    Aqua.test_stale_deps(MultiScaleArrays)
    Aqua.test_unbound_args(MultiScaleArrays)
    Aqua.test_undefined_exports(MultiScaleArrays)
end
