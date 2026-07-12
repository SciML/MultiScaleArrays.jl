using SciMLTesting, MultiScaleArrays, JET
using Test

# Aqua `ambiguities` and `deps_compat` are tracked-broken in
# https://github.com/SciML/MultiScaleArrays.jl/issues/142 (5 `construct`/`ldiv!`
# ambiguities; LinearAlgebra and Random lack `[compat]` entries). They are disabled
# in `Aqua.test_all` and recorded as `@test_broken` placeholders via `aqua_broken`;
# remove them from `aqua_broken` once the issue is resolved.
#
# JET runs as a hard check (the previous `alg_needs_extra_process` finding no longer
# fires). The ExplicitImports `ignore` lists cover non-public accesses to other
# packages' names (internals that become public as those base libraries declare
# `public`); each ignored name is grouped by its source package. Every entry below
# was verified to still flag with the lists emptied.
run_qa(
    MultiScaleArrays;
    explicit_imports = true,
    api_docs_kwargs = (; rendered = true),
    aqua_broken = (:ambiguities, :deps_compat),
    ei_kwargs = (;
        all_qualified_accesses_are_public = (;
            ignore = (
                # Base.Broadcast internals
                :AbstractArrayStyle, :Broadcasted, :DefaultArrayStyle,
                :_broadcast_getindex_eltype,
                # ForwardDiff internals
                :DerivativeConfig, :Dual, :Tag,
                # FiniteDiff internals
                :GradientCache, :JacobianCache,
                # OrdinaryDiffEqCore / OrdinaryDiffEqRosenbrock cache types
                :OrdinaryDiffEqCache, :RosenbrockMutableCache,
                # OrdinaryDiffEqDifferentiation internals
                :resize_grad_config!, :resize_jac_config!,
                # Base internals
                :typename,
            ),
        ),
        all_explicit_imports_are_public = (;
            ignore = (
                # non-public RecursiveArrayTools name used in level_iterations.jl
                :chain,
            ),
        ),
    ),
)
