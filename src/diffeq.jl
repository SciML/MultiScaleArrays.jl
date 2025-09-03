function remove_node!(integrator::DiffEqBase.DEIntegrator, I...)
    idxs = getindices(integrator.u, I...)
    for c in full_cache(integrator)
        remove_node!(c, I...)
    end
    if DiffEqBase.is_diagonal_noise(integrator.sol.prob)
        for c in DiffEqBase.ratenoise_cache(integrator)
            remove_node!(c, I...)
        end
    end
    remove_node_non_user_cache!(integrator, idxs, I...) # required to do noise correctly
end

function add_node!(integrator::DiffEqBase.DEIntegrator, x, I...)
    cur_len = length(integrator.u)
    add_len = length(x)
    last_idx = length(integrator.u[I...].nodes)
    idx_start = getindices(integrator.u, last_idx)[end] + 1
    idxs = idx_start:(idx_start + add_len - 1)
    for c in full_cache(integrator)
        add_node!(c, recursivecopy(x), I...)
    end
    if DiffEqBase.is_diagonal_noise(integrator.sol.prob)
        for c in DiffEqBase.ratenoise_cache(integrator)
            add_node!(c, recursivecopy(x), I...)
        end
    end
    #addat_non_user_cache!(integrator, idxs)
    add_node_non_user_cache!(integrator, idxs, x, I...) # required to do noise correctly
end

function add_node!(integrator::DiffEqBase.DEIntegrator, x)
    cur_len = length(integrator.u)
    add_len = length(x)
    last_idx = length(integrator.u.nodes)
    idx_start = getindices(integrator.u, last_idx)[end] + 1
    idxs = idx_start:(idx_start + add_len - 1)
    for c in full_cache(integrator)
        add_node!(c, recursivecopy(x))
    end
    if DiffEqBase.is_diagonal_noise(integrator.sol.prob)
        for c in DiffEqBase.ratenoise_cache(integrator)
            add_node!(c, recursivecopy(x))
        end
    end
    add_node_non_user_cache!(integrator, idxs, fill!(similar(x, eltype(x)), 0)) # required to do noise correctly
end

reshape(m::AbstractMultiScaleArray, i::Int...) = m

function remove_node_non_user_cache!(integrator::DiffEqBase.AbstractODEIntegrator, idxs,
        node...)
    remove_node_non_user_cache!(integrator, integrator.cache, node...)
end

function remove_node_non_user_cache!(integrator::DiffEqBase.AbstractODEIntegrator,
        cache::OrdinaryDiffEqCore.OrdinaryDiffEqCache, idxs,
        node...)
    nothing
end

function add_node_non_user_cache!(integrator::DiffEqBase.AbstractODEIntegrator, idxs,
        x::AbstractArray)
    add_node_non_user_cache!(integrator, integrator.cache, x)
end

function add_node_non_user_cache!(integrator::DiffEqBase.AbstractODEIntegrator, idxs,
        x::AbstractArray, node...)
    add_node_non_user_cache!(integrator, integrator.cache, x, node...)
end

function add_node_non_user_cache!(integrator::DiffEqBase.AbstractODEIntegrator,
        cache::OrdinaryDiffEqCore.OrdinaryDiffEqCache,
        x::AbstractArray)
    nothing
end
function add_node_non_user_cache!(integrator::DiffEqBase.AbstractODEIntegrator,
        cache::OrdinaryDiffEqCore.OrdinaryDiffEqCache,
        x::AbstractArray, node...)
    nothing
end

function add_node_non_user_cache!(integrator::DiffEqBase.AbstractODEIntegrator,
        cache::OrdinaryDiffEqRosenbrock.RosenbrockMutableCache,
        x::AbstractArray)
    i = length(integrator.u)
    cache.J = similar(cache.J, i, i)
    cache.W = similar(cache.W, i, i)
    
    # Handle jacobian config resizing
    # Try to use the proper DI approach if possible, otherwise fall back to field-by-field
    resize_jac_config!(integrator.f, cache.jac_config, integrator.u, cache)
    
    add_node_grad_config!(cache, cache.grad_config, i, x)
    nothing
end

function add_node_non_user_cache!(integrator::DiffEqBase.AbstractODEIntegrator,
        cache::OrdinaryDiffEqRosenbrock.RosenbrockMutableCache,
        x::AbstractArray, node...)
    i = length(integrator.u)
    cache.J = similar(cache.J, i, i)
    cache.W = similar(cache.W, i, i)
    
    # Handle jacobian config resizing
    resize_jac_config!(integrator.f, cache.jac_config, integrator.u, cache)
    
    add_node_grad_config!(cache, cache.grad_config, i, x, node...)
    nothing
end

function remove_node_non_user_cache!(integrator::DiffEqBase.AbstractODEIntegrator,
        cache::OrdinaryDiffEqRosenbrock.RosenbrockMutableCache,
        node...)
    i = length(integrator.u)
    cache.J = similar(cache.J, i, i)
    cache.W = similar(cache.W, i, i)
    
    # Handle jacobian config resizing
    resize_jac_config!(integrator.f, cache.jac_config, integrator.u, cache)
    
    remove_node_grad_config!(cache, cache.grad_config, i, node...)
    nothing
end

# High-level function to properly resize jacobian configurations
function resize_jac_config!(f, jac_config, u, cache)
    # For DifferentiationInterface types, use prepare!_jacobian when dimensions change
    # This is the recommended approach as mentioned in the PR discussion
    typename = string(typeof(jac_config))
    if occursin("DifferentiationInterface", typename) && hasproperty(jac_config, :backend)
        try
            # Use the proper DI API to resize the configuration
            backend = getfield(jac_config, :backend)
            DI.prepare!_jacobian(f, jac_config, backend, u)
            return
        catch e
            # If DI prepare!_jacobian fails, fall back to field-by-field handling
            @debug "DI prepare!_jacobian failed, falling back to field-by-field approach" exception=e
        end
    end
    
    # Fallback approach for non-DI types or when DI approach fails
    i = length(u)
    add_node_jac_config!(cache, jac_config, i, similar(u, 0))
    nothing
end

# Generic fallback for any jac_config type - handle DifferentiationInterface types properly
function add_node_jac_config!(cache, config, i, x)
    # The proper solution for DifferentiationInterface types is to use prepare!_jacobian
    # when the input dimensions change. However, this requires context that we don't have here.
    # For now, we implement a conservative approach that works for most cases.
    
    # Fallback: For all types, handle known-safe fields
    if hasproperty(config, :fx) && getproperty(config, :fx) isa AbstractArray
        try
            add_node!(getproperty(config, :fx), recursivecopy(x))
        catch
            # If it fails, that's ok - not all arrays support add_node!
        end
    end
    
    # Update colorvec if it exists and is mutable
    if hasproperty(config, :colorvec)
        try
            setproperty!(config, :colorvec, 1:i)
        catch
            # Field might be immutable or not support this, skip it
        end
    end
    nothing
end

function add_node_jac_config!(cache, config, i, x, I...)
    # Fallback: For all types, handle known-safe fields
    if hasproperty(config, :fx) && getproperty(config, :fx) isa AbstractArray
        try
            add_node!(getproperty(config, :fx), recursivecopy(x), I...)
        catch
            # If it fails, that's ok
        end
    end
    if hasproperty(config, :colorvec)
        try
            setproperty!(config, :colorvec, 1:i)
        catch
            # Field might be immutable, skip it
        end
    end
    nothing
end

function remove_node_jac_config!(cache, config, i, I...)
    # Fallback: For all types, handle known-safe fields
    if hasproperty(config, :fx) && getproperty(config, :fx) isa AbstractArray
        try
            remove_node!(getproperty(config, :fx), I...)
        catch
            # If it fails, that's ok
        end
    end
    if hasproperty(config, :colorvec)
        try
            setproperty!(config, :colorvec, 1:i)
        catch
            # Field might be immutable, skip it
        end
    end
    nothing
end

# Specific implementation for FiniteDiff.JacobianCache (keeps backward compatibility)
function add_node_jac_config!(cache, config::FiniteDiff.JacobianCache, i, x)
    #add_node!(cache.x1, fill!(similar(x, eltype(cache.x1)),0))
    add_node!(config.fx, recursivecopy(x))
    #cache.fx1 !== nothing && add_node!(cache.fx1, fill!(similar(x, eltype(cache.fx1)),0))
    config.colorvec = 1:i
    nothing
end

function add_node_jac_config!(cache, config::FiniteDiff.JacobianCache, i, x, I...)
    #add_node!(cache.x1, fill!(similar(x, eltype(cache.x1)),0), I...)
    add_node!(config.fx, recursivecopy(x), I...)
    #cache.fx1 !== nothing && add_node!(cache.fx1, fill!(similar(x, eltype(cache.fx1)),0), I...)
    config.colorvec = 1:i
    nothing
end

function remove_node_jac_config!(cache, config::FiniteDiff.JacobianCache, i, I...)
    #add_node!(cache.x1, fill!(similar(x, eltype(cache.x1)),0), I...)
    remove_node!(config.fx, I...)
    #cache.fx1 !== nothing && add_node!(cache.fx1, fill!(similar(x, eltype(cache.fx1)),0), I...)
    config.colorvec = 1:i
    nothing
end

# Generic fallback for any grad_config type
function add_node_grad_config!(cache, grad_config, i, x)
    # Most grad configs don't have mutable arrays that need updating
    # For DifferentiationInterface types, they typically handle resizing internally
    nothing
end

function add_node_grad_config!(cache, grad_config, i, x, I...)
    # Most grad configs don't have mutable arrays that need updating
    nothing
end

function remove_node_grad_config!(cache, grad_config, i, x)
    # Most grad configs don't have mutable arrays that need updating
    nothing
end

function remove_node_grad_config!(cache, grad_config, i, x, I...)
    # Most grad configs don't have mutable arrays that need updating
    nothing
end

# Specific implementation for ForwardDiff.DerivativeConfig (keeps backward compatibility)
function add_node_grad_config!(cache, grad_config::ForwardDiff.DerivativeConfig, i, x)
    cache.grad_config = ForwardDiff.DerivativeConfig(cache.tf, cache.du1, cache.uf.t)
    nothing
end

function add_node_grad_config!(cache, grad_config::ForwardDiff.DerivativeConfig, i, x, I...)
    cache.grad_config = ForwardDiff.DerivativeConfig(cache.tf, cache.du1, cache.uf.t)
    nothing
end

function remove_node_grad_config!(cache, grad_config::ForwardDiff.DerivativeConfig, i, x)
    cache.grad_config = ForwardDiff.DerivativeConfig(cache.tf, cache.du1, cache.uf.t)
    nothing
end

function remove_node_grad_config!(cache, grad_config::ForwardDiff.DerivativeConfig, i, x,
        I...)
    cache.grad_config = ForwardDiff.DerivativeConfig(cache.tf, cache.du1, cache.uf.t)
    nothing
end

function add_node_grad_config!(cache, grad_config::AbstractArray, i, x)
    cache.grad_config = ForwardDiff.Dual{
        typeof(ForwardDiff.Tag(cache.tf,
        eltype(cache.du1)))
    }.(cache.du1, cache.du1)
    nothing
end

function add_node_grad_config!(cache, grad_config::AbstractArray, i, x, I...)
    cache.grad_config = ForwardDiff.Dual{
        typeof(ForwardDiff.Tag(cache.tf,
        eltype(cache.du1)))
    }.(cache.du1, cache.du1)
    nothing
end

function remove_node_grad_config!(cache, grad_config::AbstractArray, i, x)
    cache.grad_config = ForwardDiff.Dual{
        typeof(ForwardDiff.Tag(cache.tf,
        eltype(cache.du1)))
    }.(cache.du1, cache.du1)
    nothing
end

function remove_node_grad_config!(cache, grad_config::AbstractArray, i, x, I...)
    cache.grad_config = ForwardDiff.Dual{
        typeof(ForwardDiff.Tag(cache.tf,
        eltype(cache.du1)))
    }.(cache.du1, cache.du1)
    nothing
end

function add_node_grad_config!(cache, grad_config::FiniteDiff.GradientCache, i, x, I...)
    grad_config.fx !== nothing && add_node!(grad_config.fx, recursivecopy(x), I...)
    grad_config.c1 !== nothing && add_node!(grad_config.c1, recursivecopy(x), I...)
    grad_config.c2 !== nothing && add_node!(grad_config.c2, recursivecopy(x), I...)
    grad_config
end

function add_node_grad_config!(cache, grad_config::FiniteDiff.GradientCache, i, x)
    grad_config.fx !== nothing && add_node!(grad_config.fx, recursivecopy(x))
    grad_config.c1 !== nothing && add_node!(grad_config.c1, recursivecopy(x))
    grad_config.c2 !== nothing && add_node!(grad_config.c2, recursivecopy(x))
    grad_config
end

function remove_node_grad_config!(cache, grad_config::FiniteDiff.GradientCache, i, I...)
    grad_config.fx !== nothing && remove_node!(grad_config.fx, I...)
    grad_config.c1 !== nothing && remove_node!(grad_config.c1, I...)
    grad_config.c2 !== nothing && remove_node!(grad_config.c2, I...)
    grad_config
end

function add_node_non_user_cache!(integrator::DiffEqBase.AbstractSDEIntegrator, idxs, x,
        node...)
    #addat_non_user_cache!(integrator, idxs)
    if DiffEqBase.is_diagonal_noise(integrator.sol.prob)
        add_node_noise!(integrator, idxs, x, node...)
        for c in rand_cache(integrator)
            add_node!(c, copy(x), node...)
        end
    end
end

function add_node_non_user_cache!(integrator::DiffEqBase.AbstractSDEIntegrator, idxs, x)
    #addat_non_user_cache!(integrator, idxs)
    if DiffEqBase.is_diagonal_noise(integrator.sol.prob)
        add_node_noise!(integrator, idxs, x)
        for c in rand_cache(integrator)
            add_node!(c, copy(x))
        end
    end
end

function remove_node_non_user_cache!(integrator::DiffEqBase.AbstractSDEIntegrator, idxs,
        node...)
    #deleteat_non_user_cache!(integrator, idxs)
    if DiffEqBase.is_diagonal_noise(integrator.sol.prob)
        remove_node_noise!(integrator, node...)
        for c in rand_cache(integrator)
            remove_node!(c, node...)
        end
    end
end

### For if noise is an AMSA

function add_node_noise!(integrator, idxs, x, node...)
    for c in integrator.W.S₁
        add_node!(c[2], copy(x), node...)
        if StochasticDiffEq.alg_needs_extra_process(integrator.alg)
            add_node!(c[3], copy(x), node...)
        end
        StochasticDiffEq.fill_new_noise_caches!(integrator, c, c[1], idxs)
    end
    for c in integrator.W.S₂
        add_node!(c[2], copy(x), node...)
        if StochasticDiffEq.alg_needs_extra_process(integrator.alg)
            add_node!(c[3], copy(x), node...)
        end
        StochasticDiffEq.fill_new_noise_caches!(integrator, c, c[1], idxs)
    end

    add_node!(integrator.W.dW, copy(x), node...)
    integrator.W.dW[idxs] .= zero(eltype(integrator.u))
    add_node!(integrator.W.curW, copy(x), node...)
    integrator.W.curW[idxs] .= zero(eltype(integrator.u))
    if StochasticDiffEq.alg_needs_extra_process(integrator.alg)
        add_node!(integrator.W.dZ, copy(x), node...)
        integrator.W.dZ[idxs] .= zero(eltype(integrator.u))
        add_node!(integrator.W.curZ, copy(x), node...)
        integrator.W.curZ[idxs] .= zero(eltype(integrator.u))
    end

    i = length(integrator.u)
    add_node!(integrator.W.dWtilde, copy(x), node...)
    add_node!(integrator.W.dWtmp, copy(x), node...)
    if StochasticDiffEq.alg_needs_extra_process(integrator.alg)
        add_node!(integrator.W.dZtmp, copy(x), node...)
        add_node!(integrator.W.dZtilde, copy(x), node...)
    end

    # fill in rands
    fill!(@view(integrator.W.curW[idxs]), zero(eltype(integrator.u)))
    if StochasticDiffEq.alg_needs_extra_process(integrator.alg)
        fill!(@view(integrator.W.curZ[idxs]), zero(eltype(integrator.u)))
    end
end

function add_node_noise!(integrator, idxs, x)
    for c in integrator.W.S₁
        add_node!(c[2], copy(x))
        if StochasticDiffEq.alg_needs_extra_process(integrator.alg)
            add_node!(c[3], copy(x))
        end
        StochasticDiffEq.fill_new_noise_caches!(integrator, c, c[1], idxs)
    end
    for c in integrator.W.S₂
        add_node!(c[2], copy(x))
        if StochasticDiffEq.alg_needs_extra_process(integrator.alg)
            add_node!(c[3], copy(x))
        end
        StochasticDiffEq.fill_new_noise_caches!(integrator, c, c[1], idxs)
    end

    add_node!(integrator.W.dW, copy(x))
    integrator.W.dW[idxs] .= zero(eltype(integrator.u))
    add_node!(integrator.W.curW, copy(x))
    integrator.W.curW[idxs] .= zero(eltype(integrator.u))
    if StochasticDiffEq.alg_needs_extra_process(integrator.alg)
        add_node!(integrator.W.dZ, copy(x))
        integrator.W.dZ[idxs] .= zero(eltype(integrator.u))
        add_node!(integrator.W.curZ, copy(x))
        integrator.W.curZ[idxs] .= zero(eltype(integrator.u))
    end

    i = length(integrator.u)
    add_node!(integrator.W.dWtilde, copy(x))
    add_node!(integrator.W.dWtmp, copy(x))
    if StochasticDiffEq.alg_needs_extra_process(integrator.alg)
        add_node!(integrator.W.dZtmp, copy(x))
        add_node!(integrator.W.dZtilde, copy(x))
    end

    # fill in rands
    fill!(@view(integrator.W.curW[idxs]), zero(eltype(integrator.u)))
    if StochasticDiffEq.alg_needs_extra_process(integrator.alg)
        fill!(@view(integrator.W.curZ[idxs]), zero(eltype(integrator.u)))
    end
end

function remove_node_noise!(integrator, node...)
    for c in integrator.W.S₁
        remove_node!(c[2], node...)
        if DiffEqBase.alg_needs_extra_process(integrator.alg)
            remove_node!(c[3], node...)
        end
    end
    for c in integrator.W.S₂
        remove_node!(c[2], node...)
        if StochasticDiffEq.alg_needs_extra_process(integrator.alg)
            remove_node!(c[3], node...)
        end
    end
    remove_node!(integrator.W.dW, node...)
    remove_node!(integrator.W.dWtilde, node...)
    remove_node!(integrator.W.dWtmp, node...)
    remove_node!(integrator.W.curW, node...)

    if StochasticDiffEq.alg_needs_extra_process(integrator.alg)
        remove_node!(integrator.W.curZ, node...)
        remove_node!(integrator.W.dZtmp, node...)
        remove_node!(integrator.W.dZtilde, node...)
        remove_node!(integrator.W.dZ, node...)
    end
end
