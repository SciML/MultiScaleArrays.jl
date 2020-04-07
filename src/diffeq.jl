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
    remove_node_non_user_cache!(integrator, I...) # required to do noise correctly
end

function add_node!(integrator::DiffEqBase.DEIntegrator, x, I...)
    cur_len = length(integrator.u)
    add_len = length(x)
    last_idx = length(integrator.u[I...].nodes)
    idx_start = getindices(integrator.u, last_idx)[end] + 1
    idxs = idx_start:idx_start+add_len-1
    for c in full_cache(integrator)
        add_node!(c, fill!(similar(x, eltype(c)),0), I...)
    end
    if DiffEqBase.is_diagonal_noise(integrator.sol.prob)
      for c in DiffEqBase.ratenoise_cache(integrator)
          add_node!(c, fill!(similar(x, eltype(c)),0), I...)
      end
    end
    add_node_non_user_cache!(integrator, x, I...) # required to do noise correctly
end

function add_node!(integrator::DiffEqBase.DEIntegrator, x)
    cur_len = length(integrator.u)
    add_len = length(x)
    last_idx = length(integrator.u.nodes)
    idx_start = getindices(integrator.u, last_idx)[end] + 1
    idxs = idx_start:idx_start+add_len-1
    for c in full_cache(integrator)
        add_node!(c, fill!(similar(x, eltype(c)),0))
    end
    if DiffEqBase.is_diagonal_noise(integrator.sol.prob)
      for c in DiffEqBase.ratenoise_cache(integrator)
          add_node!(c, fill!(similar(x, eltype(c)),0))
      end
    end
    add_node_non_user_cache!(integrator, idxs, fill!(similar(x, eltype(x)),0)) # required to do noise correctly
end

reshape(m::AbstractMultiScaleArray, i::Int...) = m

function remove_node_non_user_cache!(integrator::DiffEqBase.AbstractODEIntegrator,node...)
  remove_node_non_user_cache!(integrator,integrator.cache,node...)
end

remove_node_non_user_cache!(integrator::DiffEqBase.AbstractODEIntegrator,
                         cache::OrdinaryDiffEq.OrdinaryDiffEqCache,node...) = nothing

function remove_node_non_user_cache!(integrator::DiffEqBase.AbstractSDEIntegrator,node...)
  if DiffEqBase.is_diagonal_noise(integrator.sol.prob)
    remove_node_noise!(integrator,node...)
    for c in rand_cache(integrator)
      remove_node!(c,node...)
    end
  end
end

function remove_node_noise!(integrator,node)
  for c in integrator.W.S₁
    remove_node!(c[2],node...)
    if DiffEqBase.alg_needs_extra_process(integrator.alg)
      remove_node!(c[3],node...)
    end
  end
  for c in integrator.W.S₂
    remove_node!(c[2],node...)
    if StochasticDiffEq.alg_needs_extra_process(integrator.alg)
      remove_node!(c[3],node...)
    end
  end
  remove_node!(integrator.W.dW,node...)
  remove_node!(integrator.W.dWtilde,node...)
  remove_node!(integrator.W.dWtmp,node...)
  remove_node!(integrator.W.curW,node...)

  if StochasticDiffEq.alg_needs_extra_process(integrator.alg)
    remove_node!(integrator.W.curZ,node...)
    remove_node!(integrator.W.dZtmp,node...)
    remove_node!(integrator.W.dZtilde,node...)
    remove_node!(integrator.W.dZ,node...)
  end
end

function add_node_non_user_cache!(integrator::DiffEqBase.AbstractODEIntegrator,x::AbstractArray)
  add_node_non_user_cache!(integrator,integrator.cache,x)
end

function add_node_non_user_cache!(integrator::DiffEqBase.AbstractODEIntegrator,x::AbstractArray,node...)
  add_node_non_user_cache!(integrator,integrator.cache,x,node...)
end

add_node_non_user_cache!(integrator::DiffEqBase.AbstractODEIntegrator,
                      cache,idxs,x) = nothing
add_node_non_user_cache!(integrator::DiffEqBase.AbstractODEIntegrator,
                      cache,idxs,x,node...) = nothing

function add_node_non_user_cache!(integrator::DiffEqBase.AbstractODEIntegrator,
                      cache::OrdinaryDiffEq.RosenbrockMutableCache,idxs,x)
  i = length(integrator.u)
  cache.J = similar(cache.J,i,i)
  cache.W = similar(cache.W,i,i)
  add_node_jac_config!(cache.jac_config,i,x)
  add_node_grad_config!(cache.grad_config,i,x)
  nothing
end

function add_node_non_user_cache!(integrator::DiffEqBase.AbstractODEIntegrator,
                      cache::OrdinaryDiffEq.RosenbrockMutableCache,x,node...)
  i = length(integrator.u)
  cache.J = similar(cache.J,i,i)
  cache.W = similar(cache.W,i,i)
  add_node_jac_config!(cache.jac_config,i,x,node...)
  add_node_grad_config!(cache.grad_config,i,x,node...)
  nothing
end

function remove_node_non_user_cache!(integrator::DiffEqBase.AbstractODEIntegrator,
                      cache::OrdinaryDiffEq.RosenbrockMutableCache,node...)
  i = length(integrator.u)
  cache.J = similar(cache.J,i,i)
  cache.W = similar(cache.W,i,i)
  remove_node_jac_config!(cache.jac_config,i,node...)
  remove_node_grad_config!(cache.grad_config,i,node...)
  nothing
end

function add_node_jac_config!(cache::FiniteDiff.JacobianCache,i,x)
    #add_node!(cache.x1, fill!(similar(x, eltype(cache.x1)),0))
    add_node!(cache.fx, fill!(similar(x, eltype(cache.fx)),0))
    #cache.fx1 !== nothing && add_node!(cache.fx1, fill!(similar(x, eltype(cache.fx1)),0))
    cache.colorvec = 1:i
    nothing
end

function add_node_jac_config!(cache::FiniteDiff.JacobianCache,i,x,I...)
    #add_node!(cache.x1, fill!(similar(x, eltype(cache.x1)),0), I...)
    add_node!(cache.fx, fill!(similar(x, eltype(cache.fx)),0), I...)
    #cache.fx1 !== nothing && add_node!(cache.fx1, fill!(similar(x, eltype(cache.fx1)),0), I...)
    cache.colorvec = 1:i
    nothing
end

function remove_node_jac_config!(cache::FiniteDiff.JacobianCache,i,I...)
    #add_node!(cache.x1, fill!(similar(x, eltype(cache.x1)),0), I...)
    remove_node!(cache.fx, I...)
    #cache.fx1 !== nothing && add_node!(cache.fx1, fill!(similar(x, eltype(cache.fx1)),0), I...)
    cache.colorvec = 1:i
    nothing
end

function add_node_grad_config!(grad_config::ForwardDiff.DerivativeConfig,i,x,I...)
  resize!(grad_config.duals, i)
  grad_config
end

function add_node_grad_config!(grad_config::FiniteDiff.GradientCache,i,x,I...)
  grad_config.fx !== nothing && add_node!(grad_config.fx, fill!(similar(x, eltype(grad_config.fx)),0), I...)
  grad_config.c1 !== nothing && add_node!(grad_config.c1, fill!(similar(x, eltype(grad_config.c1)),0), I...)
  grad_config.c2 !== nothing && add_node!(grad_config.c2, fill!(similar(x, eltype(grad_config.c2)),0), I...)
  grad_config
end

function add_node_grad_config!(grad_config::FiniteDiff.GradientCache,i,x)
  grad_config.fx !== nothing && add_node!(grad_config.fx, fill!(similar(x, eltype(grad_config.fx)),0))
  grad_config.c1 !== nothing && add_node!(grad_config.c1, fill!(similar(x, eltype(grad_config.c1)),0))
  grad_config.c2 !== nothing && add_node!(grad_config.c2, fill!(similar(x, eltype(grad_config.c2)),0))
  grad_config
end

function remove_node_grad_config!(grad_config::FiniteDiff.GradientCache,i,I...)
  grad_config.fx !== nothing && remove_node!(grad_config.fx, I...)
  grad_config.c1 !== nothing && remove_node!(grad_config.c1, I...)
  grad_config.c2 !== nothing && remove_node!(grad_config.c2, I...)
  grad_config
end

function add_node_non_user_cache!(integrator::DiffEqBase.AbstractSDEIntegrator,idxs,x,node...)
  if DiffEqBase.is_diagonal_noise(integrator.sol.prob)
    add_node_noise!(integrator,idxs,x,node...)
    for c in rand_cache(integrator)
      add_node!(c,copy(x),node...)
    end
  end
end
function add_node_non_user_cache!(integrator::DiffEqBase.AbstractSDEIntegrator,idxs,x)
  if DiffEqBase.is_diagonal_noise(integrator.sol.prob)
    add_node_noise!(integrator,idxs,x)
    for c in rand_cache(integrator)
      add_node!(c,copy(x))
    end
  end
end

function add_node_noise!(integrator,idxs,x,node...)
  for c in integrator.W.S₁
    add_node!(c[2],copy(x),node...)
    if StochasticDiffEq.alg_needs_extra_process(integrator.alg)
      add_node!(c[3],copy(x),node...)
    end
    StochasticDiffEq.fill_new_noise_caches!(integrator,c,c[1],idxs)
  end
  for c in integrator.W.S₂
    add_node!(c[2],copy(x),node...)
    if StochasticDiffEq.alg_needs_extra_process(integrator.alg)
      add_node!(c[3],copy(x),node...)
    end
    StochasticDiffEq.fill_new_noise_caches!(integrator,c,c[1],idxs)
  end

  add_node!(integrator.W.dW,copy(x),node...)
  integrator.W.dW[idxs] .= zero(eltype(integrator.u))
  add_node!(integrator.W.curW,copy(x),node...)
  integrator.W.curW[idxs] .= zero(eltype(integrator.u))
  if StochasticDiffEq.alg_needs_extra_process(integrator.alg)
    add_node!(integrator.W.dZ,copy(x),node...)
    integrator.W.dZ[idxs] .= zero(eltype(integrator.u))
    add_node!(integrator.W.curZ,copy(x),node...)
    integrator.W.curZ[idxs] .= zero(eltype(integrator.u))
  end

  i = length(integrator.u)
  add_node!(integrator.W.dWtilde,copy(x),node...)
  add_node!(integrator.W.dWtmp,copy(x),node...)
  if StochasticDiffEq.alg_needs_extra_process(integrator.alg)
    add_node!(integrator.W.dZtmp,copy(x),node...)
    add_node!(integrator.W.dZtilde,copy(x),node...)
  end

  # fill in rands
  fill!(@view(integrator.W.curW[idxs]),zero(eltype(integrator.u)))
  if StochasticDiffEq.alg_needs_extra_process(integrator.alg)
    fill!(@view(integrator.W.curZ[idxs]),zero(eltype(integrator.u)))
  end
end

function add_node_noise!(integrator,idxs,x)
  for c in integrator.W.S₁
    add_node!(c[2],copy(x))
    if StochasticDiffEq.alg_needs_extra_process(integrator.alg)
      add_node!(c[3],copy(x))
    end
    StochasticDiffEq.fill_new_noise_caches!(integrator,c,c[1],idxs)
  end
  for c in integrator.W.S₂
    add_node!(c[2],copy(x))
    if StochasticDiffEq.alg_needs_extra_process(integrator.alg)
      add_node!(c[3],copy(x))
    end
    StochasticDiffEq.fill_new_noise_caches!(integrator,c,c[1],idxs)
  end

  add_node!(integrator.W.dW,copy(x))
  integrator.W.dW[idxs] .= zero(eltype(integrator.u))
  add_node!(integrator.W.curW,copy(x))
  integrator.W.curW[idxs] .= zero(eltype(integrator.u))
  if StochasticDiffEq.alg_needs_extra_process(integrator.alg)
    add_node!(integrator.W.dZ,copy(x))
    integrator.W.dZ[idxs] .= zero(eltype(integrator.u))
    add_node!(integrator.W.curZ,copy(x))
    integrator.W.curZ[idxs] .= zero(eltype(integrator.u))
  end

  i = length(integrator.u)
  add_node!(integrator.W.dWtilde,copy(x))
  add_node!(integrator.W.dWtmp,copy(x))
  if StochasticDiffEq.alg_needs_extra_process(integrator.alg)
    add_node!(integrator.W.dZtmp,copy(x))
    add_node!(integrator.W.dZtilde,copy(x))
  end

  # fill in rands
  fill!(@view(integrator.W.curW[idxs]),zero(eltype(integrator.u)))
  if StochasticDiffEq.alg_needs_extra_process(integrator.alg)
    fill!(@view(integrator.W.curZ[idxs]),zero(eltype(integrator.u)))
  end
end
