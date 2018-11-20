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
    deleteat_non_user_cache!(integrator, idxs) # required to do noise correctly
end

function add_node!(integrator::DiffEqBase.DEIntegrator, x, I...)
    cur_len = length(integrator.u)
    add_len = length(x)
    last_idx = length(integrator.u[I...].nodes)
    idx_start = getindices(integrator.u, last_idx)[end] + 1
    idxs = idx_start:idx_start+add_len-1
    for c in full_cache(integrator)
        add_node!(c, similar(x, eltype(c)), I...)
    end
    if DiffEqBase.is_diagonal_noise(integrator.sol.prob)
      for c in DiffEqBase.ratenoise_cache(integrator)
          add_node!(c, similar(x, eltype(c)), I...)
      end
    end
    addat_non_user_cache!(integrator, idxs) # required to do noise correctly
end

function add_node!(integrator::DiffEqBase.DEIntegrator, x)
    cur_len = length(integrator.u)
    add_len = length(x)
    last_idx = length(integrator.u.nodes)
    idx_start = getindices(integrator.u, last_idx)[end] + 1
    idxs = idx_start:idx_start+add_len-1
    @show idxs
    for c in full_cache(integrator)
        add_node!(c, similar(x, eltype(c)))
    end
    if DiffEqBase.is_diagonal_noise(integrator.sol.prob)
      for c in DiffEqBase.ratenoise_cache(integrator)
          add_node!(c, similar(x, eltype(c)))
      end
    end
    addat_non_user_cache!(integrator, idxs) # required to do noise correctly
end


reshape(m::AbstractMultiScaleArray, i::Int...) = m
