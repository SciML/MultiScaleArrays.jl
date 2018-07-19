function remove_node!(integrator::DiffEqBase.DEIntegrator, I...)
    idxs = getindices(integrator.u, I...)
    for c in user_cache(integrator)
        remove_node!(c, I...)
    end
    deleteat_non_user_cache!(integrator, idxs)
end

function add_node!(integrator::DiffEqBase.DEIntegrator, x, I...)
    cur_len = length(integrator.u)
    add_len = length(x)
    for c in user_cache(integrator)
        add_node!(c, similar(x, eltype(c)), I...)
    end
    last_idx = length(integrator.u[I...].nodes)
    idxs = getindices(integrator.u, I..., last_idx)
    addat_non_user_cache!(integrator, idxs)
end

reshape(m::AbstractMultiScaleArray, i::Int...) = m
