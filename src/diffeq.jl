function remove_daughter!(integrator::DEIntegrator,I...)
  remove_len = length(integrator.u[I...])
  for c in user_cache(integrator)
    remove_daughter!(c,I...)
  end
  for c in full_cache(integrator)
    if !(typeof(c) <: AbstractMultiScaleArray)
      resize!(c,length(c)-remove_len)
    end
  end
end

function add_daughter!(integrator::DEIntegrator,x,I...)
  for c in user_cache(integrator)
    add_daughter!(c,deepcopy(x),I...)
  end
  add_len = length(x)
  for c in full_cache(integrator)
    if !(typeof(c) <: AbstractMultiScaleArray)
      resize!(c,length(c)+add_len)
    end
  end
end
