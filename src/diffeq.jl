function remove_daughter!(integrator::DEIntegrator,I...)
  cur_len = length(integrator.u)
  remove_len = length(integrator.u[I...])
  for c in user_cache(integrator)
    remove_daughter!(c,I...)
  end
  resize_non_user_cache!(integrator,cur_len-remove_len)
end

function add_daughter!(integrator::DEIntegrator,x,I...)
  cur_len = length(integrator.u)
  add_len = length(x)
  for c in user_cache(integrator)
    add_daughter!(c,deepcopy(x),I...)
  end
  resize_non_user_cache!(integrator,cur_len+add_len)
end

reshape(m::AbstractMultiScaleArray,i::Int...) = m
