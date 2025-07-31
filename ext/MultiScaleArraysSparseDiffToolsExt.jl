module MultiScaleArraysSparseDiffToolsExt

using MultiScaleArrays, SparseDiffTools

function MultiScaleArrays.add_node_jac_config!(
        cache, config::SparseDiffTools.ForwardColorJacCache, i, x,
        node...)
    @assert cache.jac_config.colorvec isa UnitRange
    cache.jac_config = SparseDiffTools.ForwardColorJacCache(cache.uf, cache.uprev,
        config.chunksize)
    nothing
end

function MultiScaleArrays.add_node_jac_config!(cache, config::SparseDiffTools.ForwardColorJacCache, i, x)
    @assert cache.jac_config.colorvec isa UnitRange
    cache.jac_config = SparseDiffTools.ForwardColorJacCache(cache.uf, cache.uprev,
        config.chunksize)
    nothing
end

function MultiScaleArrays.remove_node_jac_config!(
        cache, config::SparseDiffTools.ForwardColorJacCache, i, x,
        node...)
    @assert cache.jac_config.colorvec isa UnitRange
    cache.jac_config = SparseDiffTools.ForwardColorJacCache(cache.uf, cache.uprev,
        config.chunksize)
    nothing
end

function MultiScaleArrays.remove_node_jac_config!(cache, config::SparseDiffTools.ForwardColorJacCache, i, x)
    @assert cache.jac_config.colorvec isa UnitRange
    cache.jac_config = SparseDiffTools.ForwardColorJacCache(cache.uf, cache.uprev,
        config.chunksize)
    nothing
end

end
