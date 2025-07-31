length(m::AbstractMultiScaleArrayLeaf) = length(m.values)
length(m::AbstractMultiScaleArray) = m.end_idxs[end]
Base.isempty(m::AbstractMultiScaleArray) = isempty(m.nodes) && isempty(m.values)
Base.isempty(m::AbstractMultiScaleArrayLeaf) = isempty(m.values)
num_nodes(m::AbstractMultiScaleArrayLeaf) = 0
num_nodes(m::AbstractMultiScaleArray) = size(m.nodes, 1)
ndims(m::AbstractMultiScaleArray) = 1
size(m::AbstractMultiScaleArray, i::Int) = i == 1 ? length(m) : 0
size(m::AbstractMultiScaleArray) = (length(m),)

parameterless_type(T::Type) = Base.typename(T).wrapper
parameterless_type(x) = parameterless_type(typeof(x))

@generated function similar(m::AbstractMultiScaleArrayLeaf, ::Type{T} = eltype(m)) where {T}
    assignments = [s == :values ? :(similar(m.values, T)) :
                   (sq = Meta.quot(s); :(deepcopy(getfield(m, $sq))))
                   for s in fieldnames(m)[2:end]] # 1 is values
    :(construct(parameterless_type(m), similar(m.values, T), $(assignments...)))
end

@generated function similar(m::AbstractMultiScaleArray, ::Type{T} = eltype(m)) where {T}
    assignments = [s == :values ? :(similar(m.values, T)) :
                   (sq = Meta.quot(s); :(deepcopy(getfield(m, $sq))))
                   for s in fieldnames(m)[4:end]] # 1:3 is nodes,values,end_idxs
    :(construct(parameterless_type(m), recursive_similar(m.nodes, T), similar(m.values, T),
        $(assignments...)))
end

Base.zero(A::AbstractMultiScaleArray) = fill!(similar(A), 0)

recursive_similar(x, T) = [similar(y, T) for y in x]
recursive_similar(x::Tuple, T) = tuple((similar(y, T) for y in x)...)

construct(::Type{T}, args...) where {T <: AbstractMultiScaleArrayLeaf} = T(args...)

function __construct(T, nodes, values, args...)
    vallen = length(values)
    end_idxs = Vector{Int}(undef, length(nodes) + ifelse(vallen == 0, 0, 1))
    off = 0
    @inbounds for i in 1:length(nodes)
        end_idxs[i] = (off += length(nodes[i]))
    end
    vallen == 0 || (end_idxs[end] = off + vallen)
    T(nodes, values, end_idxs, args...)
end

(construct(::Type{T},
    nodes::AbstractVector{<:AbstractMultiScaleArray},
    args...)
    where {T <: AbstractMultiScaleArray}) = __construct(T, nodes, eltype(T)[], args...)

(construct(::Type{T}, nodes::AbstractVector{<:AbstractMultiScaleArray}, values,
    args...)
    where {T <: AbstractMultiScaleArray}) = __construct(T, nodes, values, args...)

(construct(::Type{T},
    nodes::Tuple{Vararg{AbstractMultiScaleArray}},
    args...)
    where {T <: AbstractMultiScaleArray}) = __construct(T, nodes, eltype(T)[], args...)

(construct(::Type{T}, nodes::Tuple{Vararg{AbstractMultiScaleArray}}, values,
    args...)
    where {T <: AbstractMultiScaleArray}) = __construct(T, nodes, values, args...)

function vcat(m1::AbstractMultiScaleArray, m2::AbstractMultiScaleArray)
    error("AbstractMultiScaleArrays cannot be concatenated")
end

function hcat(m1::AbstractMultiScaleArray, m2::AbstractMultiScaleArray)
    error("AbstractMultiScaleArrays cannot be concatenated")
end

==(m1::AbstractMultiScaleArray, m2::AbstractMultiScaleArray) = (m1 === m2)

function recursivecopy!(b::AbstractMultiScaleArrayLeaf, a::AbstractMultiScaleArrayLeaf)
    @inbounds copyto!(b, a)
end

function recursivecopy!(b::AbstractMultiScaleArray, a::AbstractMultiScaleArray)
    @inbounds for i in eachindex(a.nodes)
        recursivecopy!(b.nodes[i], a.nodes[i])
    end
    recursivecopy!(b.values, a.values)
end

function recursivecopy(a::AbstractMultiScaleArray)
    out = similar(a)
    recursivecopy!(out, a)
    out
end
