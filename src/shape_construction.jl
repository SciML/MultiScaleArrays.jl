length(m::MultiScaleModelLeaf) = length(m.x)
length(m::AbstractMultiScaleModel) = m.end_idxs[end]
num_daughters(m::MultiScaleModelLeaf) = 0
num_daughters(m::AbstractMultiScaleModel) = size(m.x)[1]
ndims(m::AbstractMultiScaleModel) = 1
size(m::AbstractMultiScaleModel,i::Int) = Int((i==1))*length(m)
size(m::AbstractMultiScaleModel) = (length(m),)


function similar(m::AbstractMultiScaleModel)
  m_new = construct(typeof(m),deepcopy(m.x),copy(m.y)) # Copy because y is a vector!
  for i in eachindex(m.x) # Is this necessary?
    m_new.x[i] = similar(m.x[i])
  end
  m_new
end

similar(m::MultiScaleModelLeaf) = construct(typeof(m),similar(m.x))

#=
function print_matrix(IO,m::AbstractMultiScaleModel,str1,str2,str3)
  print_matrix(IO,m.x,str1,str2,str3)
end
=#

function construct{T<:MultiScaleModelLeaf,T2}(::Type{T},x::Vector{T2})
  T(x)
end

function construct{T<:AbstractMultiScaleModel,T2<:AbstractMultiScaleModel,T3<:Number}(::Type{T},x::Vector{T2},y::Vector{T3})
  end_idxs = Vector{Int}(length(x))
  end_idxs[1] = length(x[1])
  for i in 2:length(x)
    end_idxs[i] = end_idxs[i-1] + length(x[i])
  end
  if !isempty(y)
    push!(end_idxs,end_idxs[end] + length(y))
  end
  m = T(x,y,end_idxs)
end

construct{T<:AbstractMultiScaleModel,T2<:AbstractMultiScaleModel}(::Type{T},x::Vector{T2}) = construct(T,x,Float64[])

function vcat(m1::AbstractMultiScaleModel,m2::AbstractMultiScaleModel)
  error("AbstractMultiScaleModels cannot be concatenated")
end

function hcat(m1::AbstractMultiScaleModel,m2::AbstractMultiScaleModel)
  error("AbstractMultiScaleModels cannot be concatenated")
end

function ==(m1::AbstractMultiScaleModel,m2::AbstractMultiScaleModel)
  m1 === m2
end

function recursivecopy!(b::MultiScaleModelLeaf,a::MultiScaleModelLeaf)
  @inbounds copy!(b,a)
end

function recursivecopy!(b::AbstractMultiScaleModel,a::AbstractMultiScaleModel)
  @inbounds for i in eachindex(a.x)
    recursivecopy!(b.x[i],a.x[i])
  end
  @inbounds for i in eachindex(a.y)
    recursivecopy!(b.y[i],a.y[i])
  end
end
