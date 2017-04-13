length(m::AbstractMultiScaleArrayLeaf) = length(m.x)
length(m::AbstractMultiScaleArray) = m.end_idxs[end]
num_daughters(m::AbstractMultiScaleArrayLeaf) = 0
num_daughters(m::AbstractMultiScaleArray) = size(m.x)[1]
ndims(m::AbstractMultiScaleArray) = 1
size(m::AbstractMultiScaleArray,i::Int) = Int((i==1))*length(m)
size(m::AbstractMultiScaleArray) = (length(m),)


function similar(m::AbstractMultiScaleArray)
  m_new = construct(typeof(m),deepcopy(m.x),deepcopy(m.y))
  #for i in eachindex(m.x) # Is this necessary?
  #  m_new.x[i] = similar(m.x[i])
  #end
  m_new
end

similar(m::AbstractMultiScaleArrayLeaf) = construct(typeof(m),similar(m.x))

similar(m::AbstractMultiScaleArrayLeaf,T::Type) = construct(typeof(m).name.primary,similar(m.x,T))

function similar(m::AbstractMultiScaleArray,T::Type)
  new_x = [similar(v,T) for v in m.x]
  new_y = similar(m.y,T)
  construct(typeof(m).name.primary,new_x,new_y)
end

function construct{T<:AbstractMultiScaleArrayLeaf}(::Type{T},x)
  T(x)
end

function construct{T<:AbstractMultiScaleArray,T2<:AbstractMultiScaleArray}(::Type{T},x::Vector{T2},y=Float64[])
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

function vcat(m1::AbstractMultiScaleArray,m2::AbstractMultiScaleArray)
  error("AbstractMultiScaleArrays cannot be concatenated")
end

function hcat(m1::AbstractMultiScaleArray,m2::AbstractMultiScaleArray)
  error("AbstractMultiScaleArrays cannot be concatenated")
end

function ==(m1::AbstractMultiScaleArray,m2::AbstractMultiScaleArray)
  m1 === m2
end

function recursivecopy!(b::AbstractMultiScaleArrayLeaf,a::AbstractMultiScaleArrayLeaf)
  @inbounds copy!(b,a)
end

function recursivecopy!(b::AbstractMultiScaleArray,a::AbstractMultiScaleArray)
  @inbounds for i in eachindex(a.x)
    recursivecopy!(b.x[i],a.x[i])
  end
  @inbounds for i in eachindex(a.y)
    recursivecopy!(b.y[i],a.y[i])
  end
end
