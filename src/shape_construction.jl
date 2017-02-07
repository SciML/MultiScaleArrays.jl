length(m::MultiScaleArrayLeaf) = length(m.x)
length(m::AbstractMultiScaleArray) = m.end_idxs[end]
num_daughters(m::MultiScaleArrayLeaf) = 0
num_daughters(m::AbstractMultiScaleArray) = size(m.x)[1]
ndims(m::AbstractMultiScaleArray) = 1
size(m::AbstractMultiScaleArray,i::Int) = Int((i==1))*length(m)
size(m::AbstractMultiScaleArray) = (length(m),)


function similar(m::AbstractMultiScaleArray)
  m_new = construct(typeof(m),deepcopy(m.x),copy(m.y)) # Copy because y is a vector!
  #for i in eachindex(m.x) # Is this necessary?
  #  m_new.x[i] = similar(m.x[i])
  #end
  m_new
end

similar(m::MultiScaleArrayLeaf) = construct(typeof(m),similar(m.x))

#=
function print_matrix(IO,m::AbstractMultiScaleArray,str1,str2,str3)
  print_matrix(IO,m.x,str1,str2,str3)
end
=#

function construct{T<:MultiScaleArrayLeaf,T2}(::Type{T},x::Vector{T2})
  T(x)
end

function construct{T<:AbstractMultiScaleArray,T2<:AbstractMultiScaleArray,T3<:Number}(::Type{T},x::Vector{T2},y::Vector{T3}=Float64[])
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

function recursivecopy!(b::MultiScaleArrayLeaf,a::MultiScaleArrayLeaf)
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
