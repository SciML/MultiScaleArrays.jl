### Necessary until broadcast changes

if VERSION < v"0.6-"
  function .*(m::AbstractMultiScaleArray,y::Number)
    new_m = similar(m)
    for i in eachindex(m)
      new_m[i] = m[i]*y
    end
    new_m
  end

  function .*(m::AbstractMultiScaleArray,m2::AbstractMultiScaleArray)
    new_m = similar(m)
    for i in eachindex(m)
      new_m[i] = m[i]*m2[i]
    end
    new_m
  end

  function .+(m::AbstractMultiScaleArray,y::Number)
    new_m = similar(m)
    for i in eachindex(m)
      new_m[i] = m[i]+y
    end
    new_m
  end

  function .+(m::AbstractMultiScaleArray,m2::AbstractMultiScaleArray)
    new_m = similar(m)
    for i in eachindex(m)
      new_m[i] = m[i]+m2[i]
    end
    new_m
  end

  function ./(m::AbstractMultiScaleArray,y::Number)
    new_m = similar(m)
    for i in eachindex(m)
      new_m[i] = m[i]/y
    end
    new_m
  end

  function ./(y::Number,m::AbstractMultiScaleArray)
    new_m = similar(m)
    for i in eachindex(m)
      new_m[i] = y/m[i]
    end
    new_m
  end

  function ./(m::AbstractMultiScaleArray,m2::AbstractMultiScaleArray)
    new_m = similar(m)
    for i in eachindex(m)
      new_m[i] = m[i]/m2[i]
    end
    new_m
  end

  function .-(m::AbstractMultiScaleArray,y::Number)
    new_m = similar(m)
    for i in eachindex(m)
      new_m[i] = m[i]-y
    end
    new_m
  end

  function .-(y::Number,m::AbstractMultiScaleArray)
    new_m = similar(m)
    for i in eachindex(m)
      new_m[i] = y-m[i]
    end
    new_m
  end

  function .-(m::AbstractMultiScaleArray,m2::AbstractMultiScaleArray)
    new_m = similar(m)
    for i in eachindex(m)
      new_m[i] = m[i]-m2[i]
    end
    new_m
  end

else # Only define broadcasts on v0.6+
  #=
  add_idxs(x,expr) = expr
  add_idxs{T<:AbstractMultiScaleArray}(::Type{T},expr) = :($(expr).x[i])

  add_y(x,expr) = expr
  add_y{T<:AbstractMultiScaleArray}(::Type{T},expr) = :($(expr).y)

  @generated function Base.broadcast!(f,A::AbstractMultiScaleArray,B::Union{Number,AbstractMultiScaleArray}...)
    exs = ((add_idxs(B[i],:(B[$i])) for i in eachindex(B))...)
    exs_y = ((add_y(B[i],:(B[$i])) for i in eachindex(B))...)
    quote
      for i in eachindex(A.x)
        broadcast!(f,A.x[i],$(exs...))
      end
      if !(typeof(A)<:AbstractMultiScaleArrayLeaf) && !isempty(y)
        broadcast!(f,A.y,$(exs_y...))
      end
    end
  end
  =#
end

*(m::AbstractMultiScaleArray,y::Number) = m.*y
*(y::Number,m::AbstractMultiScaleArray) = m.*y
+(m::AbstractMultiScaleArray,y::AbstractMultiScaleArray) = m.+y
+(m::AbstractMultiScaleArray,y::Number) = m.+y
+(y::Number,m::AbstractMultiScaleArray) = m.+y

/(m::AbstractMultiScaleArray,y::Number) = m./y
/(y::Number,m::AbstractMultiScaleArray) = y./m

-(m::AbstractMultiScaleArray,y::AbstractMultiScaleArray) = m.-y
-(m::AbstractMultiScaleArray,y::Number) = m.-y
-(y::Number,m::AbstractMultiScaleArray) = y.-m
