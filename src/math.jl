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
