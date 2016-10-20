### Necessary until broadcast changes

function .*(m::AbstractMultiScaleModel,y::Number)
  new_m = similar(m)
  for i in eachindex(m)
    new_m[i] = m[i]*y
  end
  new_m
end

function .*(m::AbstractMultiScaleModel,m2::AbstractMultiScaleModel)
  new_m = similar(m)
  for i in eachindex(m)
    new_m[i] = m[i]*m2[i]
  end
  new_m
end

*(m::AbstractMultiScaleModel,y::AbstractMultiScaleModel) = m.*y
*(m::AbstractMultiScaleModel,y::Number) = m.*y
*(y::Number,m::AbstractMultiScaleModel) = m.*y

function .+(m::AbstractMultiScaleModel,y::Number)
  new_m = similar(m)
  for i in eachindex(m)
    new_m[i] = m[i]+y
  end
  new_m
end

function .+(m::AbstractMultiScaleModel,m2::AbstractMultiScaleModel)
  new_m = similar(m)
  for i in eachindex(m)
    new_m[i] = m[i]+m2[i]
  end
  new_m
end

+(m::AbstractMultiScaleModel,y::AbstractMultiScaleModel) = m.+y
+(m::AbstractMultiScaleModel,y::Number) = m.+y
+(y::Number,m::AbstractMultiScaleModel) = m.+y

function ./(m::AbstractMultiScaleModel,y::Number)
  new_m = similar(m)
  for i in eachindex(m)
    new_m[i] = m[i]/y
  end
  new_m
end

function ./(y::Number,m::AbstractMultiScaleModel)
  new_m = similar(m)
  for i in eachindex(m)
    new_m[i] = y/m[i]
  end
  new_m
end

function ./(m::AbstractMultiScaleModel,m2::AbstractMultiScaleModel)
  new_m = similar(m)
  for i in eachindex(m)
    new_m[i] = m[i]/m2[i]
  end
  new_m
end

/(m::AbstractMultiScaleModel,y::Number) = m./y
/(y::Number,m::AbstractMultiScaleModel) = y./m



function .-(m::AbstractMultiScaleModel,y::Number)
  new_m = similar(m)
  for i in eachindex(m)
    new_m[i] = m[i]-y
  end
  new_m
end

function .-(y::Number,m::AbstractMultiScaleModel)
  new_m = similar(m)
  for i in eachindex(m)
    new_m[i] = y-m[i]
  end
  new_m
end

function .-(m::AbstractMultiScaleModel,m2::AbstractMultiScaleModel)
  new_m = similar(m)
  for i in eachindex(m)
    new_m[i] = m[i]-m2[i]
  end
  new_m
end

-(m::AbstractMultiScaleModel,y::AbstractMultiScaleModel) = m.-y
-(m::AbstractMultiScaleModel,y::Number) = m.-y
-(y::Number,m::AbstractMultiScaleModel) = y.-m
