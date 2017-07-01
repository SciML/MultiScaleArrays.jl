#=
function Base.broadcast!(f, A::AbstractMultiScaleArray, B::AbstractMultiScaleArray, N::Number)
    for i in eachindex(A.x)
        broadcast!(f, A.x[i], B.x[i], N)
    end
    A
end

function Base.broadcast!(f, A::AbstractMultiScaleArray, N::Number, B::AbstractMultiScaleArray)
    for i in eachindex(A.x)
        broadcast!(f, A.x[i], N, B.x[i])
    end
    A
end

function Base.broadcast!(f, A::AbstractMultiScaleArray,
                         B::AbstractMultiScaleArray, C::AbstractMultiScaleArray)
    for i in eachindex(A.x)
        broadcast!(f, A.x[i], B.x[i], C.x[i])
    end
    A
end

function Base.broadcast(f, B::AbstractMultiScaleArray, N::Number)
    A = similar(B)
    for i in eachindex(A.x)
        broadcast!(f, A.x[i], B.x[i], N)
    end
    A
end

function Base.broadcast(f, N::Number, B::AbstractMultiScaleArray)
    A = similar(B)
    for i in eachindex(A.x)
        broadcast!(f, A.x[i], N, B.x[i])
    end
    A
end

function Base.broadcast(f, B::AbstractMultiScaleArray, C::AbstractMultiScaleArray)
    A = similar(B)
    for i in eachindex(A.x)
        broadcast!(f, A.x[i], B.x[i], C.x[i])
    end
    A
end
=#

+(m::AbstractMultiScaleArray, y::AbstractMultiScaleArray) = m .+ y
+(m::AbstractMultiScaleArray, y::Number) = m .+ y
+(y::Number, m::AbstractMultiScaleArray) = m .+ y

-(m::AbstractMultiScaleArray, y::AbstractMultiScaleArray) = m .- y
-(m::AbstractMultiScaleArray, y::Number) = m .- y
-(y::Number, m::AbstractMultiScaleArray) = y .- m

*(m::AbstractMultiScaleArray, y::Number) = m .* y
*(y::Number, m::AbstractMultiScaleArray) = m .* y

/(m::AbstractMultiScaleArray, y::Number) = m ./ y
/(y::Number, m::AbstractMultiScaleArray) = y ./ m
