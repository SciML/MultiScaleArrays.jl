import Base.getindex
import Base.length
import Base.==
import Base.println
import Unicode.isnumeric
import Base.isnan
# SUPR
# import Base.string
# string(::Nothing) = "nothing"
function isnan(::String) false end

stringNothingInVec(a::Vector) = replace(string([x==nothing ? "nothing" : x for x in a]),  "\"nothing\""=>"nothing")
stringNothingInVec(a::Vector) = replace(replace(replace(replace(string([x==nothing ? "nothing" : x for x in a]),  "\"nothing\""=>"nothing"), "Array{Float64,1}"=>""), "Any"=>""),"Array{Int64,1}"=>"")
stringNothingInVec(a::Nothing) = "nothing"

function isnumeric(s::AbstractString)
    tryparse(Float64, s) isa Number
end

function isinside(val::Union{Number,Symbol,Char},vect)
    if isa(vect,Number) | isa(vect,Symbol) | isa(vect,Char)
        vect = [vect]
    end
    out::Bool = false
    i=0
    while !(out | (i == length(vect)))
        i += 1 
        out = val == vect[i]
    end
    return out
end

function printTable(c::Array)
    a = (println(c[i,:]) for i in 1:size(c)[1])
    for aa in a end
end

function printTable(c::Tuple)
    a = (println(c[i]) for i in eachindex(c)[1])
    for aa in a end
end

function seq(min,max,len) ; collect(min:(max-min)/len:max) end


# flat : equivalent of the R function "unlist"
function grep!(v::Array,rst::Vector{Any})
    for x in v
        if isa(x, Array) | isa(x, Tuple) | isa(x, NamedTuple) grep!(x,rst) else push!(rst, x) end
end;end
function flat(arr::Array)
    rst = Any[]
    grep!(arr,rst)
    rst
end

function grep!(v::NamedTuple,rst::Vector{Any})
    for x in v
        if isa(x, Array) | isa(x, Tuple) | isa(x, NamedTuple) grep!(x,rst) else push!(rst, x) end
end;end
function flat(arr::NamedTuple)
    rst = Any[]
    grep!(arr,rst)
    rst
end

function grep!(v::Tuple,rst::Vector{Any})
    for x in v
        if isa(x, Array) | isa(x, Tuple) | isa(x, NamedTuple) grep!(x,rst) else push!(rst, x) end
end;end

function flat(arr::Tuple)
    rst = Any[]
    grep!(arr,rst)
    rst
end

MaxIntValues = (Int8=Int8(2)^7-1, Int16=Int16(2)^15-1, Int32=(2)^31-1, Int64=Int64(2)^63-1, Int128=Int128(2)^127-1, BigInt=Inf)

function Int_(x::Integer,Max::Integer)
    if (Max < (2^7-1)) return(Int8(x)) elseif (Max < (2^15-1)) return(Int16(x)) elseif (Max < (2^31-1)) return(Int32(x)) elseif (Max < (2^63-1)) return(Int64(x)) elseif (Max < (2^127-1))  return(Int128(x)) else return(BigInt(x)) end
end
Int_(x::Integer) = Int_(x,x)

function Int_(x::Vector{<:Integer})
    Max = max(maximum(x),abs(minimum(x)))
    if (Max < (2^7-1)) T = Int8 elseif (Max < (2^15-1)) T = Int16 elseif (Max < (2^31-1)) T = Int32 elseif (Max < (2^63-1)) T = Int64 elseif (Max < (2^127-1))  T = Int128 else T = BigInt end
    return(T.(x))
end

function Int_(x::NTuple{N,<:Integer} where N)
    Max = max(maximum(x),abs(minimum(x)))
    if (Max < (2^7-1)) T = Int8 elseif (Max < (2^15-1)) T = Int16 elseif (Max < (2^31-1)) T = Int32 elseif (Max < (2^63-1)) T = Int64 elseif (Max < (2^127-1))  T = Int128 else T = BigInt end
    return(tuple(T.(x)...))
end

function Int_(x::UnitRange{<:Integer})
    Max = max(maximum(x),abs(minimum(x)))
    if (Max < (2^7-1)) T = Int8 elseif (Max < (2^15-1)) T = Int16 elseif (Max < (2^31-1)) T = Int32 elseif (Max < (2^63-1)) T = Int64 elseif (Max < (2^127-1))  T = Int128 else T = BigInt end
    Max = maximum(x)
    Min = minimum(x)
    return(T(Min):T(Max))
end

function Int_(x::StepRange{<:Integer,<:Integer})
    Max = maximum(x)
    if (Max < (2^7-1)) T = Int8 elseif (Max < (2^15-1)) T = Int16 elseif (Max < (2^31-1)) T = Int32 elseif (Max < (2^63-1)) T = Int64 elseif (Max < (2^127-1))  T = Int128 else T = BigInt end
    xx = tryparse.(Int128,split(string(x),":"))
    return(T(xx[1]):T(xx[2]):T(xx[3]))
end

function Int_(x::Integer,T::Type)
    return(T(x))
end

function Int_(x::Vector{<:Integer},T::Type)
    return(T.(x))
end

function Int_(x::NTuple{N,<:Integer} where N,T::Type)
    return(tuple(T.(x)...))
end

function Int_(x::UnitRange{<:Integer},T::Type)
    Max = maximum(x)
    Min = minimum(x)
    return(T(Min):T(Max))
end

function Int_(x::StepRange{<:Integer,<:Integer},T::Type)
    Max = maximum(x)
    xx = tryparse.(Int,split(string(x),":"))
    return(T(xx[1]):T(xx[2]):T(xx[3]))
end


# SUPR
# function F0_1_to_0_∞(x) 1/(-x+1)-1 end

function Power10__0_1_to_min_max(x,Min,Max)
    (((x-Min)/(Max-Min))^10)*(Max-Min)+Min
end

function POW10_01_MinMax(Min,Max,N;Zero=true,integer=false)
    if (!Zero & integer) error("In the way it is implemented, this function does not work if Zero==false & integer==true)") end
    if Zero
        out = Power10__0_1_to_min_max.(seq(Min,Max,N-1),Min,Max)[1: end   ]
    else
        out = Power10__0_1_to_min_max.(seq(Min,Max,N  ),Min,Max)[2:(end-1)]
    end
    if integer
        return Int.(round.(out))
    else
        return out
    end
end



# SUPR!
# function SetExponetialInteger(Nalleles,Max,base)
#     ExponetialInteger = []
#     n=Nalleles
#     while length(ExponetialInteger)<Nalleles
#         ExponetialInteger = unique(round.(base.^(seq(0,log(base,Max),(n-1)))))
#         n+=1
#     end
#     return ExponetialInteger
# end


################################# MeanBeta
@doc """MeanBeta(alpha, beta): Get the mean of a Beta distribution from its alpha and beta parameters
# Examples (uniform distribution)
```jldoctest
julia> a = 0.5 ; b = 0.5
julia> MeanBeta(a,b)
0.5
```
# See also:
```jldoctest
VarBeta
GetBeta_a
GetBeta_b
```
""" -> 
function MeanBeta(a,b)
    a/(a+b)
end

################################# VarBeta
@doc """VarBeta(alpha, beta): Get the variance of a Beta distribution from its alpha and beta parameters
# Examples (uniform distribution)
```jldoctest
julia> a = 0.5 ; b = 0.5
julia> VarBeta(a,b)
0.125
```
# See also:
```jldoctest
MeanBeta
GetBeta_a
GetBeta_b
```
""" -> 
function VarBeta(a,b)
    a*b/((a+b)^2*(a+b+1))
end

@doc """GetBeta_a(Mean, Variance): Get the alpha of a Beta distribution from its mean and variance""" ->
function GetBeta_a(M,V) 1*(-M^3+M^2-M*V)/V end

@doc """GetBeta_b(Mean, Variance): Get the beta of a Beta distribution from its mean and variance""" ->
function GetBeta_b(M,V) 1*(1/M-1)*(-M^3+M^2-M*V)/V end

# a,b = rand(2)
# 
# a - GetBeta_a(MeanBeta(a,b) , VarBeta(a,b))
# b - GetBeta_b(MeanBeta(a,b) , VarBeta(a,b))

function CheckMeanVarOfBeta(ShapeOfTheCostOfBeingAggressed)
    a = GetBeta_a(ShapeOfTheCostOfBeingAggressed[:Mean],ShapeOfTheCostOfBeingAggressed[:Var]) ; if (a<=0) error("This choice of mean and variance are not possible for the Beta distribution.\nChoose a lower value for the variance or a less extrem value for the mean") end
    b = GetBeta_a(ShapeOfTheCostOfBeingAggressed[:Mean],ShapeOfTheCostOfBeingAggressed[:Var]) ; if (a<=0) error("This choice of mean and variance are not possible for the Beta distribution.\nChoose a lower value for the variance or a less extrem value for the mean") end
    return("Ok, your choice of the mean and variance of CostOfBeingAggressed can be achieved with the Beta distribution !")
end


function AddKeyToDict!(d,key,val,KeyCheck=:Iam)
    for k in keys(d)
        if typeof(d[k])<:Dict
            AddKeyToDict!(d[k],key,val,KeyCheck)
        end
    end
    if typeof(key)<:Symbol & any(keys(d).==KeyCheck)
        d[key] = val
    else
        (length(key)==length(val)) ? nothing : error("!length(key)==length(val)")
        for i in eachindex(key)
            if any(keys(d[i]).==KeyCheck)
                d[key[i]] = val[i]
end;end;end;end


function AddKeyToDict(d,key,val,KeyCheck=:Iam)
    d=deepcopy(d)
    for k in keys(d)
        if typeof(d[k])<:Dict
            AddKeyToDict!(d[k],key,val,KeyCheck)
        end
    end
    if typeof(key)==Symbol & any(keys(d).==KeyCheck)
        d[key] = val
    else
        (length(key)==length(val)) ? nothing : error("!length(key)==length(val)")
        for i in eachindex(key)  
            if any(keys(dkey[i]).==KeyCheck) & typeof(key[i])==Symbol # Exclud keys that are strings
                d[key[i]] = val[i]
            end
        end
    end
    return d
end

function ChangeValuesInDict!(d,val,rm=nothing)
    convert(Dict{Any,Any}, d)
    # remove rm
    if typeof(rm)<:Symbol
        delete!(d,rm)
    else
        for i in eachindex(rm)
            delete!(d,rm[i])
    end;end
    # Change other values to val
    # for k in filter(x -> typeof(x)<:Symbol,keys(d))
    for k in keys(d)
        if typeof(d[k])<:Dict
            ChangeValuesInDict!(d[k],val,rm)
        else
            d[k] = val
end;end;end

function ChangeValuesInDict(d,val,rm=nothing)
    d=deepcopy(d)
    convert(Dict{Any,Any}, d)
    # remove rm
    if typeof(rm)<:Symbol
        delete!(d,rm)
    else
        for i in eachindex(rm)
            delete!(d,rm[i])
    end;end
    # Change other values to val
    # for k in filter(x -> typeof(x)<:Symbol,keys(d))
    for k in keys(d)
        if typeof(d[k])<:Dict
            ChangeValuesInDict!(d[k],val,rm)
        else
            d[k] = val
        end
    end
    return d
end

function ChangeValuesInDict!(d,val)
    convert(Dict{Any,Any}, d)
    # Change other values to val
    # for k in filter(x -> typeof(x)<:Symbol,keys(d))
    for k in keys(d)
        if typeof(d[k])<:Dict
            ChangeValuesInDict!(d[k],val)
        else
            d[k] = val
end;end;end

function ChangeValuesInDict(d,val)
    d=deepcopy(d)
    convert(Dict{Any,Any}, d)
    # Change other values to val
    # for k in filter(x -> typeof(x)<:Symbol,keys(d))
    for k in keys(d)
        if typeof(d[k])<:Dict
            ChangeValuesInDict!(d[k],val)
        else
            d[k] = val
        end
    end
    return d
end


function DuplicateTraitsValuesINIValuesForEachSp!(d::Dict,TraitAffectHostsSymbiontsCells,TraitsValues,MuMu,SimParam)
    Ref = string.(keys(TraitAffectHostsSymbiontsCells))
    TraitsValues_, TraitAffectHostsSymbiontsCells_, D, MuMu_ = Dict(), Dict(), Dict(), Dict()
    for k in keys(d)
#       println(k)
        if occursin("_",string(k)) | occursin("!",string(k))
            ref = split(string(k),"_")[end]
            ref = split(ref,"!")[1]
            if !any(ref .== Ref) error("!any(ref .== Ref)") else ref = Symbol(ref) end
        else
            ref = k
        end
#       println("ifelse 1    ", TraitAffectHostsSymbiontsCells[k])
        if (any(TraitAffectHostsSymbiontsCells[ref] .== :OnlyHostInd) | (TraitAffectHostsSymbiontsCells[ref] == [:HostCells,:HostIndOrCell])) | ((sort(collect(TraitAffectHostsSymbiontsCells[ref])) == sort([:HostCells,:HostIndOrCell,:TraitByComp])) | (sort(collect(TraitAffectHostsSymbiontsCells[ref])) == sort(collect((:HostCells,:HostInd,:HostIndOrCell))))) | (sort(collect(TraitAffectHostsSymbiontsCells[ref])) == sort([:HostCells,:HostInd,:HostIndOrCell,:TraitByComp]))
            # host ony
#             println("_if 1")
                TraitAffectHostsSymbiontsCells_[Symbol(string(k)*"_sp"*string(1))] = TraitAffectHostsSymbiontsCells[ref]
                                              D[Symbol(string(k)*"_sp"*string(1))] = d[k]
                                  TraitsValues_[Symbol(string(k)*"_sp"*string(1))] = TraitsValues[ref]
                if occursin("Mu°",string(ref))   MuMu_[Symbol(string(k)*"_sp"*string(1))] = MuMu[ref]  end
        elseif (TraitAffectHostsSymbiontsCells[ref] == [:Symbiont]) | (TraitAffectHostsSymbiontsCells[ref] == [:Symbiont,:TraitByComp])
            # Symbiont only
#             println("_elseif 2")
            for sp in 2:SimParam.Nsp
                TraitAffectHostsSymbiontsCells_[Symbol(string(k)*"_sp"*string(sp))] = TraitAffectHostsSymbiontsCells[ref]
                                              D[Symbol(string(k)*"_sp"*string(sp))] = d[k]
                                  TraitsValues_[Symbol(string(k)*"_sp"*string(sp))] = TraitsValues[ref]
                if occursin("Mu°",string(ref))   MuMu_[Symbol(string(k)*"_sp"*string(sp))] = MuMu[ref] end
            end
            d[k] = [nothing,[deepcopy(d[k]) for i in 1:SimParam.NsymbiontSp]...]
        elseif (any(TraitAffectHostsSymbiontsCells[ref] .== :Symbiont) & any(TraitAffectHostsSymbiontsCells[ref] .== :HostCells)) | (any(TraitAffectHostsSymbiontsCells[ref] .== :Symbiont) & any(TraitAffectHostsSymbiontsCells[ref] .== :HostInd))
            # all
#             println("_elseif 3")
            for sp in 1:SimParam.Nsp
                TraitAffectHostsSymbiontsCells_[Symbol(string(k)*"_sp"*string(sp))] = TraitAffectHostsSymbiontsCells[ref]
                                              D[Symbol(string(k)*"_sp"*string(sp))] = d[k]
                                  TraitsValues_[Symbol(string(k)*"_sp"*string(sp))] = TraitsValues[ref]
                if occursin("Mu°",string(ref))   MuMu_[Symbol(string(k)*"_sp"*string(sp))] = MuMu[ref] end
            end
        else
#             println("_else 1")
            error("It seems a case that happen in TraitAffectHostsSymbiontsCells is not handeled by DuplicateTraitsValuesINIValuesForEachSp!, k=:$k")
        end
    end
    empty!(d)
    empty!(TraitAffectHostsSymbiontsCells)
    empty!(TraitsValues)
    empty!(MuMu)
    for k in keys(D)
        d[k]=D[k]
        TraitAffectHostsSymbiontsCells[k]=TraitAffectHostsSymbiontsCells_[k]
        TraitsValues[k]=TraitsValues_[k]
        if occursin("Mu°",string(k))  MuMu[k]=MuMu_[k] end
    end
end

mutable struct NamedMatrix
    M::Union{Matrix,Array{Any,1}}
    rownames::Union{Vector{Symbol},Symbol}
    colnames::Union{Vector{Symbol},Symbol}
    Nrow::Int
    Ncol::Int 
end

NamedMatrix(M::Matrix,rownames::Vector{Symbol},colnames::Vector{Symbol}) = 
    NamedMatrix(M::Matrix,rownames::Vector{Symbol},colnames::Vector{Symbol},
    size(M)[1] == length(rownames) ? size(M)[1] : error("size(Matrix)[1] !== length(rownames)")
   ,size(M)[2] == length(colnames) ? size(M)[2] : error("size(Matrix)[2] !== length(colnames)"))

NamedMatrix(M::Matrix,rownames::Vector{Any},colnames::Vector{Any}) = 
    NamedMatrix(M::Matrix
        ,convert(Vector{Symbol},rownames)
        ,convert(Vector{Symbol},colnames)
        ,size(M)[1] == length(rownames) ? size(M)[1] : error("size(Matrix)[1] !== length(rownames)")
        ,size(M)[2] == length(colnames) ? size(M)[2] : error("size(Matrix)[2] !== length(colnames)"))

NamedMatrix(M::Matrix,rownames::Union{Vector{Symbol},Symbol},colnames::Symbol) = 
    NamedMatrix(M::Matrix,rownames::Vector{Symbol},colnames::Vector{Symbol},
    size(M)[1] == length(rownames) ? size(M)[1] : error("size(Matrix)[1] !== length(rownames)")
   ,1 == length(colnames) ? 1 : error("size(Matrix)[2] !== length(colnames)"))

NamedMatrix(M::Matrix,rownames::Symbol,colnames::Union{Vector{Symbol},Symbol}) = 
    NamedMatrix(M::Matrix,rownames::Vector{Symbol},colnames::Vector{Symbol},
    1 == length(rownames) ? 1 : error("size(Matrix)[1] !== length(rownames)")
   ,size(M)[2] == length(colnames) ? size(M)[2] : error("size(Matrix)[2] !== length(colnames)"))

function length(::Number) 1 end
function length(::Symbol) 1 end

function getindex(M::NamedMatrix,R::Union{Vector{Symbol},Symbol},C::Union{Vector{Symbol},Symbol})
    @inbounds Rint = filter(i -> isinside(M.rownames[i],R), 1:M.Nrow)
    @inbounds Cint = filter(j -> isinside(M.colnames[j],C), 1:M.Ncol)
    return NamedMatrix(M.M[Rint,Cint],R,C,length(R),length(C))
end

function getindex(M::NamedMatrix,R::Union{Vector{Int},Int,UnitRange{Int}},C::Union{Vector{Int},Int,UnitRange{Int}})
    return NamedMatrix(M.M[R,C],M.rownames[R],M.colnames[C],length(R),length(C))
end

function getindex(x::NamedTuple, i::Integer)
   collect(x)[i]
end

function getindex(x::Base.RefValue, i::Integer)
   x[][i]
end





#                                                 AlleleINI = deepcopy(TraitsValuesINI)
function SetAlleleIniToTheClosestTraitsValuesIni!(d        ,Diff,TraitsValues,SpeciesCanMutateTo,SimParam)
    IsVec = NaN
    for k in keys(d)
#     println(k)
        if typeof(TraitsValues[k][1]) <: Vector
            ii = [findmin(map(i -> abs(d[k][l] - TraitsValues[k][i][l]), eachindex(TraitsValues[k]))) for l in 1:length(TraitsValues[k][1])]
            l = findmin(map(i -> abs(ii[i][2] - mean([ii[l][2] for l in eachindex(ii)])), 1:2))[2]
            ii = (ii[l][1], ii[l][2])
            if ((ii[1]/d[k][l])>1e2) error("For the parameter $k, Diff = $(ii[1])") end
            if (ii[1]>0) & ((ii[2]==1)|(ii[2]==length(TraitsValues[k]))) error("For the parameter $k, the initial value is out side the range of parameters") end
            push!(Diff,d[k][l]==0 ? ii[1] : ii[1]/d[k][l])
            d[k] = Int(ii[2])
        else
            ii = findmin(map(i -> abs(d[k] - TraitsValues[k][i]), eachindex(TraitsValues[k])))
            if ((ii[1]/d[k])>1e2) error("For the parameter $k, Diff = $(ii[1])") end
            if (ii[1]>0) & ((ii[2]==1)|(ii[2]==length(TraitsValues[k]))) error("For the parameter $k, the initial value is out side the range of parameters") end
            push!(Diff,d[k]==0 ? ii[1] : ii[1]/d[k])
            d[k] = Int(ii[2])
        end
    end
end

# Storage Type Conversion
# Recursive array to tuple conversion
function RecursiveArrayToTuple(a,i=0)
    A = deepcopy(a)
    a = Array{Any,1}()
    for aa in A push!(a,aa) end 
    for i in eachindex(a)
#         println(i)
        if typeof(a[i])<:Vector
            a[i] = RecursiveArrayToTuple(a[i])
        end
#         println(i/length(a))
    end
    a = tuple(a...)
    return a
end

# Make NamedTuple from dictionary
# Simple
namedtuple(d::Dict{Any,T}) where {T} =
    NamedTuple{(Symbol.(collect(keys(d)))...,)}((collect(values(d))...,))


# Recursive
function MakeNamedTupleRecursive(d,i=0)
    D=deepcopy(d)
    d=Dict{Any,Any}()
    for k in keys(D)
#         println(k," _")
        d[k] = D[k]
    end
#     println("Ok")
    for k in keys(d)
#         println(k)
        if typeof(d[k])<:Dict
#             println("if $k")
            d[k] = MakeNamedTupleRecursive(d[k],k)
        elseif typeof(d[k])<:Vector
#             println("else $k")
            d[k] = RecursiveArrayToTuple(d[k],k)
        end
    end
    d=NamedTuple{(Symbol.(collect(keys(d)))...,)}((collect(values(d))...,))
    return d
end


# TraitsValuesINI = UnfoldDirectory(temp,"!",(:Preaction,:DeadDecomp))

function UnfoldDirectory!(d::Dict,d_::Dict,BaseKey::Union{Symbol,String},sep::String,KeepFold::Array{Symbol})
    for k in keys(d)
#         println(k)
        if typeof(BaseKey) == Symbol
            kk_ = Symbol(BaseKey,"_",k)
        else
            kk_ = BaseKey*"_"*k
        end
        if (length(d[k])>1) & !any(KeepFold.==k)
            for i in 1:length(d[k])
                if (typeof(d[k]) <: Dict) | (typeof(d[k]) <: NamedTuple)
#                     println("\n")
#                     println("----------------",k)
                    UnfoldDirectory!(d[k],d_,kk_,sep,KeepFold)
                else
                    if typeof(k) == Symbol
                        kk = Symbol(kk_,sep,i)
                    else
                        kk = kk_*sep*string(i)
                    end
#                     println(d[k])
#                     println(kk)
                    d_[kk] = d[k][i]
                end
            end
        else
            d_[kk_] = d[k]
        end
    end
end

function UnfoldDirectory(d::Dict,sep,KeepFold::Array{Symbol})
    d_ = Dict{Any,Any}()
    for k in keys(d)
        if (length(d[k])>1) & !any(KeepFold.==k)
            for i in 1:length(d[k])
                if (typeof(d[k]) <: Dict) | (typeof(d[k]) <: NamedTuple)
#                     println("\n")
#                     println("----------------",k)
                    UnfoldDirectory!(d[k],d_,k,sep,KeepFold)
                else
                    if typeof(k) == Symbol
                        kk = Symbol(k,sep,i)
                    else
                        kk = k*sep*string(i)
                    end
#                     println(d[k])
#                     println(kk)
                    d_[kk] = d[k][i]
                end
            end
        else
            d_[k] = d[k]
        end
    end
    return(d_)
end



function PrintHierarchicalStruct(d ; SpaceAddedAtEachLevels=10,add = "")
    for k in keys(d)
        print("\n",add,k," => ")
        if typeof(d[k])<:Dict
            PrintHierarchicalStruct(d[k],SpaceAddedAtEachLevels=SpaceAddedAtEachLevels,add=add*" "^SpaceAddedAtEachLevels)
        else
            print(string(d[k]))
        end
    end
end

function PrintHierarchicalStruct(d              ,file::String,opened=false,FILE=nothing ; SpaceAddedAtEachLevels=10,add = "",heading="")
    function f(d,FILE)
#         save("/supr/supr.jld", "d",d)
        first = true
        for k in collect(keys(d))[sortperm(string.(collect(keys(d))))]
#         println(k)
            if typeof(k)<:Symbol
#                 println("if")
                if first
                    write(FILE, "\n"*add*" :"*string(k)*" => ")
                else
                    write(FILE, "\n"*add*",:"*string(k)*" => ")
                end
            elseif typeof(k)<:String
#                 println("elseif")
                if first
                    write(FILE, "\n"*add*" \""*string(k)*"\""*" => ")
                else
                    write(FILE, "\n"*add*",\""*string(k)*"\""*" => ")
                end
            else
#                 println("else")
                if first
                    write(FILE, "\n"*add*" "*string(k)*" => ")
                else
                    write(FILE, "\n"*add*","*string(k)*" => ")
                end
            end
            if typeof(d[k])<:Dict
#                 println("Dict")
                write(FILE, "Dict( ")
                PrintHierarchicalStruct(d[k],file,true,FILE;SpaceAddedAtEachLevels=SpaceAddedAtEachLevels,add=add*" "^SpaceAddedAtEachLevels)
                write(FILE, "\n"*add*")")
            elseif length(d[k])>1
#                println("Val")
                write(FILE, stringNothingInVec(d[k]))
            else
                write(FILE, string(d[k]))
            end
            first = false
        end
    end
    if !opened
        open(file, "w") do FILE
            write(FILE, heading)
            f(d,FILE)
            write(FILE, "\n)")
        end
    else
        f(d,FILE)
    end
end


function MakeUniformProportions(Nalleles,k)
    # Nalleles = n^k
    # n = Nalleles^(1/k)
    n = round(Nalleles^(1/k))
    U = seq(0,1,n-1)
    P = seq(0,1,n-1)
    while length(size(P))==1 || size(P)[end] < (k)
        add = flat([ [deepcopy(u) for i in 1:size(P)[1]] for u in U])
        P = cat([P for i in 1:n]...;dims=1)
        P = hcat(P,add)
        # println(size(P)[end])
    end
    P = map(i -> P[i,:] = P[i,:]./sum(P[i,:]), 1:size(P)[1])
    P = P[ filter(i -> !isnan(P[i][1]) & !isinf(P[i][1]), 1:size(P)[1]) ,:]
    Names = map(pp -> string(pp), P)
    Names_ = unique(Names)
    P = P[map(name -> filter(i -> name == Names[i]  , eachindex(Names))[1] ,Names_),:]
    
    return P
end    


function StandardizeAllelesSpecificityMatrix!(SpecificityMatrix::Matrix{Float})
    SD_Tot = std(SpecificityMatrix[1:end])
    while (var([mean(SpecificityMatrix[i,:]) for i in 1:NallelesSusceptibilities]) > 1e-30)
        M = [mean(SpecificityMatrix[i,:]) for i in 1:NallelesSusceptibilities]
        for i  in 1:NallelesSusceptibilities
            SpecificityMatrix[i,:] = SpecificityMatrix[i,:]./M[i]
        end
        M = [mean(SpecificityMatrix[:,j]) for j in 1:NallelesAggressiveness]
        for j in 1:NallelesAggressiveness
            SpecificityMatrix[:,j] = SpecificityMatrix[:,j]./M[j]
        end
        SpecificityMatrix[1:end] = SpecificityMatrix[1:end]./(std(SpecificityMatrix[1:end])/SD_Tot)
    end
end

function IndexesFromList(x::AbstractMultiScaleArray, Indexes::Union{Tuple,Vector})
    for i in Indexes
        x=x.nodes[i]
    end
    return(x)
end

function IndexesFromList(x::Union{Tuple,NamedTuple,Vector}, Indexes::Union{Tuple,Vector})
    for i in Indexes
        x=x[i]
    end
    return x
end

;

# SUPR
# using JLD
# d="gf"
# d = load("/supr/supr.jld")
# d=d["d"]

#########################################  
############### println(X::AbstractMultiScaleArray
function addSeparator!(firstElement,ToPrint,level)
    if firstElement
        ToPrint[level][end] = "+"*ToPrint[level][end]
#         println(ToPrint)
    end
    return false
end

function Take1stChar(s::Union{String,SubString{String}},N::Int64)     s[1:min(length(s),N)]     end

function ToPrintAbstractMultiScaleArray!(ToPrint,X::AbstractMultiScaleArray{Float},fields,LevelMax,NcharPerName,level = 1)
    if length(ToPrint) < level
        push!(ToPrint,Vector{String}())
    end
    firstElement = true
    if level<LevelMax
        for x in X.nodes
            if fieldname(typeof(x),1) == :nodes
                push!(ToPrint[level],  "|"*Take1stChar(split(string(typeof(x)),"{")[1],NcharPerName)  )
                firstElement = addSeparator!(firstElement,ToPrint,level)
                ToPrintAbstractMultiScaleArray!(ToPrint,x,fields,LevelMax,NcharPerName,level+1)
            else
                if fields == nothing
                    push!(ToPrint[level],  Take1stChar(split(string(typeof(x)),"{")[1],NcharPerName)  )
                    firstElement = addSeparator!(firstElement,ToPrint,level)
                else
                    if !(typeof(x) <: AbstractMultiScaleArrayLeaf{Float})
                        push!(ToPrint[level],  Take1stChar(split(string(typeof(x)),"{")[1],NcharPerName)  )
                        firstElement = addSeparator!(firstElement,ToPrint,level)
                    else
                        push!(ToPrint[level],  join([Take1stChar(string(field),NcharPerName)*": "*string(getfield(x,Symbol(field))) for field in collect(fields)],", ")  )
                        firstElement = addSeparator!(firstElement,ToPrint,level)
                    end
                end
            end
        end
    end
end

function print_human_readable(X::AbstractMultiScaleArray ; NcharPerName = 6, fields = nothing, LevelMax=Inf, NItemMaxPerLevels = Inf) # fields = nothing OR [field1,field2,...]
#     if typeof(X) <: AbstractMultiScaleArrayLeaf
#         println(X)
#     else
        ToPrint = Vector{Vector{String}}([[]]) # one vector per row, one string per element
        ToPrintAbstractMultiScaleArray!(ToPrint,X,fields,LevelMax,NcharPerName)
        toprint = map(x -> join(x,"; "),ToPrint)
    #     toprint[1] = join(split(toprint[1]," "),"; ")
        Done = false
        while !Done
            Done = true # temporary
            for level1 in length(toprint):-1:2
                toprint0 = split(toprint[level1-1],"|")
                toprint1 = split(toprint[level1  ],"+")
                toprint0_ , toprint1_ = "", ""
                for (x0, x1) in zip(toprint0,toprint1)
                    if length(x0) > length(x1)
                        x1=join([x1,repeat(" ",length(x0)-length(x1))],"")
                        Done = false
                    elseif length(x0) < length(x1)
                        x0=join([x0,repeat(" ",length(x1)-length(x0))],"")
                        Done = false
                    end
                    toprint0_ = join([toprint0_, x0],"|")
                    toprint1_ = join([toprint1_, x1],"+")
                end
                toprint0_ = toprint0_[2:end]
                toprint1_ = toprint1_[2:end]
                toprint[level1-1] = toprint0_
                toprint[level1  ] = toprint1_
            end
        end
        # handle "ula; +         |Popu" cases
        #              ^^^^^^^^^^^
        for level in length(toprint):-1:1
            while occursin("+ ", toprint[level])
                Move = findfirst("+ ", toprint[level])[1]
                To   = findnext(" |", toprint[level], Move[end])[1]
                toprint[level] = join([toprint[level][1:(Move-1)],  toprint[level][(Move+1):To] , toprint[level][Move]  , toprint[level][(To+1):end]])
                for level_ in 1:(level-1)
                    toprint[level_] = join([toprint[level_][1:(Move-1)],  toprint[level_][(To+1):(To+To-Move)],  toprint[level_][Move:To],  toprint[level_][(To+To-Move+1):end]])
                end
            end
        end
        println(join(toprint,"\n"))
#     end
end

