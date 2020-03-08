# MapOfMutationalPathSp = param



















function Make_AllelesIniBySp(sp::Int,cellType::Symbol,AlleleINI::NamedTuple,TraitAffect::typeof(AllelesParam.TraitAffectHostsSymbiontsCells),alleles=Dict())
    sp_ = "sp"*string(sp)
    for k in filter(k -> split(string(k),"_")[end]==sp_ ,collect(keys(AlleleINI)))
#         println(k)
#         if (typeof(AlleleINI[k])<:NamedTuple) & ((cellType !== :Symbiont) | (k !== :TraitByComp))
        if any(TraitAffect[k] .== cellType) & (!isnan(AlleleINI[k]))
#                 println("AlleleINI[k] = ",AlleleINI[k])
                alleles[k] = Int_(AlleleINI[k], typeof(maximum([maximum(AllelesParam.Nalleles[k]) for k in keys(AllelesParam.Nalleles[k])])))
#                 println("done")
        elseif !any(TraitAffect[k] .== cellType)
            nothing
        else
            error("In the function 'Make_AllelesIniBySp', a key k=$k is not well handeled")
        end
    end
#     println("alleles: ",alleles)
    return MakeNamedTupleRecursive(alleles)
end

allelesHostSymb = [Make_AllelesIniBySp(sp,[:HostIndOrCell,fill(:Symbiont,SimParam.NsymbiontSp)...][sp],  AllelesParam.AlleleINI,  AllelesParam.TraitAffectHostsSymbiontsCells) for sp in 1:SimParam.Nsp];


function GetMuEachTrait(alleles::NamedTuple, sp::typeof(SimParam.Nsp), cell::Val{true}, MuMu::typeof(AllelesParam.MuMu), TraitsValues, TraitAffectHostsSymbiontsCells)
    MuEachTrait = Dict()
    for k in keys(alleles)
#         println(k)
        ref = split(string(k),"_")[end-1]
        ref = split(ref,"!")[1]
        MuSymbol = Symbol(:Mu°,ref,:_sp,sp)
        if isinside(:Symbiont,TraitAffectHostsSymbiontsCells[k]) | isinside(:HostCells,TraitAffectHostsSymbiontsCells[k]) # cell::Val{true}
            if string(k)[1:3] == "Mu°"
                MuEachTrait[k] = MuMu[k]
            else
                MuEachTrait[k] = TraitsValues[MuSymbol][alleles[MuSymbol]]
            end
        end
    end
    return MakeNamedTupleRecursive(MuEachTrait);
end

function GetMuEachTrait(alleles::NamedTuple, sp::typeof(SimParam.Nsp), cell::Val{false}, MuMu::typeof(AllelesParam.MuMu), TraitsValues, TraitAffectHostsSymbiontsCells)
    MuEachTrait = Dict()
    for k in keys(alleles)
#         println(k)
        ref = split(string(k),"_")[end-1]
        ref = split(ref,"!")[1]
        MuSymbol = Symbol(:Mu°,ref,:_sp,sp)
        if isinside(:OnlyHostInd, TraitAffectHostsSymbiontsCells[k]) # cell::Val{false}
            if string(k)[1:3] == "Mu°"
                MuEachTrait[k] = MuMu[k]
            else
                MuEachTrait[k] = TraitsValues[MuSymbol][alleles[MuSymbol]]
            end
        end
    end
    return MakeNamedTupleRecursive(MuEachTrait);
end

# 
# x = GetMuEachTrait(allelesHostSymb[1], convert(typeof(SimParam.Nsp),1) , Val(false), AllelesParam.MuMu, TraitsValues, AllelesParam.TraitAffectHostsSymbiontsCells)

for sp in typeof(SimParam.Nsp)(1):SimParam.Nsp
    str = "function GetMuTot(x::typeof( GetMuEachTrait(allelesHostSymb[$sp], convert(typeof(SimParam.Nsp),$sp) , Val(true ) , AllelesParam.MuMu, TraitsValues, AllelesParam.TraitAffectHostsSymbiontsCells) ))
            return 1.0-prod(map(X -> 1-X, x))
        end"
    eval(Meta.parse(str))
end

for sp in typeof(SimParam.Nsp)(1):SimParam.Nsp
    str = "function GetMuTot(x::typeof( GetMuEachTrait(allelesHostSymb[$sp], convert(typeof(SimParam.Nsp),$sp) , Val(false) , AllelesParam.MuMu, TraitsValues, AllelesParam.TraitAffectHostsSymbiontsCells) ))
            return 1.0-prod(map(X -> 1-X, x))
        end"
    eval(Meta.parse(str))
end


function KeysToPathAlleles(Model::Union{typeof(allelesHostSymb[1]),typeof(allelesHostSymb[2])},Keys::Union{Vector{Union{Symbol,Integer}},NTuple{N,Union{Symbol,Integer}} where N})
    m=Model
    Path=Vector{Int8}()
    for k in Keys
        push!(Path,filter(i -> keys(m)[i]==k,Int8(1):length(m))[1])
        m = m[k]
    end
    return Path
end


# Model,key,PrevVal,NewVal,pathAll = 
# X.alleles, :EnterMminusCell,  0x01,  0x04 ,  MapMutatingTraits.PathAll[7]

# Mutate_Alleles(X.alleles, :prodR,  0x01,  0x04 ,  MapMutatingTraits.PathAll[53]) # sp > 1
# Mutate_Alleles(X.alleles, :prodR,  0x01,  0x04 ,  MapMutatingTraits.PathAll[120]) # sp = 1

# # Mutate_Alleles(X.alleles, :GetOutMplusCell,  0x01,  0x04 ,  MapMutatingTraits.PathAll[20]) #  sp = 1 & sp > 1
# Mutate_Alleles(X.alleles, :GetOutMplusHost,  0x01,  0x04 ,  MapMutatingTraits.PathAll[7]) #  sp = 1 & sp > 1



str="""
function Mutate_Alleles(Model::Union{typeof(allelesHostSymb[1]),typeof(allelesHostSymb[2])},key::Symbol,PrevVal::Integer,NewVal::Integer,pathAll::Union{Vector{Union{Symbol,Integer}},NTuple{N,Union{Symbol,Integer}} where N})
    if length(pathAll) == 1
        PrevVal = split(   split(string([PrevVal]),"[")[end]   ,"]" )[1]
        NewVal  = split(   split(string([NewVal]) ,"[")[end]   ,"]" )[1]
        #
        x1 = string(key)*" = "*PrevVal
        x2 = string(key)*" = "*NewVal
        return( eval(Meta.parse(    replace(string(Model), x1=>x2)   ))   )
    else
        posiKey = findlast(K -> K==key,pathAll)[1]
        PathToKey = KeysToPathAlleles(Model,pathAll[1:posiKey])
        #
        if posiKey == length(pathAll)
            PrevVal = split(   split(string([PrevVal]),"[")[end]   ,"]" )[1]
            NewVal  = split(   split(string([NewVal]) ,"[")[end]   ,"]" )[1]
            x1 = string(key)*" = "*PrevVal
            x2 = string(key)*" = "*NewVal
        else
            x1 = Model
            kReached = false
            for p in pathAll
                if !kReached
                    x1 = x1[p]
                end
                if kReached
                    x2 = collect(x2)
                    if (x2[p]!==PrevVal) error("In Mutate_Alleles: x2[p]!==PrevVal. pathAll = "*string(pathAll)) end
                    x2[p] = NewVal
                    x2 = tuple(x2...)
                end
                if (p == key) kReached=true ; x2 = deepcopy(x1) end
            end
            x2 = string(key)*" = "*string(x2)
            x1 = string(key)*" = "*string(x1)
        end
        if posiKey == 1
            return( eval(Meta.parse(    replace(string(Model), x1=>x2)   ))   )
        else
            Start = Vector{String}([])
            End = Vector{String}([])
            K = Vector{Symbol}([])
            m = Model
            for p in PathToKey[1:(end-1)]
                push!(K,keys(m)[p])
                push!(Start,join([typeof(m[i])<:Symbol ? string(keys(m)[i])*" = :"*string(m[i]) : string(keys(m)[i])*" = "*string(m[i]) for i in 1:(p-1)],", "))
                push!(End  ,join([typeof(m[i])<:Symbol ? string(keys(m)[i])*" = :"*string(m[i]) : string(keys(m)[i])*" = "*string(m[i]) for i in (p+1):length(m)],", "))
                m = m[Int(p)]
            end
            m = replace(string(m), x1=>x2)
            Reordered = Vector{String}(["("])
            for (s,k) in zip(Start,K)
                if ((s!=="") & (s!=="(")) push!(Reordered,s,", ") end
                push!(Reordered,string(k)," = (")
            end
            push!(Reordered,m[2:(end-1)])
            for e in End[end:-1:1]
                push!(Reordered,"), ")
                push!(Reordered,e)
            end
            push!(Reordered,")")
            $PRINT println(join(Reordered))
            return( eval(Meta.parse(    join(Reordered)   ))   )
    #         Alleles = join(Reordered)
        end
    end
end"""
eval(Meta.parse(str))

# SUPR
# Model,key,PrevVal,NewVal,PathToKey = 
# X.alleles, k, PrevVal,NewVal ,KeysToPathAlleles(X.alleles,path) 

########################

### Where a given cell genotype is mutating when a mutation happen
# str = """#                      alleles
# function MakeMapMutatingTraits!(x::NamedTuple,cell::Bool, out::Vector{Vector{Any}}, posi::Vector{Union{Int8,Symbol}}, Mu, K::Union{Int,Symbol}, Floor::Int)
#         for k in keys(x)
# #         println(k)
#         if length(posi) == Floor
#             posi[Floor] = k
#         elseif length(posi) > Floor
#             posi = posi[1:Floor]
#             posi[Floor] = k
#         else
#             push!(posi,k)
#             if length(posi) !== Floor
#                 error("length(posi) !== Floor. There is an error in the code of the function 'MakeMapMutatingTraits'.")
#             end
#         end
#         if isa(x[k] , Integer)
#             OnlyHostInd = isinside(:OnlyHostInd,$(AllelesParam.TraitAffectHostsSymbiontsCells)[k])
#             if !isnothing(Mu[k]) & (cell ? !OnlyHostInd : OnlyHostInd)
#                 push!(out[1],$(AllelesParam.TraitsValues.Mu)[Mu[k]])
#                 push!(out[2],deepcopy(posi))
#                 push!(out[3],deepcopy(posi))
#             end
#         else
# #             println(k,"   cell=", cell,"   out=", out,"   posi=", posi,"   k=",  k,"   Floor=", Floor+1)
#             MakeMapMutatingTraits!(x[k], cell, out, posi, Mu, k, Floor+1)
#         end
#     end
# end"""
# eval(Meta.parse(str))
# 
# str = """#                      alleles
# function MakeMapMutatingTraits!(x::Tuple,cell::Bool, out::Vector{Vector{Any}}, posi::Vector{Union{Int8,Symbol}}, Mu, K::Union{Int,Symbol}, Floor::Int)
#         for k in 1:length(x)
# #         println(k)
#         if length(posi) == Floor
#             posi[Floor] = Int8(k)
#         elseif length(posi) > Floor
#             posi = posi[1:Floor]
#             posi[Floor] = Int8(k)
#         else
#             push!(posi,Int8(k))
#             if length(posi) !== Floor
#                 error("length(posi) !== Floor. There is an error in the code of the function 'MakeMapMutatingTraits'.")
#             end
#         end
#         if isa(x[k] , Integer)
#             OnlyHostInd = isinside(:OnlyHostInd,$(AllelesParam.TraitAffectHostsSymbiontsCells)[k])
#             if !isnothing(Mu[k]) & (cell ? !OnlyHostInd : OnlyHostInd)
#                 push!(out[1],$(AllelesParam.TraitsValues.Mu)[Mu[K]])
#                 push!(out[2],deepcopy(posi))
#                 push!(out[3],deepcopy(posi[1:(end-1)]))
#             end
#         else
# #             println(k,"   cell=", cell,"   out=", out,"   posi=", posi,"   k=",  k,"   Floor=", Floor+1)
#             MakeMapMutatingTraits!(x[k], cell, out, posi, Mu, k, Floor+1)
#         end
#     end
# end"""
# eval(Meta.parse(str))

# str = """#                        alleles
# function MakeMapMutatingTraits(x::NamedTuple,cell::Bool)
#     out = Vector{Vector{Any}}[[],[],[]]
#     posi = Vector{Union{Int8,Symbol}}(undef,0)
#     Floor = 1
#     Mu=x.Mu
#     for k in keys(x)
# #     #     println(k)
#         if length(posi) == Floor
#             posi[Floor] = k
#         elseif length(posi) > Floor
#             posi = posi[1:Floor]
#             posi[Floor] = k
#         else
#             push!(posi,k)
#             if length(posi) !== Floor
#                 error("length(posi) !== Floor. There is an error in the code of the function 'MakeMapMutatingTraits'.")
#             end
#         end
#         if isa(x[k] , Integer)
#             OnlyHostInd = isinside(:OnlyHostInd,$(AllelesParam.TraitAffectHostsSymbiontsCells)[k])
#             if !isnothing(Mu[k]) & (cell ? !OnlyHostInd : OnlyHostInd)
#                 push!(out[1],$(AllelesParam.TraitsValues.Mu)[Mu[k]])
#                 push!(out[2],deepcopy(posi))
#                 push!(out[3],deepcopy(posi))
#             end
#         else
# #             println(k,"   cell=", cell,"   out=", out,"   posi=", posi,"   k=",  k,"   Floor=", Floor+1)
#             MakeMapMutatingTraits!(x[k], cell, out, posi, Mu, k, Floor+1)
#         end
#     end
#     out = ( Prob = [x/sum(out[1]) for x in out[1]]
#                ,PathAll = RecursiveArrayToTuple(out[2])
#                ,PathToTrait = RecursiveArrayToTuple(out[3])
#                ,traits = tuple([filter(xx -> isa(xx, Symbol) ,x)[end] for x in out[2]]...)
#                )
#     return(out)
# end"""
# eval(Meta.parse(str))

# SUPR
# MakeMapMutatingTraits(allelesHostSymb[1], true)


# PathToKey = KeysToPathAlleles(allelesHostSymb[1], [ :TraitByComp  , :FemaleGerma  ,:OffspringComp])
# Model,key,PrevVal,NewVal = allelesHostSymb[1], :SizeInocula,126,500
# allelesHost2 = Mutate_Alleles(allelesHostSymb[1], :SizeInocula,126,500,PathToKey)

## Estimate the distance between pairs of genotypes START
DistMaxHost = Int_(sum([any(AllelesParam.TraitAffectHostsSymbiontsCells[k].==:HostIndOrCell) ? length(AllelesParam.TraitsValues) : 0 for k in keys(AllelesParam.TraitsValues)]))
DistMaxSymb = Int_(sum([any(AllelesParam.TraitAffectHostsSymbiontsCells[k].==:Symbiont     ) ? length(AllelesParam.TraitsValues) : 0 for k in keys(AllelesParam.TraitsValues)]))
DistBetweenGenotypes = ntuple(sp -> Vector{Tuple}([(0.0,)]) , SimParam.Nsp)

# # SUPR
# function IterativeComparaison!(X::NamedTuple,Y::NamedTuple,Tot::Vector{typeof(DistMaxHost)})
#     for (x,y,k) in zip(X,Y,keys(X))
#         if (typeof(x)<:Number)
#             Tot[1] += abs(x - y)
#         elseif (typeof(x)<:Tuple) | (typeof(x)<:Array) # | all(map(xx -> typeof(xx)<:Number,x))
#             Tot[1] += sum(map(i -> abs(x[i] - y[i]), eachindex(x)))
#         elseif (typeof(x)<:NamedTuple) | (typeof(x)<:Dict)
#             IterativeComparaison!(x,y,Tot)
#         elseif typeof(x)<:Symbol
#             Tot[1] += x == y ? 0 : 1
#         else
#             error("A case happen that is not handeled by IterativeComparaison. This is with x, y, k=$x, $y, $k")
#         end
#     end
# end

ListOfGenotypes = tuple([Array{NamedTuple}(undef,1) for sp in 1:SimParam.Nsp]...)
for sp in 1:SimParam.Nsp
    ListOfGenotypes[sp][1] = allelesHostSymb[sp]
end

ListOfPhenotypes = tuple([Array{NamedTuple}(undef,1) for sp in 1:SimParam.Nsp]...)
for sp in typeof(SimParam.Nsp)(1):SimParam.Nsp
    ListOfPhenotypes[sp][1] = MakePhenotypesFromAlleles(allelesHostSymb[sp],sp)
end

str = "function EstimateDistBetweenGenotypes(a1::typeof(allelesHostSymb[1]),a2::typeof(allelesHostSymb[1]),sp::Val{$(typeof(SimParam.Nsp))(1)})
    # Here there is an approximtion, because mutations of traits that affect only host level phenotypes are only simulated when it reproduces. But this save a lot of time.
    return sum(map((x,y) -> abs(x-y), a1,a2)) / $DistMaxHost
end"
eval(Meta.parse(str))

for sp in typeof(SimParam.Nsp)(2):SimParam.Nsp
    str = "function EstimateDistBetweenGenotypes(a1::typeof(allelesHostSymb[$sp]),a2::typeof(allelesHostSymb[$sp]),sp::Val{$(typeof(SimParam.Nsp))($sp)})
        return sum(map((x,y) -> abs(x-y), a1,a2))/$DistMaxSymb
    end"
    eval(Meta.parse(str))
end

# EstimateDistBetweenGenotypes(allelesHostSymb[2],allelesHostSymb[2],Val(typeof(SimParam.Nsp)(2)))
# EstimateDistBetweenGenotypes(allelesHostSymb[1],allelesHostSymb[1],Val(typeof(SimParam.Nsp)(1)))

## Estimate the distance between pairs of genotypes END

# Estimate Map of Mutational Path
# function Make_AllelesIniBySp_TEMP(sp::Int,cellType::Symbol,AlleleINI::NamedTuple,TraitAffect::typeof(AllelesParam.TraitAffectHostsSymbiontsCells),alleles=Dict(),out=true)
#     for k in keys(AlleleINI)
# #         if (typeof(AlleleINI[k])<:NamedTuple) & ((cellType !== :Symbiont) | (k !== :TraitByComp))
#         if (typeof(AlleleINI[k])<:NamedTuple) # & (k !== :TraitByComp)
#             alleles[k] = Dict()
#             Make_AllelesIniBySp_TEMP(sp,cellType,AlleleINI[k],TraitAffect,alleles[k],false)
#             if alleles[k] == Dict()
#                 delete!(alleles,k)
#             end
#         elseif any(TraitAffect[k] .== cellType) & (!isnothing(AlleleINI[k][sp]))
#                 alleles[k] = Int_(AlleleINI[k][sp],Int)
#         elseif !any(TraitAffect[k] .== cellType)
#             nothing
#         else
#             error("In the function 'Make_AllelesIniBySp', a key k=$k is not well handeled")
#         end
#     end
#     if out return MakeNamedTupleRecursive(alleles) end
# end
# 
# allelesHostSymb_TEMP = [Make_AllelesIniBySp_TEMP(sp,[:HostIndOrCell,fill(:Symbiont,SimParam.NsymbiontSp)...][sp],  AllelesParam.AlleleINI,  AllelesParam.TraitAffectHostsSymbiontsCells) for sp in 1:SimParam.Nsp]


# # NamedTuple(for each mutating trait Array(for each genoype.ID two Array(mutating by deacreasing.increasing trait value; position 1, 2, 3, ... = ID of the genotypes)) ))
N=maximum(flat(AllelesParam.Nalleles))
T = Int
temp   = Array{T}(undef,N).*0
temp2  = Array{typeof(tuple(deepcopy(temp), deepcopy(temp)))}(undef,1)
temp2[1] = tuple(deepcopy(temp), deepcopy(temp))
MapOfMutationalPath = Dict()
for sp in 1:SimParam.Nsp
    for k in keys(allelesHostSymb[sp])
        MapOfMutationalPath[k] = deepcopy(temp2)
    end
end
MapOfMutationalPath = NamedTuple{(Symbol.(collect(keys(MapOfMutationalPath)))...,)}((collect(values(MapOfMutationalPath))...,))

struct Parameters{B<:Float} <: AbstractMultiScaleArrayLeaf{B} #  AbstractMultiScaleArrayLeaf    *** Leaf ***
    values::Vector{B}
    Iam::Symbol
    SimParam::typeof(SimParam)
    AllelesParam::typeof(AllelesParam)
    DistBetweenGenotypes::NTuple{Int64(SimParam.Nsp),Array{Tuple}}
     ListOfGenotypes::Tuple{[Array{NamedTuple} for sp in 1:SimParam.Nsp]...}
    ListOfPhenotypes::Tuple{[Array{NamedTuple} for sp in 1:SimParam.Nsp]...}
    NeighbourhoodActual::Array{Float} # internal storage used by some functions
    NeighbourhoodEstim::Array{Float} # internal storage used by some functions
    Agressions::Array{Float} # internal storage used by some functions
    Eabsorb::Ref{Float} # internal storage used by some functions
    UsableE::Ref{Float} # # internal storage used by some functions
    Prepro::Ref{Float} # # internal storage used by some functions
    NMu::Ref{Int} # # internal storage used by some functions
    NMu_expect::Ref{Float} # # internal storage used by some functions
    MapOfMutationalPath::typeof(MapOfMutationalPath)
    # 
    # Mutations
    IDnew::Ref{Int}
    nMu::Ref{Int}
    CurPosi::Ref{Int}
    NewPosi::Ref{Int}
    NewVal::Ref{Integer}
    CellGenotMuEnd::Ref{Union{Int,Nothing}}
end

#                      |  values | Iam  | SimParam| AllelesParam| DistBetweenGenotypes| ListOfGenotypes| ListOfPhenotypes|       NeighbourhoodActual     |NeighbourhoodEstim:sp,kin,sp,sp, SP|               Agressions                |Eabsorb |UsableE | Prepro |     NMu    |  NMu_expect  |MapOfMutationalPath|    IDnew   |      NMu    |   CurPosi  |   NewPosi  |    NewVal  |        CellGenotMuEnd      |
parameters = Parameters([0.0]    ,:Param, SimParam, AllelesParam, DistBetweenGenotypes, ListOfGenotypes, ListOfPhenotypes, [0.0 for sp in 1:SimParam.Nsp], [0.0 for sp in 1:(SimParam.Nsp+1)], [0.0 for sp in 1:NallelesAggressiveness],Ref(0.0),Ref(0.0),Ref(0.0),Ref{Int}(0),Ref{Float}(0),MapOfMutationalPath,Ref(Int(0)),Ref(Int(0)),Ref(Int(0)),Ref(Int(0)),Ref{Integer}(0),Ref{Union{Int,Nothing}}(0))

struct CellHostGenot{B<:Float} <: CellGenot{B} #  AbstractMultiScaleArrayLeaf    *** Leaf ***
    values::Vector{B} # N: Number of individuals of that genotype and E, Amount of Energy cells contain, Total amount of reproduction ; to be modifiable, "values" must be contained in a vector/array
    Iam::Symbol
    ID::Int # Unique identifing number of this cell genotype in Host. This is the position of this genotype in ListOfGenotypes and DistBetweenGenotypes
    sp::typeof(SimParam.Nsp)
    alleles::NamedTuple
    pheno::NamedTuple
    MuEachTraitCell::NamedTuple
    MuEachTraitHost::NamedTuple
    MuTotCell::Float
    MuTotHostOnly::Float # 1-P(no mutations at all)
#     Removed, this would need to be updated each time a new cell genotype appears
#     MapMutatingTraitsCell::typeof(MakeMapMutatingTraits(allelesHostSymb[1], true)) # Where is this particular cell genotype is mutating when a mutation happen
#     MapMutatingTraitsHostOnly::typeof(MakeMapMutatingTraits(allelesHostSymb[1], false)) # Where is this particular cell genotype is mutating when a mutation happen
end

struct CellSymbGenot{B<:Float} <: CellGenot{B} #  AbstractMultiScaleArrayLeaf    *** Leaf ***
    values::Vector{B} # N: Number of individuals of that genotype ; to be modifiable, "values" must be contained in a vector/array
    Iam::Symbol
    ID::Int # Unique identifing number of this cell genotype in each species. This is the position of this genotype in ListOfGenotypes and DistBetweenGenotypes
    sp::typeof(SimParam.Nsp)
    alleles::NamedTuple
    pheno::NamedTuple
    MuEachTraitCell::NamedTuple
    MuTotCell::Float # 1-P(no mutations at all)
#     Removed, this would need to be updated each time a new cell genotype appears
#      MapMutatingTraitsCell::typeof(MakeMapMutatingTraits(allelesHostSymb[2], true)) # Where is this particular cell genotype mutating when a mutation happen
end

struct SubEnvir{T<:AbstractMultiScaleArray,B<:Float}<:AbstractMultiScaleArray{B}
    nodes::Vector{T}
    values::Vector{B} # [Amount of ressources]
    end_idxs::Vector{Int}
end

struct ParamSubEnvir{B<:Float} <: AbstractMultiScaleArrayLeaf{B} #  AbstractMultiScaleArrayLeaf    *** Leaf ***
    values::Vector{B} # Volume
    Iam::Symbol
    Flow::Float
    ConcArriving::NTuple{Int64(SimParam.Nresources),Float}
    CurrentAgress::Array{Float}
    IndexNodesBySp::Tuple{[NamedTuple{(:Posi, :ID),Tuple{Array{Int,1},Array{Int,1}}} for sp in 1:SimParam.Nsp]...}
    PosiEnvir::Vector{Int64}
end

struct MainEnvir{T<:AbstractMultiScaleArray,B<:Float}<:AbstractMultiScaleArray{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
end

struct Container{T<:AbstractMultiScaleArray,B<:Float}<:AbstractMultiScaleArrayHead{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
end

SubEnvir(X::AbstractMultiScaleArrayLeaf{Float},n::Integer) = construct(SubEnvir,[deepcopy(X) for i in 1:n])
SubEnvir(X::AbstractMultiScaleArrayLeaf{Float}) = SubEnvir(X,1)
SubEnvir(X::Vector) = construct(SubEnvir,[deepcopy(x) for x in X])

SubEnvir(X::AbstractMultiScaleArrayLeaf{Float},n::Integer,val::Vector{Float}) = construct(SubEnvir,[deepcopy(X) for i in 1:n],val)
SubEnvir(X::AbstractMultiScaleArrayLeaf{Float},val::Vector{Float}) = SubEnvir(X,1,val)
SubEnvir(X::Vector,val::Vector{Float}) = construct(SubEnvir,[deepcopy(x) for x in X],val)

MainEnvir(X::AbstractMultiScaleArray{Float},n::Integer) = construct(MainEnvir,[deepcopy(X) for i in 1:n])
MainEnvir(X::AbstractMultiScaleArray{Float}) = MainEnvir(X,1)
MainEnvir(X::Vector) = construct(MainEnvir,[deepcopy(x) for x in X])

MainEnvir(X::SubEnvir,n::Integer,val::Vector{Float}) = construct(MainEnvir,[deepcopy(X) for i in 1:n],val)
MainEnvir(X::SubEnvir,val::Vector{Float}) = MainEnvir(X,1,val)
MainEnvir(X::Vector,val::Vector{Float}) = construct(MainEnvir,[deepcopy(x) for x in X],val)

Container(X::MainEnvir,n::Integer) = construct(Container,[parameters,[deepcopy(X) for i in 1:n]...])
Container(X::MainEnvir) = Container(X,1)
Container(X::Vector) = construct(Container,[parameters,[deepcopy(x) for x in X]...])

Container(X::MainEnvir,n::Integer,val::Vector{Float}) = construct(Container,[parameters,[deepcopy(X) for i in 1:n]...],val)
Container(X::MainEnvir,val::Vector{Float}) = Container(X,1,val)
Container(X::Vector,val::Vector{Float}) = construct(Container,[parameters,[deepcopy(x) for x in X]...],val)

# sp = 2
# sp_ = 1
# i = 2 ; posi = param.IndexNodesBySp[sp].Posi[i]


# Kill genotype N <= 0.0     a sub environement
str = """ 
function RemoveGenotypesN_0!(XX::AbstractMultiScaleArray, Head::Container{AbstractMultiScaleArray{Float64},Float64}, dHead::Container{AbstractMultiScaleArray{Float64},Float64})
    $PRINT println("Removing N=0 genotypes, typeof(XX) = ",typeof(XX))
    $PRINT println("RemoveGenotypesN_0    ",1)
    param = XX.nodes[1]
    if isa(param,AbstractMultiScaleArrayLeaf) && (param.Iam == :ParamSubEnvir) # | :Host
        for sp in 1:$(SimParam.Nsp)
            $PRINT println("RemoveGenotypesN_0    ",2,", sp = ",sp)
            for i in length(param.IndexNodesBySp[sp].Posi):-1:1
                posi = param.IndexNodesBySp[sp].Posi[i]
                $PRINT println("RemoveGenotypesN_0    ",2.333)
                $PRINT println("posi = ",posi, ",  length(XX.nodes) = ", length(XX.nodes))
                X = XX.nodes[posi]
                $PRINT println("RemoveGenotypesN_0    typeof(X)=",typeof(X),"  values[1] =", X.values[1],"   ",2.666)
                if isa(X,AbstractMultiScaleArrayLeaf)
                    $PRINT println("RemoveGenotypesN_0    ",2.75)
                    if (typeof(X)<:CellGenot) && (X.values[1] <= 0.0)
                        $PRINT println("RemoveGenotypesN_0    ",3, "posi = ", posi)
                        $PRINT println("N of cell genotypes in the SuvEnvir.s") ; [println([xx.values[1] for xx in x.nodes]) for x in Head.nodes[2].nodes]
                        $PRINT println("Removing N=0 genotype",IndexesFromList(Head,[param.PosiEnvir..., posi]).values)
                        $PRINT Len1 = length(XX.nodes)
                        remove_node!( Head, [param.PosiEnvir..., posi]...)
                        remove_node!(dHead, [param.PosiEnvir..., posi]...)
                        $PRINT println("N of cell genotypes in the SuvEnvir.s") ; [println([xx.values[1] for xx in x.nodes]) for x in Head.nodes[2].nodes]
                        $PRINT Len2 = length(XX.nodes)
                        $PRINT println("Len1 = ",Len1, ", Len2 = ",Len2)
                        $PRINT if Len1 == Len2
                                   $PRINT error("Len1 == Len2")
                        $PRINT end
                        # Update Posi
                        $PRINT println("RemoveGenotypesN_0    ",3.1, "posi = ", posi)
                        deleteat!(param.IndexNodesBySp[sp].Posi,i)
                        deleteat!(param.IndexNodesBySp[sp].ID  ,i)
                        $PRINT println("RemoveGenotypesN_0    ",3.2, "posi = ", posi)
                        for sp_ in 1:$(SimParam.Nsp)
                            $PRINT println("RemoveGenotypesN_0    ",3.3, "sp_ = ", sp_)
                            for ii in eachindex(param.IndexNodesBySp[sp_].Posi)
                                $PRINT println("RemoveGenotypesN_0    ",3.4, "ii = ", ii)
                                if param.IndexNodesBySp[sp_].Posi[ii] > posi
                                    $PRINT println("RemoveGenotypesN_0    ",3.5)
                                    param.IndexNodesBySp[sp_].Posi[ii] -= 1
                                    $PRINT println("RemoveGenotypesN_0    ",3.6)
                                end
                            end
                        end
                        $PRINT Max = maximum(flat([param.IndexNodesBySp[sspp].Posi for sspp in 1:$(SimParam.Nsp)]))
                        $PRINT if (Max !== length(XX.nodes))
                            $PRINT println("[param.IndexNodesBySp[sspp].Posi for sspp in 1:$(SimParam.Nsp)] = ",[param.IndexNodesBySp[sspp].Posi for sspp in 1:$(SimParam.Nsp)])
                            $PRINT println("param.IndexNodesBySp = ", param.IndexNodesBySp)
                            $PRINT println("length(XX.nodes) = ",length(XX.nodes), ", Max = ",Max, ", posi = ",posi)
                            $PRINT println("typeof.(XX.nodes) = ", typeof.(XX.nodes))
                            $PRINT error("Max!==length(XX.nodes)")
                        $PRINT end
                    end
                else
                    $PRINT println("RemoveGenotypesN_0    ",4, "posi = ", posi)
                    RemoveGenotypesN_0!(X, Head::AbstractMultiScaleArrayHead, dHead::AbstractMultiScaleArrayHead)
                end
            end
        end
    else
        $PRINT println("  - - length(XX.nodes) = ",length(XX.nodes), " typeof(XX) = ",typeof(XX))
        $PRINT ijii = 0
        for X in XX.nodes
            $PRINT ijii += 1
            $PRINT println(ijii)
            if !isa(X,AbstractMultiScaleArrayLeaf)
                $PRINT println("RemoveGenotypesN_0   ",5, ", typeof(X) = ", typeof(X))
                RemoveGenotypesN_0!(X, Head::AbstractMultiScaleArrayHead, dHead::AbstractMultiScaleArrayHead)
            end
        end
    end
end"""
eval(Meta.parse(str))

RemoveGenotypesN_0!(Head,dHead) = RemoveGenotypesN_0!(Head,Head,dHead)



function ChekGenotypesN_0(X::AbstractMultiScaleArray,str="")
    for x in X.nodes
        if isa(x,AbstractMultiScaleArrayLeaf)
            if isa(x,CellGenot) && x.values[1]<=0
                error("Cell genotype with N="*string(x.values[1])*"\n"*str)
            end
        else
            ChekGenotypesN_0(x,str)
        end
    end
end
    



#######################
# Mutate
# Map_ID_AllGenot = deepcopy(Map_ID_AllGenot_)

N=maximum([NallelesAggressiveness, NallelesSusceptibilities, Nalleles])
T = Int
Zero = split(   split(string([T(0)]),"[")[2]  ,"]")[1]
temp = Array{Tuple{Array{T},Array{T}}}(undef,1)
temp2  = Array{T}(undef,N).*0
temp[1] = tuple(deepcopy(temp2), deepcopy(temp2))


# function FindNumberOfTraitsDiffering!(X::NamedTuple,Y::NamedTuple,Tot::Ref{Union{Integer,Symbol}},pathALL::NTuple{N,NTuple{N,Union{Symbol,Integer}} where N} where N, cell::Val{true})
#     for pathAll in pathALL
#         
#     end
#     for (x,y,k) in zip(X,Y,keys(X))
#         if (typeof(x)<:Number)
#             if x!==y Tot[] += abs(x - y)
#         elseif (typeof(x)<:Tuple) | (typeof(x)<:Array) # | all(map(xx -> typeof(xx)<:Number,x))
#             Tot[1] += sum(map(i -> abs(x[i] - y[i]), eachindex(x)))
#         elseif (typeof(x)<:NamedTuple) | (typeof(x)<:Dict)
#             IterativeComparaison!(x,y,Tot)
#         elseif typeof(x)<:Symbol
#             Tot[1] += x == y ? 0 : 1
#         else
#             error("A case happen that is not handeled by FindNumberOfTraitsDiffering. This is with x, y, k=$x, $y, $k")
#         end
#     end
# end
# 
# 
# FOR EACH PAIR OF GENOTYPES
# for each trait => dist
# if they differ only of one trait then ...


# MapOfMutationalPathSp = param
# For each genotype that differ from X.ID only for 'trait' (and that are thus in Map_ID_AllGenot), add their position in Map_ID_AllGenot[IDnew[]] and add the position of IDnew[] in other genotypes
str="""
function Update_MapOfMutationalPathSp!(MapOfMutationalPathSp::NamedTuple, ListOfGenotypesSp::Vector{NamedTuple}
                                      ,IDold::Int, IDnew::Ref{Int}, NewAlleles::NamedTuple, pathALL::NTuple{N,NTuple{N,Union{Symbol,Integer}} where N} where N, cell::Val{true})
    PathToDifferingTrait = Ref{Union{Integer,NTuple{N,Union{NTuple,Symbol,Integer}} where N }}(0)
    Break = Ref{Bool}(false)
    Diff = Ref{Integer}(0)
    dist = Ref{Integer}(0)
    dir = Ref{Integer}(0)
    $PRINT println("Update_MapOfMutationalPathSp 1")
    for (ID,genot) in enumerate(ListOfGenotypesSp[1:(end-1)]) # ListOfGenotypesSp[end] is new ID
        $PRINT println("Update_MapOfMutationalPathSp 2, ID = ",ID)
        PathToDifferingTrait[] = 0 # reset PathToDifferingTrait
        Break[] = false # reset Break
        # Find if ID and IDnew only differ for one trait, and if yes, for which one
        for pathAll in pathALL
            $PRINT println("Update_MapOfMutationalPathSp 2.333, genot = ",genot, ", NewAlleles = ",NewAlleles,", pathAll = ",pathAll)
            genot_ = IndexesFromList(genot,pathAll)
            $PRINT println("Update_MapOfMutationalPathSp 2.666")
            NewAlleles_ = IndexesFromList(NewAlleles,pathAll)
            $PRINT println("Update_MapOfMutationalPathSp 3, genot_ = ",genot_, ", NewAlleles_ = ",NewAlleles_)
            if NewAlleles_ !== genot_
                if PathToDifferingTrait[] == 0
                    $PRINT println("Update_MapOfMutationalPathSp 3.333")
                    PathToDifferingTrait[] = pathAll
                    $PRINT println("Update_MapOfMutationalPathSp 3.5")
                    Diff[] = NewAlleles_ - genot_
                else # the two genotypes differ for more than one trait and cannot be linked by only one mutational event.
                    Break[] = true
                    break
                end
            end
        end
        if !Break[]
            $PRINT if (PathToDifferingTrait[] == 0) error("There is an error in the code of 'Mutate!' or of 'Update_MapOfMutationalPathSp!':\nwe are creating a new genotype that alredy exist!") end
            $PRINT println("Update_MapOfMutationalPathSp 3.666")
            if (sign(Diff[]) == 1) dir[]=2 else dir[]=1 end
            dist[] = abs(Diff[])
            $PRINT println("Update_MapOfMutationalPathSp 4, dist = ",dist[]," dir = ",dir[])
            Map_ID_AllGenot = MapOfMutationalPathSp
            for k in PathToDifferingTrait[]
                $PRINT println("Update_MapOfMutationalPathSp 4.2, k = ",k)
                Map_ID_AllGenot = Map_ID_AllGenot[k]
                $PRINT println("Map_ID_AllGenot = ",Map_ID_AllGenot)
            end
            $PRINT println("PathToDifferingTrait[] = ",PathToDifferingTrait[])
            $PRINT println("length(Map_ID_AllGenot) = ",length(Map_ID_AllGenot), ", ID = ",ID,", IDnew = ",IDnew)
            $PRINT if length(ListOfGenotypesSp)!==length(Map_ID_AllGenot)
            $PRINT      global ListOfGenotypesSp_ = deepcopy(ListOfGenotypesSp)
            $PRINT      global MapOfMutationalPathSp_ = deepcopy(MapOfMutationalPathSp)
            $PRINT      error("length(ListOfGenotypesSp)!==length(Map_ID_AllGenot)")
            $PRINT end
            $PRINT println("Update_MapOfMutationalPathSp 4.4, length(Map_ID_AllGenot[ID   ]) = ",length(Map_ID_AllGenot[ID   ]))
            Map_ID_AllGenot[ID      ][         dir[]   ][dist[] ] = IDnew[]
            $PRINT println("Update_MapOfMutationalPathSp 4.6, length(Map_ID_AllGenot[IDnew]) = ",length(Map_ID_AllGenot[IDnew]))
            Map_ID_AllGenot[IDnew[] ][$((2,1))[dir[] ] ][dist[] ] = genot
            $PRINT println("Update_MapOfMutationalPathSp 5")
        end
    end
    $PRINT println("Update_MapOfMutationalPathSp 6")
end"""
eval(Meta.parse(str))




# SUPR
# X = HolobiontMetaPop.nodes[2].nodes[2].nodes[2]
# NMu = HolobiontMetaPop.nodes[1].NMu
# NMu_expect = HolobiontMetaPop.nodes[1].NMu_expect 
# MapOfMutationalPathSp = deepcopy(MapOfMutationalPath[X.sp])
# ListOfGenotypesSp = Param.ListOfGenotypes[X.sp]
# ListOfPhenotypesSp = Param.ListOfPhenotypes[X.sp]
# NMu[] = 120.765
# d_nodes = dXX.nodes
# IndexNodesSp = XX.nodes[1].IndexNodesBySp[X.sp]
# nMu = Param.nMu
# CurPosi = Param.CurPosi
# NewPosi = Param.NewPosi
# NewVal = Param.NewVal
# IDnew = Param.IDnew
# CellGenotMuEnd = Param.CellGenotMuEnd
#  Head = HolobiontMetaPop
# dHead = similar(HolobiontMetaPop)
# Posi = XX.nodes[1].PosiEnvir
# DistSucepSp = Param.AllelesParam.SpecificityMatrix.DistSucep[X.sp]
# DistAggressSp = Param.AllelesParam.SpecificityMatrix.DistAggress[X.sp]
# 
# 
# trait = 50
# trait = 53


#                       #      |CellGenot|   dN_E             |Container      |Container        |MapMutatingTraits               |            Posi     |      LocalComp             |    d_LocalComp              |      NMu    |         NMu_expect   |  cell         |      MapOfMutationalPathSp      |      ListOfGenotypesSp               |      ListOfPhenotypesSp               |      DistBetweenGenotypesSp    |                    DistAggress                       |                  DistSucepSp                       |         IndexNodesSp                                                    |      IDnew      |      nMu      |      CurPosi      |      NewPosi      |      NewVal         |      CellGenotMuEnd|dt
#                       Mutate!(    X    ,dX.values           ,   Head        ,  dHead          ,X.MapMutatingTraitsCell/HostOnly,XX.nodes[1].PosiEnvir,         XX                 ,        dXX                  ,Param.NMu    ,   Param.NMu_expect   ,Val(true)      ,Param.MapOfMutationalPath[X.sp]  ,Param.ListOfGenotypes[X.sp]           ,Param.ListOfPhenotypes[X.sp]           ,Param.DistBetweenGenotypes[X.sp],Param.AllelesParam.SpecificityMatrix.DistAggress[X.sp],Param.AllelesParam.SpecificityMatrix.DistSucep[X.sp],XX.nodes[1].IndexNodesBySp[X.sp]                                         ,Param.IDnew      ,Param.nMu      ,Param.CurPosi      ,Param.NewPosi      ,Param.NewVal         ,Param.CellGenotMuEnd,dt)
  str = """ function Mutate!(X::CellGenot, dN_E::Vector{Float},Head::Container, dHead::Container,MapMutatingTraits::NamedTuple   , Posi::Vector{Int64} , XX::AbstractMultiScaleArray, dXX::AbstractMultiScaleArray,NMu::Ref{Int},NMu_expect::Ref{Float},cell::Val{true},MapOfMutationalPathSp::NamedTuple, ListOfGenotypesSp::Vector{NamedTuple}, ListOfPhenotypesSp::Vector{NamedTuple}, DistBetweenGenotypesSp::Vector , DistAggressSp::Tuple                                 , DistSucepSp::Tuple                                 , IndexNodesSp::NamedTuple{(:Posi, :ID),Tuple{Array{Int,1},Array{Int,1}}}, IDnew::Ref{Int}, nMu::Ref{Int}, CurPosi::Ref{Int}, NewPosi::Ref{Int}, NewVal::Ref{Integer}, CellGenotMuEnd::Ref{Union{Int,Nothing}}, dt::Float)
    # tuple(for each sp :NamedTuple(for each mutating trait Array(for each genoype.ID two Array(mutating by deacreasing.increasing trait value; position 1, 2, 3, ... = ID of the genotypes)) )))
        # Choose the number of traits
        $PRINT println("Mu -2")
        NMu_expect[] = X.values[1]*X.MuTotCell*dt
        $PRINT println("Mu -1")
        if NMu_expect[] < 1e16
            NMu[]=Int(rand(Poisson(NMu_expect[])))
        else
            NMu[]=Int(round(rand(Normal(NMu_expect[],sqrt(NMu_expect[])),1)[1]))
        end
        if !iszero(NMu[])
            $PRINT println("Mu 0")
            # Choose the mutations events
            if cell
                MutatingTraits     = rand(Multinomial(NMu[],MapMutatingTraits.Prob)) # number of mutations in each MutatingTraits
            else
                MutatingTraits     = rand(Multinomial(NMu[],MapMutatingTraits.Prob)) # number of mutations in each MutatingTraits
            end
            DireMutatingTraits = [nMu == 0 ? nothing : rand(Multinomial(nMu,[0.5,0.5])) for (i,nMu) in enumerate(MutatingTraits)] # Number of mutations that decrease or that increase trait value
            # une serie de distance pour chaque directions
            DistMutatingTraits = [isnothing(DireMutatingTraits[i]) ? nothing : [rand(Multinomial(nMu,$(AllelesParam.DistrMut)[MapMutatingTraits.traits[i]])) for nMu in DireMutatingTraits[i]] for i in eachindex(MutatingTraits)]
            #
            $PRINT println("Mu 1")
            for trait in Iterators.filter(T -> MutatingTraits[T]!==0 ,1:length(MutatingTraits))
                $PRINT println("Mu 2, trait = ",trait)
                pathAll = MapMutatingTraits.PathAll[trait]
                path = MapMutatingTraits.PathToTrait[trait]
                # tuple(for each sp :NamedTuple(for each mutating trait Array(for each genoype.ID two Array(mutating by deacreasing.increasing trait value; position 1, 2, 3, ... = ID of the genotypes)) )))
                Map_ID = MapOfMutationalPathSp
                PrevVal = X.alleles
                println(path)
                println(pathAll)
                for k in pathAll
                    println("k = ",k)
                    PrevVal = PrevVal[k]
                    Map_ID = Map_ID[k]
                    println("Map_ID = ", Map_ID)
                end
                $PRINT println("Mu 2.25")
                println("length(Map_ID) = ",length(Map_ID))
                println("X.ID = ",X.ID)
                Map_ID = Map_ID[X.ID]
                $PRINT println("Mu 2.5")
                k = path[end]
                $PRINT println("Mu 2.75")
                $PRINT println("Mu 3, path = ",path,", k = ",k)
                for dir in Iterators.filter(D -> DireMutatingTraits[trait][D]!==0 ,1:length(DistMutatingTraits[trait]))
                    for dist in Iterators.filter(D -> DistMutatingTraits[trait][dir][D]!==0 ,1:length(DistMutatingTraits[trait][dir]))
                        $PRINT println("Mu 4, dir = ",dir," dist = ",dist)
                        nMu[] = DistMutatingTraits[trait][dir][dist]
                        $PRINT println("Mu 4.33")
                        if any(pathAll .== :Mu)
                            CanMutateTo = $(AllelesParam.SpeciesCanMutateTo).Mu[X.sp]
                        else
                            CanMutateTo = $(AllelesParam.SpeciesCanMutateTo)[k][X.sp]
                        end
                        $PRINT println("Mu 4.66")
                        println("k = ",k)
                        println("PrevVal = ", PrevVal)
                        println("CanMutateTo = ", CanMutateTo)
                        println(findfirst(x -> PrevVal == x,  CanMutateTo))
                        CurPosi[] = findfirst(x -> PrevVal == x,  CanMutateTo)
                        $PRINT println("Mu 5")
                        if (dir == 2) NewPosi[] = CurPosi[]+dist else NewPosi[] = CurPosi[]-dist end
                        if (NewPosi[] > 0) & (NewPosi[] <= $(AllelesParam.Nalleles)[k][X.sp] ) # if we are inside of the range of possible allele value according to the simulation parameters
                            $PRINT println("Mu 6")
                            if (!((k==:TypeAggress) | (k==:TypeSucepti))) | any(pathAll .== :Mu)
                                println("NewPosi = ",NewPosi[])  ;   println("length(CanMutateTo) = ",length(CanMutateTo))  ;   println("Nalleles = ",AllelesParam.Nalleles[k][X.sp])
                                NewVal[] = CanMutateTo[NewPosi[]]
                            elseif (k==:TypeSucepti)
                                println("CurPosi = ",CurPosi[])  ;   println("PrevVal = ", PrevVal)  ;   println("CanMutateTo = ", CanMutateTo)  ;   println("k = ",k)  ;   println("length(DistSucepSp) = ",length(DistSucepSp))  ;   println("Nalleles = ",AllelesParam.Nalleles[k][X.sp])  ;   println("dist = ",dist)  ;   println("length(DistSucepSp[CurPosi[]]) = ",length(DistSucepSp[CurPosi[]]))
                                NewVal[] = DistSucepSp[CurPosi[]][dist]
                            else
                                println("CurPosi = ",CurPosi[])  ;   println("PrevVal = ", PrevVal)  ;   println("CanMutateTo = ", CanMutateTo)  ;   println("k = ",k)  ;   println("length(DistAggressSp) = ",length(DistAggressSp))  ;   println("Nalleles = ",AllelesParam.Nalleles[k][X.sp])  ;   println("dist = ",dist)  ;   println("length(DistAggressSp[CurPosi[]]) = ",length(DistAggressSp[CurPosi[]]))
                                NewVal[] = DistAggressSp[CurPosi[]][dist]
                            end
                            $PRINT println("Mu 7")
                            println(trait)  ;   println(k)  ;   println(path)  ;   println(pathAll)  ;   println(dir)  ;   println(dist)  ;   println(length(Map_ID))
                            #                             println(Map_ID)
                            println(length(Map_ID[dir]))
                            if iszero(Map_ID[dir][dist]) # we need to create the new genotype and update everything (ListOfGenotypesSp, MapOfMutationalPathSp and IndexNodesBySp)
                                $PRINT println("Mu 8")
                                NewAlleles = Mutate_Alleles(X.alleles, k, PrevVal, NewVal[] ,pathAll)
                                $PRINT println("Mu 9")
                                push!(ListOfGenotypesSp, NewAlleles )
                                if (length(ListOfGenotypesSp) > (MaxIntValues[Symbol(typeof(X.ID))] - 2)) error("The simulation has reached a too high number of genotypes.\nWhen defining 'struct CellHostGenot' and 'struct CellSymbGenot', you need to use an other type of integer.\nRefere to the object 'MaxIntValues' for choosing.") end
                                $PRINT println("Mu 9.333")
                                push!(ListOfPhenotypesSp, MakePhenotypesFromAlleles(NewAlleles,X.sp) )
                                $PRINT println("Mu 9.666, sp = ",X.sp)
                                #           Update MapOfMutationalPathSp
                                IDnew[] = length(ListOfGenotypesSp)
 # For each genotype that differ from X.ID only for 'trait' (and that are thus in Map_ID_AllGenot), add their position in Map_ID_AllGenot[IDnew[]] and add the position of IDnew[] in other genotypes
                                Update_MapOfMutationalPathSp!(MapOfMutationalPathSp,ListOfGenotypesSp,X.ID, IDnew, NewAlleles, MapMutatingTraits.PathAll, cell)
                                $PRINT println("Mu 10")
                                # Estimate the distance between the new genotype and the other genotypes
#                                 println(tuple(map(genot -> EstimateDistBetweenGenotypes(NewAlleles,genot,Val(X.sp)) ,ListOfGenotypesSp)...))
                                push!(DistBetweenGenotypesSp,
                                                             tuple(map(genot -> EstimateDistBetweenGenotypes(NewAlleles,genot,Val(X.sp)) ,ListOfGenotypesSp)...) )
                                $PRINT println("Mu 10.5")
                            else
                                $PRINT println("Mu 11")
                                IDnew[] = Map_ID[dir][dist]
                                NewAlleles = ListOfGenotypesSp[IDnew[]]
                                $PRINT println("Mu 11")
                            end
                            #
                            $PRINT println("Mu 11.333")
                            $PRINT println(isnothing(CellGenotMuEnd[]))
                            if !isnothing(CellGenotMuEnd[]) println(CellGenotMuEnd[]) end
                            CellGenotMuEnd[] = findfirst(id -> id == IDnew[],IndexNodesSp.ID)
                            $PRINT println("Mu 11.666")
                            $PRINT println(isnothing(CellGenotMuEnd[]))
                            if !isnothing(CellGenotMuEnd[]) println(CellGenotMuEnd[]) end
                            if isnothing(CellGenotMuEnd[]) # is the exiting genotype absent in the local compartment
                                $PRINT println("Mu 12")
                                $PRINT println("length(ListOfPhenotypesSp) = ",length(ListOfPhenotypesSp), ", IDnew[] = ",IDnew[])
                                if X.sp == 1
                                    $PRINT          GetMuTot(NewAlleles, X.sp , Val(true ))
                                    $PRINT println("GetMuTot(NewAlleles, X.sp , Val(true ))")
                                    $PRINT          GetMuTot(NewAlleles, X.sp , Val(false))
                                    $PRINT println("GetMuTot(NewAlleles, X.sp , Val(false))")
                                    $PRINT          MakeMapMutatingTraits(NewAlleles, true)
                                    $PRINT println("MakeMapMutatingTraits(NewAlleles, true)")
                                    $PRINT          MakeMapMutatingTraits(NewAlleles, false)
                                    $PRINT println("MakeMapMutatingTraits(NewAlleles, false)")
                                    cg = CellHostGenot([1.0 ,X.values[2:3]...] # N_E_Aging
                                        ,:aCellHostGenot # Iam
                                        ,IDnew[]
                                        ,X.sp
                                        ,NewAlleles # allele
                                        ,ListOfPhenotypesSp[IDnew[]] # pheno
                                        ,GetMuTot(NewAlleles, X.sp , Val(true )) # MuTotCell
                                        ,GetMuTot(NewAlleles, X.sp , Val(false)) # MuTotHostOnly
                                        ,MakeMapMutatingTraits(NewAlleles, true) # MapMutatingTraitsCell
                                        ,MakeMapMutatingTraits(NewAlleles, false) # MapMutatingTraitsHostOnly
                                        )
                                else
                                    $PRINT          GetMuTot(NewAlleles, X.sp , Val(true ))
                                    $PRINT println("GetMuTot(NewAlleles, X.sp , Val(true ))")
                                    $PRINT          MakeMapMutatingTraits(NewAlleles, true)
                                    $PRINT println("MakeMapMutatingTraits(NewAlleles, true)")
                                    cg = CellSymbGenot([1.0 ,X.values[2]] # N_E
                                        ,:aSymbGeno
                                        ,IDnew[]
                                        ,X.sp
                                        ,NewAlleles # allele
                                        ,ListOfPhenotypesSp[IDnew[]] # pheno
                                        ,GetMuTot(NewAlleles, X.sp , Val(true )) # MuTotCell
                                        ,MakeMapMutatingTraits(NewAlleles, true) # MapMutatingTraitsCell
                                        )
                                end
                                $PRINT println("Mu 13")
                                $PRINT PREV = deepcopy(XX)
                                add_node!(dHead, cg, Posi...) # add_node! the new cell genotype 
                                add_node!( Head, cg, Posi...) # add_node! the new cell genotype 
                                #
                                push!(IndexNodesSp.Posi,length(XX.nodes))
                                push!(IndexNodesSp.ID  ,IDnew[])
                                $PRINT if length(unique(IndexNodesSp.Posi)) < length(IndexNodesSp.Posi)
                                    $PRINT println("typeof.(PREV.nodes) = ",typeof.(PREV.nodes),", typeof.(XX.nodes) = ",typeof.(XX.nodes))
                                    $PRINT println("length(PREV.nodes) = ",length(PREV.nodes),", length(XX.nodes) = ",length(XX.nodes))
                                    $PRINT println("PREV.nodes[1].IndexNodesBySp[X.sp].Posi = ",PREV.nodes[1].IndexNodesBySp[X.sp].Posi,", IndexNodesSp.Posi = ",IndexNodesSp.Posi)
                                    $PRINT error("In Mutate!,   length(unique(IndexNodesSp.Posi)) < length(IndexNodesSp.Posi)")
                                $PRINT end
                                $PRINT println("Mu 14")
                            else
                                $PRINT println("Mu 15")
                                dN_E_ = dXX.nodes[IndexNodesSp.Posi[CellGenotMuEnd[]]].values
                                N_E_ =  XX.nodes[IndexNodesSp.Posi[CellGenotMuEnd[]]].values
                                #
                                dN_E_[1] += nMu[]
                                dN_E_[2] += #                  E(t+dt)                                -                       E(t) = 
                                            (N_E_[1]*N_E_[2] + nMu[]*X.values[2]) / (N_E_[1] + nMu[])  -  XX.nodes[IndexNodesSp.Posi[CellGenotMuEnd[]]].values[2]
                                $PRINT println("Mu 16")
                            end
                            dN_E[1] -= nMu[] 
                            $PRINT println("Mu 17")
                        end
                    end
                end
            end
        end
end"""
eval(Meta.parse(str))




###############
str = 
"""
function ResetsToZero!(X::AbstractMultiScaleArray)
    for i in eachindex(X.values) X.values[i] = 0.0 end
    for x in X.nodes
        if typeof(x)<:AbstractMultiScaleArrayLeaf
            for i in eachindex(x.values) x.values[i] = 0.0 end
        else
            ResetsToZero!(x)
        end
    end
    $PRINT if any(isnan.(flat([x.values for x in X.nodes]))) println("Which as NaN ?", flat([any(isnan.(x.values)) ? typeof(x) : "" for x in X.nodes])) ; println("flat([x.values for x in X.nodes]) = ", flat([x.values for x in X.nodes])) ; error("output flat([x.values for x in X.nodes]) has NaN (ResetsToZero)") end
end"""
eval(Meta.parse(str))

