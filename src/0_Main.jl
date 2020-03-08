# /DISQUE/FAC/Julia/julia-1.0.0-linux-x86_64/julia-1.0.0/bin/julia
# 
# exit()
/DISQUE/FAC/Julia/julia-1.2.0-linux-x86_64/julia-1.2.0/bin/julia
# include("/DISQUE/0_DefensiveSymbiosis/Model_2/The_model_in_julia/Main.jl")


# using DifferentialEquations
using MultiScaleArrays

# using Plots; plotly() # Using the Plotly backend
using Distributions
using StatsBase
using MultivariateStats

using LinearAlgebra
using StatsPlots # density plot of Beta
using Plots # Using the Plotly backend

# using StaticArrays
# using RecursiveArrayTools

Float = Float64 # precision of the floats
Int = Int64 # precision of the integers

# BOUND = "@inbounds"  # @inbounds eliminates array bounds checking within expressions. This faster the code
# or
BOUND = ""  # BOUND = "" ensure the normal secure behaviour

PRINT = "" # <= make verbose function (useful for debugging)
# or
# PRINT = "#"  # <= make silence function


# PlotToChooseParameters = true
PlotToChooseParameters = false

cd("/DISQUE/0_DefensiveSymbiosis/Model_2/The_model_in_julia/")

include("1_Functions_Helping_Setting_Parameters.jl")
include("2_Set_Parameters.jl")
include("3_Set_Alleles.jl")
include("4_Make_Phenotypes.jl") # this take a while when the number of allele per trait or and the number of species are large

include("5_Make_initial_state.jl") # Create the initial holobiont (meta)population

# You now have the initial holobiont (meta)population which is set according to your choices (and default values).
#### The hierarchical structure :
#
#              Container
# Parameters__/     |
#              MainEnvironment
#                /     \
#           SubEnvir  SubEnvir
#              /         \
#         CellGenotype  Hosts      # in the function evolve, some of these object only need to have been an obj and some an OBJ. This is mean that the 'Environment' needs to have been an OBJ. This is achived by the use of the level 'Container'.
#                         |
#                        ...
#                         |
#                   CellGenotype
#
####################################
# You can visualize using the function print_human_readable
# ex 1 :
println(HolobiontMetaPop ; NcharPerName = 5, LevelMax = 4)
# ex 2 :
println(HolobiontMetaPop ; NcharPerName = 1, LevelMax = 4)
# And you can tune it more finelly using the functions of the library "MultiScaleArrays". For more details on these functions, please refer to https://github.com/JuliaDiffEq/MultiScaleArrays.jl

include("6_evolve.jl")

Time = 0.0:1.0:10.0
#                                  W        | dt  |Duration|SaveAt|
EvoTraject = EvolveOverTime(HolobiontMetaPop, 0.01,  10.0  , Time )

sp = 2
length(HolobiontMetaPop.nodes[1].MapOfMutationalPath[sp])

for i in 1:length(X.nodes[2].nodes)
    plot!(Time, [X.nodes[2].nodes[i].nodes[2].values[1] for X in EvoTraject])
end
plot!()

plot()
for i in 1:length(X.nodes[2].nodes)
    plot!(Time, [X.nodes[2].nodes[i].nodes[2].values[2] for X in EvoTraject])
end
plot!()

plot()
for i in 1:SimParam.SubEnvir.NsubenvironmentsBySpecies
    plot!(Time, [X.nodes[2].nodes[i].nodes[2].values[3] for X in EvoTraject])
end
plot!()



# include("3_GeneralFunctions.jl")
# include("addition_deletion_IndexesAsList.jl")











































###################################################################
# Solve
# evolve(HolobiontMetaPop,2.0,1.0)

# NO discreteEvents
tspan = (0.0,8.0)
SaveAt=seq(tspan[1],tspan[2],10)
prob = ODEProblem(evolve, HolobiontMetaPop, tspan)


Sol = solve(prob
    #ODE Solvers
#     ,BS3()
    ,KuttaPRK2p5() # A 5 parallel, 2 processor explicit Runge-Kutta method of 5th order.

#     ,callback=Callbacks
    ,dense=false
    ,saveat=SaveAt
    ,save_everystep=false
    ,abstol=1 # 1e-2
    ,dt=0.01 # Sets the initial stepsize. This is also the stepsize for fixed timestep methods. Defaults to an automatic choice if the method is adaptive.
    ,dtmax=1000 # Maximum dt for adaptive timestepping. Defaults are package-dependent.
    ,dtmin=1e-2 # Minimum dt for adaptive timestepping. Defaults are package-dependent.
# force_dtmin = false # Declares whether to continue, forcing the minimum dt usage. Default is false, which has the solver throw a warning and exit early when encountering the minimum dt. Setting this true allows the solver to continue, never letting dt go below dtmin (and ignoring error tolerances in those cases). Note that true is not compatible with most interop packages.
)


[ [Sol.u[i][ii] for ii in eachindex(Sol.u[i]) ] for i in eachindex(Sol.u)]


[ Sol.u[i].nodes for i in eachindex(Sol.u)]

#                  Main SubE     Cell       1=N 2=E 3=Aging
plot([Sol.u[i].nodes[2].nodes[1  ].nodes[2].values[1] for i in eachindex(Sol.u)]) 
for sub in 2:length(Sol.u[1].nodes[2].nodes)
plot!([Sol.u[i].nodes[2].nodes[sub].nodes[2].values[1] for i in eachindex(Sol.u)]) 
end
plot!()

#                  Main SubE     Cell       1=N 2=E 3=Aging
plot([Sol.u[i].nodes[2].nodes[1  ].nodes[2].values[2] for i in eachindex(Sol.u)]) 
for sub in 2:length(Sol.u[1].nodes[2].nodes)
plot!([Sol.u[i].nodes[2].nodes[sub].nodes[2].values[2] for i in eachindex(Sol.u)]) 
end
plot!()

typeof(Sol.u[i].nodes[2].nodes[1])













integrator = init(prob, BS3())


# tuple([u for u in  integrator_1.u]...)
# integrator = deepcopy(integrator_1)

integrator.dt
integrator_1 = deepcopy(integrator)
solve!(integrator; alg_hints=[:additive # Denotes that the underlying SDE has additive noise.
                             ,
                             ,
                             ])
step!(integrator)
tuple([integrator.u[i]<0 ? i : "" for i in eachindex(integrator.u)]...)
tuple([u for u in  integrator.u]...)

#                   Envir   SubEnvir  Cell
typeof(integrator.u.nodes[2].nodes[1])
       integrator.u.nodes[2].nodes[1].values
typeof(integrator.u.nodes[2].nodes[1].nodes[2])
       integrator.u.nodes[2].nodes[1].nodes[2].values

















### ### ### ### ### ### ### ### ### ### ### ### ### ###
### CREATE THE HOLOBIONT POPULATION (the object that will evolve)










# Create an female host cell genotype
afemalehostcellgenotype = CreateACellGenotype(1,Nresources,NsymbiontSp,true)
afemalehostcellgenotype.SexF
# 1-element Array{Bool,1}:
#  true

# Create a male host cell genotype
amalehostcellgenotype = CreateACellGenotype(1,Nresources,NsymbiontSp,false)
amalehostcellgenotype.SexF
# 1-element Array{Bool,1}:
# false

# Create some symbiont genotype of species 1 and 2
asymbiont2cellgenotype_a = CreateACellGenotype(2,Nresources,NsymbiontSp)
asymbiont2cellgenotype_a.Vmax[2][1] = 2
asymbiont2cellgenotype_a.Km[2][1]   = 2

println(asymbiont2cellgenotype_a)

asymbiont2cellgenotype_b = CreateACellGenotype(2,Nresources,NsymbiontSp)
asymbiont2cellgenotype_b.Vmax[1][1] = 2
asymbiont2cellgenotype_b.Km[1][1]   = 2

asymbiont3cellgenotype = CreateACellGenotype(3,Nresources,NsymbiontSp)

[asymbiont2cellgenotype_a.Species, asymbiont2cellgenotype_a.Iam]
# 2-element Array{Any,1}:
#  UInt16[0x0002]
#  :CellSymbiontGenotype

[asymbiont2cellgenotype_b.Species, asymbiont2cellgenotype_b.Iam]
# 2-element Array{Any,1}:
#  UInt16[0x0002]
#  :CellSymbiontGenotype

[asymbiont3cellgenotype.Species, asymbiont3cellgenotype.Iam]
# 2-element Array{Any,1}:
#  UInt16[0x0003]
#  :CellSymbiontGenotype

# Create the property of the external environment 1
PropEnvirHosts = PropEnvir(
    cat(1,ConcArriving...,[0.0 for i in 1:Nalleles.Allelo]) # values: vector of length Nresources: Current NUMBER of molecules in the environment. A meaningful default starting value: ConcArriving
    ,:PropEnvir # Iam
    ,1e7 # Ve: volume of the environment
    ,1e7 / 1e3 * 1.0 # VeInOut: rate volume renewing in the of the main environment   make it proportionnal to the main environment
    ,SVector(ConcArriving) # ConcArriving: Concentration(s) of the resource(s) arriving in the environment (in a volume VeInOut per time step)
)
# Create the property of the external environment 2
PropEnvirSymbiont_2 = deepcopy(PropEnvirHosts) # at start, they are the same
# Create the property of the external environment 3
PropEnvirSymbiont_3 = deepcopy(PropEnvirHosts) # at start, they are the same

MutationRates = # Tuple of length 1 + nrb of symbiont species
    # The list of mutation rate must have the very same structure as "CellGenotype"s (exepte for the values). Non mutating traits (such as Iam or Species) must have the value "nothing"
    (
        (# Mutations rate of the host.
          nothing # values
        , nothing # I am
        , nothing # rE
        , nothing # d
        , nothing # Vmax
        , nothing # Km
        , nothing # PrE
        , nothing # ConstPr
        , nothing # SpePr
        , nothing # Detec
        , nothing # Allelo
        , nothing # Resist
        , nothing # S
        , nothing # DeadDecomp
        , nothing # Satur
        , nothing # Vc
        , nothing # VeInOutPerCell
        , nothing # nch
        , nothing # VgH
        , nothing # Aging
        , nothing # SizeDeath
        )
        ,
        (# Mutations rate of the symbiont species 1.
        nothing, # values
        nothing, # I am
        nothing , # rE
        nothing , # d
        (nothing, 1e-9) , # Vmax   Set the second value to 0 if their is only one resource, it wont be used any way
        (1e-7  ,nothing) , # Km     Set the second value to 0 if their is only one resource, it wont be used any way
        nothing , # PrE
        # Production of molecules of allelopathic compounds per individual
        nothing , # ConstPr Constitutive production of the allelopathic compound
        nothing , # SpePr   Specific production in response to a the presence of a species ; must be array Tuple, of size (NsymbiontSp+1)
        nothing , # Detec   ability to detetect a species meaningful range of values : -5, 5 ; must be array Tuple, of size (NsymbiontSp+1)
        nothing , # Allelo : a tuple of integrers in th range of 1:Nalleles.Allelo     giving the set of alleles to which this species can mutate
        nothing , # Resist : a tuple of integrers in th range of 1:Nalleles.Resist giving the set of alleles to which this species can mutate
        nothing , # S
        nothing , # D1
        nothing   # Species # 1 => host ; 2, 3, ... => Symbioint Species
        )
        ,
        (# Mutations rate of the symbiont species 2.
        nothing, # values
        nothing, # I am
        nothing , # rE
        nothing , # d
        (nothing, 1e-7) , # Vmax   Set the second value to 0 if their is only one resource, it wont be used any way
        (1e-7  ,nothing) , # Km     Set the second value to 0 if their is only one resource, it wont be used any way
        nothing , # PrE
        # Production of molecules of allelopathic compounds per individual
        nothing , # ConstPr Constitutive production of the allelopathic compound
        nothing , # SpePr   Specific production in response to a the presence of a species ; must be array Tuple, of size (NsymbiontSp+1)
        nothing , # Detec   ability to detetect a species meaningful range of values : -5, 5 ; must be array Tuple, of size (NsymbiontSp+1)
        nothing , # Allelo : a tuple of integrers in th range of 1:Nalleles.Allelo     giving the set of alleles to which this species can mutate
        nothing , # Resist : a tuple of integrers in th range of 1:Nalleles.Resist giving the set of alleles to which this species can mutate
        nothing , # S
        nothing , # D1
        nothing   # Species # 1 => host ; 2, 3, ... => Symbioint Species
        )
    )


# check if all the traits that can have a mutation rate that is not 'nothing' have several alleles
Problem = [[] for i in 1:(NsymbiontSp+1)]
for sp in 1:(NsymbiontSp+1)
    for trait in filter(i -> !(MutationRates[sp][i]==nothing) , eachindex(MutationRates[sp]))
        if length(MutationRates[sp][trait]) > 1
            for trait2 in filter(i -> !(MutationRates[sp][trait][i]==nothing) , eachindex(MutationRates[sp][trait]))
                if length(getfield( Alleles[sp],trait)[trait2]) == 1
                    push!(Problem[sp],Symbol(fieldnames( Alleles[sp])[trait],"_",trait2))
                end
            end
        else
            if length(getfield( Alleles[sp],trait)) == 1
                push!(Problem[sp],fieldnames( Alleles[sp])[trait])
end;end;end;end
if !(Problem == [[] for i in 1:(NsymbiontSp+1)])
    error("You defined that some alleles can mutate (Non 'nothing' value) but there is only one value corresponding to this trait in the objet 'Allele'.\nThe problem is for the following traits:\n"join([
        "Species "string(sp)" Trait(s): "ifelse(Problem[sp]==[],"No problem detected",join(Problem[sp],", ")) for sp in 1:(NsymbiontSp+1)
    ],"\n"))
end
# END check if all the traits that can have a mutation rate that is not 'nothing' have several alleles

TraitsIndexies = tuple([linearIndexies(MutationRates[sp]) for sp in 1:(NsymbiontSp+1)]...) # get the position of the traits in the infofields
MutationRatesLinear = LineariesMutationRates(MutationRates)

MutationRatesAllTraits = Vector{Any}([0.0 for i in 1:(NsymbiontSp+1)]) # sum of the MutationRates of all mutating traits: probability of mutating a trait
MutationRatesAsWeight = Vector{Any}([Weights([1]) for i in 1:(NsymbiontSp+1)])
for sp = 1:(NsymbiontSp+1)
    if MutationRatesLinear[sp] == ()
        MutationRatesAllTraits[sp] = 0.0
        MutationRatesAsWeight[sp] = Weights([1])
    else
        MutationRatesAllTraits[sp] = sum(MutationRatesLinear[sp])
        MutationRatesAsWeight[sp] = Weights(collect(MutationRatesLinear[sp]))
    end
end
MutationRatesAllTraits = tuple(MutationRatesAllTraits...)
MutationRatesAsWeight = tuple(MutationRatesAsWeight...)

Mutations = MutationsParam(
    MutationRates
   ,TraitsIndexies
   ,MutationRatesLinear # MutationRates of each mutating trait in the order in which the mutating trais are in TraitsIndexies
   ,MutationRatesAllTraits # sum of the MutationRates of all mutating traits: probability of mutating a trait
   ,MutationRatesAsWeight
)

#########
### Testing if there are some mutating individuals at every time step, for each genotype and each trait would be way too slow
Events = DiscreteEvents(
# N : For eah species, for each genotypes, number of individuals
 MVector{NsymbiontSp+1,Vector{Float64}}(ntuple(i->[],  NsymbiontSp+1  )...)
# Nmutating : for each species, amount mutating *individual*
,MVector{NsymbiontSp+1,Int64}()
# RemoveGenotypes : for each species, is their some empty genotypes without any individual that should be removed ?
,MVector{NsymbiontSp+1,Vector{Bool}}(ntuple(i->[],  NsymbiontSp+1  )...)
)

# Create the property of the container
PropSimul = PropSimulation(
    [0.0] # Time since the start of the simulation
    ,:PropSimulation # Iam
    ## Parameters of the simulation that are commmon to the whole hollobiont population
    ,Int8(2) # number of type of resources used in the simulation
    ,Int8(NsymbiontSp) # NsymbiontSp
    ,Int8(Nalleles.Allelo) #
    ,SVector(1.0, 2.0) # E: amount of energy contained in each kind of resources
    ,deepcopy(IndexesPerSp0)  # IndexesPerSp.Nodes :  # for each species: positions where the species is found in the AbstractMultiScaleArray
                                            # List of size 1+NsymbiontSp; for each species, where are the leaf describing the abundance of a genotype. E.g.: [1,2,4] will be the position of the leaf HolobiontPop.nodes[1].nodes[2].nodes[4]
                                            # All start by 2 because of the "Container", which 1st element is the property
                                            # Then comes the ID of the sub-environment and the position in the sub-environment. The 1st position in the sub-environment (or the host) being reserved for the PropEnvir
    ,Alleles # Alleles: Tuple containing 3 objects of type "allele"
    ,SpecificityMatrix # Specificity of the allelopathic compounds
    ,CostProductionAllelo # Cost of producing an allelopathic compounds
    ,Limitingresources
    ,Mutations
    ,MVector{NsymbiontSp+1}([MVector{length(sp),Int64}([0 for i in sp]) for sp in MutationRatesLinear]) # MutationsCounter
    ,Events
    ,MVector{1,Bool}(false)# extinction
    ,tuple([PhenFun0 for sp in 1:(NsymbiontSp+1)]...) #PhenotypesFun
    ,VoidFun # UpDateNode!
    ,VoidFun # UpDateNodeAsString
)

# Create a youngest host
AHostFemale = CreateAHost0old(afemalehostcellgenotype ,[1.1,1.1],PropSimul) # eggcell is the cell from which the host is created,  ***Concentrations*** is a vector with the all the compounds (resources or allelopathic compounds)
AHostMale   = CreateAHost0old(amalehostcellgenotype   ,[1.1,1.1],PropSimul) # eggcell is the cell from which the host is created,  ***Concentrations*** is a vector with the all the

# # # # # # # # # # # At construction all nodes should contain a leaf to be able to accept one later. SAME FOR NODES SUBTYPES
# # # # # # # # # # # If you don't want one at start, you should crerat the node with the leaf, and then delete the leaf (with remove_node!)
HolobiontPop = construct(Container,[
    PropSimul
    ,construct(Environment,[ # the host environment
#         PropEnvirHosts, AHostFemale, AHostMale
        PropEnvirHosts, afemalehostcellgenotype, amalehostcellgenotype
        ])
    ,construct(Environment,[ # the symbiont 2 environment
        PropEnvirSymbiont_2, asymbiont2cellgenotype_a, asymbiont2cellgenotype_b])
    ,construct(Environment,[ # the symbiont 3 environment
        PropEnvirSymbiont_3, asymbiont3cellgenotype])
])
# IndexesFromList(HolobiontPop, [1])
# remove_node!(HolobiontPop,1)
# IndexesFromList(HolobiontPop, [1])

include("MakePhenotypes.jl")

# Get the Nodes Indexes
GetNodeIndexes!(HolobiontPop)
# Get the Full Indexes
GetFullIndexes!(HolobiontPop)
# Make the phenotypes
MakePhenotypeValAndFunWholePop!(HolobiontPop)
# Make the function that will evolve the population
MakeUpDateNode!(HolobiontPop)

# HolobiontPop.nodes[1].UpDateNodeAsString()
# HolobiontPop.nodes[1].UpDateNode!()

println(HolobiontPop)

# XX_CASE_XX = HolobiontPop
# X  = HolobiontPop
# X  = HolobiontPop.nodes[2]
# x  = X.nodes[2]
# dX = similar(HolobiontPop)
# PropSimul = HolobiontPop.nodes[1]
# HolobiontPop.nodes[2].nodes[3][1] = 30000.0
# AbsolutPosiNode = filter(j ->  HolobiontPop[j] == 30000.0, eachindex(HolobiontPop))[1]
# HolobiontPop.nodes[2].nodes[3][1] = 30.0


# X = HolobiontPop
# X.nodes[2].nodes[3].nodes[2].Phenotypes.IndividualConsumption(  X.nodes[2].nodes[3].nodes[1][1:Nresources]  )
# X.nodes[2].nodes[3].nodes[3].Phenotypes.IndividualConsumption(  X.nodes[2].nodes[3].nodes[1][1:Nresources]  )
# X.nodes[2].nodes[4].nodes[2].Phenotypes.IndividualConsumption(  X.nodes[2].nodes[3].nodes[1][1:Nresources]  )


####################################################################
## Solve
# evolve(HolobiontPop,2.0,1.0)
#
# # NO discreteEvents
# tspan = (0.0,8.0)
# prob = ODEProblem(evolve, HolobiontPop, tspan)
# integrator = init(prob, BS3())
#
#
# # tuple([u for u in  integrator_1.u]...)
# # integrator = deepcopy(integrator_1)
#
# integrator.dt
# integrator_1 = deepcopy(integrator)
# step!(integrator)
# tuple([integrator.u[i]<0?i:"" for i in eachindex(integrator.u)]...)
# tuple([u for u in  integrator.u]...)
#
# integrator.u.nodes[1].nodes[1]
#
# integrator.dt
#
# BS3(integrator)
#
# # sol = solve(prob,dtmax = 1.0)
# @time sol = solve(prob,BS3())
# # 1.867647 seconds (3.55 M allocations: 260.097 MiB, 5.43% gc time)
# # 0.211747 seconds (537.01 k allocations: 35.516 MiB, 4.68% gc time)
#
# @time sol = solve(prob,Tsit5())
# # 3.802899 seconds (7.49 M allocations: 569.266 MiB, 5.18% gc time)
# # 0.482009 seconds (1.27 M allocations: 84.971 MiB, 7.89% gc time)

integrator = init(ODEProblem(evolve, HolobiontPop, (0.0,50.0)), BS3())

# DiscreteEvents
include("DiscreteEvents.jl")
tspan = (0.0,50.0)
tspan = (0.0,Inf)
prob = ODEProblem(evolve, HolobiontPop, tspan)
Sol = solve(prob,BS3(),callback=Callbacks,dense=false,save_everystep=false)
# Sol = solve(prob,BS3())

println(Sol.u[end],fields=["values"])

env = 4
for c in Sol.u[end].nodes[env].nodes[2:end]
    println(c.PhenotypesVal)
end

# # FreqDiscreteEvents=0.1
# # tspan = (0.0,5)
# # discreteEvents = linspace(tspan[1], tspan[2], tspan[2]/FreqDiscreteEvents)
#
# # tspan = (0,FreqDiscreteEvents)
# # prob = ODEProblem(evolve, HolobiontPop, tspan)
#
# states = [HolobiontPop]
# state  = deepcopy(HolobiontPop)
# t = 0.0
# tEnd = 5.0
# while t < tEnd
#     println(t/tEnd)
#     tspan = (t,tEnd)
#     prob = ODEProblem(evolve, state, tspan)
#   # integrator = init(prob, Tsit5())
#     sol = solve(prob,callback=CBstopSolveToMakeDiscreteEvents)
#     println(DiscreteEvents["Nmutating"])
#     println(DiscreteEvents["RemoveGenotypes"])
#     DiscreteEvents!(sol.u[end])
#     t = sol.t[end]
#     state = sol.u[end]
#     push!(states,sol.u[end])
# end
#
# state.nodes[1].nodes[1].infofields[InfoPropEnvir.IndexesPerSp]
#
#
# # for s in states
# #     for c in s.nodes[1].nodes[3].nodes
# #         println(c.infofields)
# #     end
# # end
#
# [[c.infofields for c in s.nodes[1].nodes[2].nodes] for s in states]
# [[c.infofields for c in s.nodes[1].nodes[3].nodes] for s in states]
#
# s.nodes[1].nodes
