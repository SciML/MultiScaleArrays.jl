__precompile__()

module MultiScaleArrays

    import Base: length, push!, deleteat!, getindex, setindex!, eachindex, ndims,
        size, similar, hcat, vcat,
        ==, *, +, -, show, reshape

    import RecursiveArrayTools: recursivecopy, recursivecopy!
    using Base.Iterators: flatten

    """
    AbstractMultiScaleArray{B}

    AbstractMultiScaleArray{B} <: AbstractVector{B}

    Abstract supertype for the nodes of a hierarchical array. A concrete leaf subtype stores a
    `values` field, while non-leaf subtypes store `nodes`, `values`, and `end_idxs` fields in that
    order. The head node should additionally subtype [`AbstractMultiScaleArrayHead`](@ref), so it
    can be passed to [`add_node!`](@ref), [`remove_node!`](@ref), and `getindices`.

    # Interface

    - A leaf subtype must have `values::AbstractVector` as its first field.
    - A non-leaf subtype must have `nodes`, `values`, and `end_idxs::AbstractVector{Int}` as its
      first three fields. `end_idxs` stores the cumulative linear lengths of `nodes`, followed by
      the length including `values` when `values` is nonempty.
    - Extra fields may follow the required fields and may be mutable or immutable.

    # Examples

    ```julia
    struct Cell <: AbstractMultiScaleArrayLeaf{Float64}
        values::Vector{Float64}
    end

    struct Tissue <: AbstractMultiScaleArray{Float64}
        nodes::Vector{Cell}
        values::Vector{Float64}
        end_idxs::Vector{Int}
    end
    ```

    # Defining A MultiScaleModel: The Interface

    The required interface is as follows. Leaf types must extend AbstractMultiScaleArrayLeaf, the
    highest level of the model or the head extends MultiScaleModelHead, and all
    intermediate types extend AbstractMultiScaleArray. The leaf has an array `values::Vector{B}`.
    Each type above then contains three fields:

      - `nodes::Vector{T}`
      - `values::Vector{B}`
      - `end_idxs::Vector{Int}`

    Note that the ordering of the fields matters.
    `B` is the `BottomType`, which has to be the same as the eltype for the array
    in the leaf types. `T` is another `AbstractMultiScaleArray`. Thus, at each level,
    an `AbstractMultiScaleArray` contains some information of its own (`values`), the
    next level down in the hierarchy (`nodes`), and caching for indices (`end_idxs`).
    You can add and use extra fields as you please, and you can even make the types immutable.

    ## The MultiScaleModel API

    The resulting type acts as an array. A leaf type `l` acts exactly as an array
    with `l[i] == l.values[i]`. Higher nodes also act as a linear array. If `ln` is level
    `n` in the hierarchy, then `ln.nodes` is the vector of level `n-1` objects, and `ln.values`
    are its “intrinsic values”. There is an indexing scheme on `ln`, where:

      - `ln[i,j,k]` gets the `k`th `n-3` object in the `j`th `n-2` object in the `i`th level `n-1`
        object. Of course, this recurses for the whole hierarchy.
      - `ln[i]` provides a linear index through all `.nodes` and `.values` values in every lower
        level and `ln.values` itself.

    Thus, `ln isa AbstractVector{B}`, where `B` is the eltype of its leaves and
    all `.values`'s.

    In addition, iterators are provided to make it easy to iterate through levels.
    For `h` being the head node, `level_iter(h,n)` iterates through all level objects
    `n` levels down from the top, while `level_iter_idx(h,n)` is an enumeration
    `(node,y,z)` where `node` are the `n`th from the head objects, with `h[y:z]` being
    the values it holds in the linear indexing.

    ## Indexing and Iteration

    The head node then acts as the king. It is designed to have functionality which
    mimics a vector in order for usage in DifferentialEquations or Optim. So for example

    ```julia
    embryo[12]
    ```

    returns the “12th protein”, counting by Embryo > Tissue > Population > Cell in order
    of the vectors. The linear indexing exists for every `AbstractMultiScaleArray`.
    These types act as full linear vectors, so standard operations do the sensical
    operations:

    ```julia
    embryo[10] = 4.0 # changes protein concentration 10
    embryo[2, 3, 1] # Gives the 1st cell in the 3rd population of the second tissue
    embryo[:] # generates a vector of all of the protein concentrations
    eachindex(embryo) # generates an iterator for the indices
    ```

    Continuous models can thus be written at the protein level and will work seamlessly
    with DifferentialEquations or Optim which will treat it like a vector of protein concentrations.
    Using the iterators, note that we can get each cell population by looping through
    2 levels below the top, so

    ```julia
    for cell in level_iter(embryo, 3)
        # Do something with the cells!
    end
    ```

    or the multiple-level iter, which is the one generally used in
    DifferentialEquations.jl functions:

    ```julia
    for (cell, dcell) in LevelIter(3, embryo, dembryo)
        # If these are similar structures, `cell` and `dcell` are the similar parts
        cell_ode(dcell, cell, p, t)
    end
    ```

    `LevelIterIdx` can give the indices along with iteration:

    ```julia
    for (cell, y, z) in LevelIterIdx(embryo, 3)
        # cell = embryo[y:z]
    end
    ```

    However, the interesting behavior comes from event handling. Since `embryo` will be the
    “vector” for the differential equation or optimization problem, it will be the value
    passed to the event handling. MultiScaleArrays includes behavior for changing the
    structure. For example:

    ```julia
    tissue3 = construct(Tissue, deepcopy([population, population2]))
    add_node!(embryo, tissue3) # Adds a new tissue to the embryo
    remove_node!(embryo, 2, 1) # Removes population 1 from tissue 2 of the embryo
    ```

    Combined with event handling, this allows for dynamic structures to be derived from
    low-level behaviors.

    ## Heterogeneous Nodes via Tuples

    Note that tuples can be used as well. This allows for type-stable broadcasting with
    heterogeneous nodes. This could be useful for mixing types
    inside of the nodes. For example:

    ```julia
    struct PlantSettings{T}
        x::T
    end
    struct OrganParams{T}
        y::T
    end

    struct Organ{B <: Number, P} <: AbstractMultiScaleArrayLeaf{B}
        values::Vector{B}
        name::Symbol
        params::P
    end

    struct Plant{B, S, N <: Tuple{Vararg{Organ{<:Number}}}} <: AbstractMultiScaleArray{B}
        nodes::N
        values::Vector{B}
        end_idxs::Vector{Int}
        settings::S
    end

    struct Community{B, N <: Tuple{Vararg{Plant{<:Number}}}} <: AbstractMultiScaleArray{B}
        nodes::N
        values::Vector{B}
        end_idxs::Vector{Int}
    end

    mutable struct Scenario{B, N <: Tuple{Vararg{Community{<:Number}}}} <:
        AbstractMultiScaleArrayHead{B}
        nodes::N
        values::Vector{B}
        end_idxs::Vector{Int}
    end

    organ1 = Organ([1.1, 2.1, 3.1], :Shoot, OrganParams(:grows_up))
    organ2 = Organ([4.1, 5.1, 6.1], :Root, OrganParams("grows down"))
    organ3 = Organ([1.2, 2.2, 3.2], :Shoot, OrganParams(true))
    organ4 = Organ([4.2, 5.2, 6.2], :Root, OrganParams(1 // 3))
    plant1 = construct(Plant, (deepcopy(organ1), deepcopy(organ2)), Float64[], PlantSettings(1))
    plant2 = construct(
        Plant, (deepcopy(organ3), deepcopy(organ4)), Float64[],
        PlantSettings(1.0)
    )
    community = construct(Community, (deepcopy(plant1), deepcopy(plant2)))
    scenario = construct(Scenario, (deepcopy(community),))
    ```

    (of course at the cost of mutability).
    """
    abstract type AbstractMultiScaleArray{B} <: AbstractVector{B} end

    """
        AbstractMultiScaleArrayLeaf{B}

    Abstract supertype for leaf nodes in a multiscale array hierarchy.

    # Interface

    A concrete leaf must define `values` as its first field. Its elements provide the leaf's
    linear array entries. Use [`construct`](@ref) directly for a leaf only when its ordinary
    constructor accepts the same arguments.

    # Examples

    ```julia
    struct Cell <: AbstractMultiScaleArrayLeaf{Float64}
        values::Vector{Float64}
    end

    Cell([1.0, 2.0])
    ```
    """
    abstract type AbstractMultiScaleArrayLeaf{B} <: AbstractMultiScaleArray{B} end

    """
        AbstractMultiScaleArrayHead{B}

    Abstract supertype for the top node of a multiscale array hierarchy.

    # Interface

    A head has the same required fields as [`AbstractMultiScaleArray`](@ref). It is the entry
    point for structural changes: [`add_node!`](@ref), [`remove_node!`](@ref), and
    [`getindices`](@ref) accept a head and update or inspect the whole hierarchy.

    # Examples

    ```julia
    mutable struct ModelHead <: AbstractMultiScaleArrayHead{Float64}
        nodes::Vector{Tissue}
        values::Vector{Float64}
        end_idxs::Vector{Int}
    end
    ```
    """
    abstract type AbstractMultiScaleArrayHead{B} <: AbstractMultiScaleArray{B} end

    using DiffEqBase: DiffEqBase
    using Statistics: Statistics
    using LinearAlgebra: LinearAlgebra, ldiv!
    using FiniteDiff: FiniteDiff
    import OrdinaryDiffEq, OrdinaryDiffEqCore, OrdinaryDiffEqRosenbrock, StochasticDiffEq, ForwardDiff
    import OrdinaryDiffEqDifferentiation
    import SciMLBase
    using SciMLBase: full_cache, rand_cache
    import DifferentiationInterface as DI
    using PrecompileTools: @compile_workload, @setup_workload

    Base.show(io::IO, x::AbstractMultiScaleArray) = invoke(show, Tuple{IO, Any}, io, x)
    Base.show(io::IO, ::MIME"text/plain", x::AbstractMultiScaleArray) = show(io, x)

    include("shape_construction.jl")
    include("addition_deletion.jl")
    include("indexing.jl")
    include("math.jl")
    include("level_iterations.jl")
    include("diffeq.jl")
    include("print_human_readable.jl")

    # Types
    export AbstractMultiScaleArray, AbstractMultiScaleArrayLeaf,
        AbstractMultiScaleArrayHead

    # Constructors
    export construct

    # Addition Deletion
    export add_node!, remove_node!

    # Indexing
    export num_nodes, getindices

    # Misc
    export LevelIterIdx, LevelIter, level_iter

    # print_human_readable
    export print_human_readable

    struct _PrecompileCell <: AbstractMultiScaleArrayLeaf{Float64}
        values::Vector{Float64}
    end

    mutable struct _PrecompileHead{T <: AbstractMultiScaleArray} <:
        AbstractMultiScaleArrayHead{Float64}
        nodes::Vector{T}
        values::Vector{Float64}
        end_idxs::Vector{Int}
    end

    @setup_workload begin
        @compile_workload begin
            cell = _PrecompileCell([1.0, 2.0])
            head = construct(_PrecompileHead, [cell])
            head[1]
            head[2] = 3.0
            num_nodes(head)
            collect(level_iter(head, 1))
            add_node!(head, _PrecompileCell([4.0]))
            remove_node!(head, 2)
        end
    end

end # module
