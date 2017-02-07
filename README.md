# MultiScaleArrays

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Build Status](https://travis-ci.org/JuliaDiffEq/MultiScaleArrays.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/MultiScaleArrays.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/y0mjys35k6rbntbv?svg=true)](https://ci.appveyor.com/project/ChrisRackauckas/multiscalearrays-jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaDiffEq/MultiScaleArrays.jl/badge.svg)](https://coveralls.io/github/JuliaDiffEq/MultiScaleArrays.jl)
[![codecov](https://codecov.io/gh/JuliaDiffEq/MultiScaleArrays.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaDiffEq/MultiScaleArrays.jl)

MultiScaleArrays.jl allows you to easily build multiple scale models which are
fully compatible with native Julia scientific computing packages like
DifferentialEquations.jl or Optim.jl. These models utilize
a tree structure to describe phenomena of multiple scales, but the interface allows
you to describe equations on different levels, using aggregations from lower
levels to describe complex systems. Their structure allows for complex and dynamic
models to be developed with only a small performance difference. In the end, they present
themselves as an `AbstractArray` to standard solvers, allowing them to be used
in place of a `Vector` in any appropriately made Julia package.

## Example

The usage is best described by an example. Here we build a hierarchy where
Embryos contain Tissues which contain Populations which contain Cells, and the
cells contain proteins whose concentrations are modeled as simply a vector
of numbers (it can be anything linearly indexable).

```julia
immutable Cell{B} <: MultiScaleArrayLeaf{B}
  x::Vector{B}
end
immutable Population{T<:AbstractMultiScaleArray,B<:Number} <: AbstractMultiScaleArray{B}
  x::Vector{T}
  y::Vector{B}
  end_idxs::Vector{Int}
end
immutable Tissue{T<:AbstractMultiScaleArray,B<:Number} <: AbstractMultiScaleArray{B}
  x::Vector{T}
  y::Vector{B}
  end_idxs::Vector{Int}
end
immutable Embryo{T<:AbstractMultiScaleArray,B<:Number} <: MultiScaleArrayHead{B}
  x::Vector{T}
  y::Vector{B}
  end_idxs::Vector{Int}
end
```

Using this we can directly construct leaf types:

```julia
cell1 = Cell([1.0;2.0;3.0])
cell2 = Cell([4.0;5])
```

and build types higher up in the hierarchy by using the `constuct` method. The method
is `construct(T::AbstractMultiScaleArray,x,y)`, though if `y` is not given it's
taken to be empty.

```julia
population = construct(Population,deepcopy([cell1;cell2])) # Make a Population from cells
cell3 = Cell([3.0;2.0;5.0])
cell4 = Cell([4.0;6])
population2 = construct(Population,deepcopy([cell3;cell4]))
tissue1 = construct(Tissue,deepcopy([population;population2])) # Make a Tissue from Populations
tissue2 = construct(Tissue,deepcopy([population2;population]))
embryo = construct(Embryo,deepcopy([tissue1;tissue2])) # Make an embryo from Tissues
```

The head node then acts as the king. It is designed to have functionality which
mimics a vector in order for usage in DifferentialEquations or Optim. So for example

```julia
embryo[12]
```

returns the "12th protein", counting by Embryo > Tissue > Population > Cell in order
of the vectors. The linear indexing exists for every `AbstractMultiScaleArray`.
These types act as full linear vectors, so standard operations do the sensical
operations:

```julia
embryo[10] = 4.0 # changes protein concentration 10
embryo[2,3,1] # Gives the 1st cell in the 3rd population of the second tissue
embryo[:] # generates a vector of all of the protein concentrations
eachindex(embryo) # generates an iterator for the indices
```

Continuous models can thus be written at the protein level and will work seamlessly
with DifferentialEquations or Optim which will treat it like a vector of protein concentrations.
Using the iterators, note that we can get each cell population by looping through
2 levels below the top, so

```julia
for cell in level_iter(embryo,2)
  # Do something with the cells!
end
```

To apply a function cell-by-cell, you can write a dispatch `f` on the type for the
level. Using `level_iter_idx`, we can have its changes update some other head node
`d_em` via:

```julia
for (cell,y,z) in level_iter_idx(embryo,2)
  f(t,cell,@view embryo[y:z])
end
```

(`@view` is a standard Julia construct for creating a view instead of a copy when slicing an array,
making things alittle faster when indexing to create a sub-array). Notice that this updates the top
vector cell-by-cell via the function `f` without allocating. This allows one to apply an ODE "cell-wise".

However, the interesting behavior comes from event handling. Since `em` will be the
"vector" for the differential equation or otimization problem, it will be the value
passed to the event handling. MultiScaleArrays includes behavior for changing the
structure. For example:

```julia
tissue3 = construct(Tissue,deepcopy([population;population2]))
add_daughter!(embryo,tissue3) # Adds a new tissue to the embryo
remove_daughter!(embryo,2,1) # Removes population 1 from tissue 2 of the embryo
```

Combined with event handling, this allows for dynamic structures to be derived from
low level behaviors.

## Idea

The idea behind MultiScaleArrays is simple. The `*DiffEq` solvers (OrdinaryDiffEq.jl,
StochasticDiffEq.jl, DelayDiffEq.jl, etc.) and native optimization packages like
Optim.jl in their efficient in-place form all work with any Julia-defined
`AbstractArray` which has a linear index. Thus, to define our multiscale model,
we develop a type which has an efficient linear index. One can think of representing
cells with proteins as each being an array with values for each protein. The linear
index of the multiscale model would be indexing through each protein of each cell.
With proper index overloads, one can define a type such that `a[i]` does just that,
and thus it will work in the differential equation solvers. MultiScaleArrays.jl
takes that further by allowing one to recursively define an arbitrary `n`-level
hierarchical model which has efficient indexing structures. The result is a type
which models complex behavior, but the standard differential equation solvers will
work directly and efficiently on this type, making it easy to develop novel models
without having to re-develop advanced adaptive/stiff/stochastic/etc. solving
techniques for each new model.

## Defining A MultiScaleModel: The Interface

The required interface is as follows. Leaf types must extend MultiScaleArrayLeaf, the
highest level of the model or the head extends MultiScaleModelHead, and all
intermediate types extend AbstractMultiScaleArray. The leaf has an array `x::Vector{B}`.
Each type above then contains three fields:

- `x::Vector{T}`
- `y::Vector{B}`
- `end_idxs::Vector{Int}``

`B` is the `BottomType`, which has to be the same as the eltype for the array
in the leaf types. `T` is another `AbstractMultiScaleArray`. Thus at each level,
an` AbstractMultiScaleArray` contains some information of its own (`y`), the
next level down in the heirarchy (`x`), and caching for indices (`end_idxs`).
You can add and use extra fields as you please, and even make the types immutable.

## The MultiScaleModel API

The resulting type acts as an array. A leaf type `l` acts exactly as an array
with `l[i]==l.x[i]`. Higher nodes also act as a linear array. If `ln` is level
`n` in the heirarchy, then `ln.x` is the vector of level `n-1` objects, and `ln.y`
are its "intrinsic values". There is an indexing scheme on `ln`, where:

- `ln[i,j,k]` gets the `k`th `n-3` obejct in the `j`th `n-2` object in the `i`th level `n-1`
  object. Of course, this recurses for the whole hierarchy.
- `ln[i]` provides a linear index through all `.x` and `.y` values in every lower
  level and `ln.y` itself.

Thus `typeof(ln) <: AbstarctArray{B,1}` where `B` is the eltype of its leaves and
all `.y`'s.

In addition, iterators are provided to make it easy to iterate through levels.
For `h` being the head node, `level_iter(h,n)` iterates through all level objects
`n` levels down from the top, while `level_iter_idx(h,n)` is an enumeration
`(x,y,z)` where `x` are the `n`th from the head objects, with `h[y:z]` being
the values it holds in the linear indexing.

### Extensions

Note that this only showed the most basic MultiScaleModel. These types can be
extended as one pleases. For example, we can change the definition of the cell
to have:

```julia
immutable Cell{B} <: MultiScaleArrayLeaf{B}
  x::Vector{B}
  celltype::Symbol
end
```

Then we'd construct cells with `cell3 = Cell([3.0;2.0;5.0],:BCell)`, and can
give it a cell type. This information is part of the call, so

```julia
for (cell,y,z) in level_iter_idx(embryo,2)
  f(t,cell,@view embryo[y:z])
end
```

can allow one to check the `cell.celltype` in `f` an apply a different ODE depending
on the cell type. You can add fields however you want, so you can use them
to name cells and track lineages.

Showing the use of `y`, you just pass it to the constructor. Let's pass it an array
of 3 values:

```julia
tissue = construct(Tissue,deepcopy([population;population2]),[0.0;0.0;0.0])
```

We can selectively apply some function on these `y` values via:

```julia
for (tissue,y,z) in level_iter_idx(embryo,1)
  f(t,tissue,@view embryo[y:z])
end
```

and mutatate `tis.y` in `f`. For example, we could have

```julia
function f(t,tissue::Tissue,du)
  du .+= randn(3)
end
```

applies normal random numbers to the three values. We could use this to add to the
model the fact that `tissue.y[1:3]` are the tissue's position, and `f` would then be
adding Brownian motion.

Of course, you can keep going and kind of do whatever you want. The power is yours!
