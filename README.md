# MultiScaleModels

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Build Status](https://travis-ci.org/JuliaDiffEq/MultiScaleModels.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/MultiScaleModels.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/vfi59h7s6bva5x0m?svg=true)](https://ci.appveyor.com/project/ChrisRackauckas/multiscalemodels-jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaDiffEq/MultiScaleModels.jl/badge.svg)](https://coveralls.io/github/JuliaDiffEq/MultiScaleModels.jl)
[![codecov](https://codecov.io/gh/JuliaDiffEq/MultiScaleModels.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaDiffEq/MultiScaleModels.jl)

MultiScaleModels.jl allows you to easily build multiple scale models which are
fully compatible with the DifferentialEquations.jl ecosystem. These models utilize
a tree structure to describe phenomena of multiple scales, but the interface allows
you to describe equations on different levels, using aggregations from lower
levels to describe complex systems. Their structure allows for complex and dynamic
models to be developed without losing much performance.

## Defining A MultiScaleModel: The Interface

The required interface is as follows. Leaf types must extend MultiScaleModelLeaf, the
highest level of the model or the head extends MultiScaleModelHead, and all
intermediate types extend AbstractMultiScaleModel. The leaf has an array `x::Vector{B}`.
Each type above then contains three fields:

- `x::Vector{T}`
- `y::Vector{B}`
- `end_idxs::Vector{Int}``

`B` is the `BottomType`, which has to be the same as the eltype for the array
in the leaf types. `T` is another `AbstractMultiScaleModel`. Thus at each level,
an` AbstractMultiScaleModel` contains some information of its own (`y`), the
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

## Example

The usage is best described by an example. Here we build a hierarchy where
Embryos contain Tissues which contain Populations which contain Cells, and the
cells contain proteins whose concentrations are modeled as simply a vector
of numbers (it can be anything linearly indexable).

```julia
immutable Cell{B} <: MultiScaleModelLeaf{B}
  x::Vector{B}
end
immutable Population{T<:AbstractMultiScaleModel,B<:Number} <: AbstractMultiScaleModel{B}
  x::Vector{T}
  y::Vector{B}
  end_idxs::Vector{Int}
end
immutable Tissue{T<:AbstractMultiScaleModel,B<:Number} <: AbstractMultiScaleModel{B}
  x::Vector{T}
  y::Vector{B}
  end_idxs::Vector{Int}
end
immutable Embryo{T<:AbstractMultiScaleModel,B<:Number} <: MultiScaleModelHead{B}
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
is `construct(T::AbstractMultiScaleModel,x,y)`, though if `y` is not given it's
taken to be empty.

```julia
p = construct(Population,deepcopy([cell1;cell2])) # Make a Population from cells
cell3 = Cell([3.0;2.0;5.0])
cell4 = Cell([4.0;6])
p2 = construct(Population,deepcopy([cell3;cell4]))
tis = construct(Tissue,deepcopy([p;p2])) # Make a Tissue from Populations
tis2 = construct(Tissue,deepcopy([p2;p]))
em = construct(Embryo,deepcopy([tis;tis2])) # Make an embryo from Tissues
```

The head node then acts as the king. It is designed to have functionality which
mimics a vector in order for usage in DifferentialEquations. So for example

```julia
em[12]
```

returns the "12th cell", counting by Embryo > Tissue > Population > Cell in order
of the vectors. The linear indexing exists for every `AbstractMultiScaleModel`.
These types act as full linear vectors, so standard operations do the sensical
operations:

```julia
em[10] = 4.0 # changes protein concentration 10
em[2,3,1] # Gives the 1st cell in the 3rd population of the second tissue
em[:] # generates a vector of all of the protein concentrations
eachindex(em) # generates an iterator for the indices
```

Continuous models can thus be written at the protein level and will work seamlessly
with DifferentialEquations which will treat it like a vector of protein concentrations.
Using the iterators, note that we can get each cell population by looping through
2 levels below the top, so

```julia
for cell in level_iter(em,2)
  # Do something with the cells!
end
```

To apply a function cell-by-cell, you can write a dispatch `f` on the type for the
level. Using `level_iter_idx`, we can have its changes update some other head node
`d_em` via:

```julia
for cell in level_iter_idx(em,2)
  f(t,cell,@view d_em[y:z])
end
```

Notice that this updates the top vector cell-by-cell via the function `f` without
allocating. This allows one to apply an ODE "cell-wise".

However, the interesting behavior comes from event handling. Since `em` will be the
"vector" for the differential equation, it will be the value passed to the event
handling. MultiScaleModels includes behavior for changing the structure. For example:

```julia
tis3 = construct(Tissue,deepcopy([p;p2]))
add_daughter!(em,tis3) # Adds a new tissue to the embryo
remove_daughter!(em,2,1) # Removes population 1 from tissue 2 of the embryo
```

Combined with event handling, this allows for dynamic structures to be derived from
low level behaviors.

### Extensions

Note that this only showed the most basic MultiScaleModel. These types can be
extended as one pleases. For example, we can change the definition of the cell
to have:

```julia
immutable Cell{B} <: MultiScaleModelLeaf{B}
  x::Vector{B}
  celltype::Symbol
end
```

Then we'd construct cells with `cell3 = Cell([3.0;2.0;5.0],:BCell)`, and can
give it a cell type. This information is part of the call, so

```julia
for cell in level_iter_idx(em,2)
  f(t,cell,@view d_em[y:z])
end
```

can allow one to check the `cell.celltype` in `f` an apply a different ODE depending
on the cell type. You can add fields however you want, so you can use them
to name cells and track lineages.

Showing the use of `y`, you just pass it to the constructor. Let's pass it an array
of 3 values:

```julia
tis = construct(Tissue,deepcopy([p;p2]),[0.0;0.0;0.0])
```

We can selectively apply some function on these `y` values via:

```julia
for tis in level_iter_idx(em,1)
  f(t,tis,@view d_em[y:z])
end
```

and mutatate `tis.y` in `f`. For example, we could have

```julia
function f(t,tis::Tissue,du)
  tis.y += randn(3)
end
```

applies normal random numbers to the three values. We could use this to add to the
model the fact that `tis.y[1:3]` are the tissue's position, and `f` would then be
adding Brownian motion.

Of course, you can keep going and kind of do whatever you want. The power is yours!
