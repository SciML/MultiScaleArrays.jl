# MultiScaleModels


[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Build Status](https://travis-ci.org/JuliaDiffEq/MultiScaleModels.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/MultiScaleModels.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/vfi59h7s6bva5x0m?svg=true)](https://ci.appveyor.com/project/ChrisRackauckas/multiscalemodels-jl)
[![Coverage Status](https://coveralls.io/repos/ChrisRackauckas/MultiScaleModels.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/ChrisRackauckas/MultiScaleModels.jl?branch=master)
[![codecov.io](http://codecov.io/github/ChrisRackauckas/MultiScaleModels.jl/coverage.svg?branch=master)](http://codecov.io/github/ChrisRackauckas/MultiScaleModels.jl?branch=master)

MultiScaleModels.jl allows you to easily build multiple scale models which are
fully compatible with the DifferentialEquations.jl ecosystem. These models utilize
a tree structure to describe phenomena of multiple scales, but the interface allows
you to describe equations on different levels, using aggregations from lower
levels to describe complex systems. Their structure allows for complex and dynamic
models to be developed without losing much performance.

## The MultiScaleModel Interface

The interface is as follows. Leaf types must extend MultiScaleModelLeaf, the
highest level of the model or the head extends MultiScaleModelHead, and all
intermediate types extend AbstractMultiScaleModel. Each time then contains
two fields: x::Vector{T} and end_idxs::Vector{Int} (note that leaf types do not
have end_idxs). The rest is up to you!
You can add extra fields as you please, and even make the types immutable.

## Example

The usage is best described by an example. Here we build a hierarchy where
Embryos contain Tissues which contain Populations which contain Cells, and the
cells contain proteins whose concentrations are modeled as simply a vector
of numbers (it can be anything linearly indexable).

```julia
immutable Cell{T} <: MultiScaleModelLeaf
  x::Vector{T}
end
immutable Population{T<:AbstractMultiScaleModel} <: AbstractMultiScaleModel
  x::Vector{T}
  end_idxs::Vector{Int}
end
immutable Tissue{T<:AbstractMultiScaleModel} <: AbstractMultiScaleModel
  x::Vector{T}
  end_idxs::Vector{Int}
end
immutable Embryo{T<:AbstractMultiScaleModel} <: MultiScaleModelHead
  x::Vector{T}
  end_idxs::Vector{Int}
end
```

Using this we can directly construct leaf types:

```julia
cell1 = Cell([1.0;2.0;3.0])
cell2 = Cell([4.0;5])
```

and build types higher up in the hierarchy by using the `constuct` method:

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
