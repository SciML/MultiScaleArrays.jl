# API

```@docs
AbstractMultiScaleArray
print_human_readable
```

### Extensions

Note that this only showed the most basic MultiScaleArray. These types can be
extended as one pleases. For example, we can change the definition of the cell
to have:

```julia
struct Cell{B} <: AbstractMultiScaleArrayLeaf{B}
    values::Vector{B}
    celltype::Symbol
end
```

Note that the ordering of the fields matters here: the extra fields must come
after the standard fields (so for a leaf it comes after `values`, for a standard
multiscale array it would come after `nodes,values,end_idxs`).
Then we'd construct cells with
`cell3 = Cell([3.0; 2.0; 5.0], :BCell)`, and can give it a cell type.
This information is part of the call, so

```julia
for (cell, y, z) in level_iter_idx(embryo, 2)
    f(t, cell, @view embryo[y:z])
end
```

can allow one to check the `cell.celltype` in `f` an apply a different ODE depending
on the cell type. You can add fields however you want, so you can use them
to name cells and track lineages.

Showing the use of `values`, you just pass it to the constructor. Let's pass it an array
of 3 values:

```julia
tissue = construct(Tissue, deepcopy([population; population2]), [0.0; 0.0; 0.0])
```

We can selectively apply some function on these `values` via:

```julia
for (tissue, y, z) in level_iter_idx(embryo, 1)
    f(t, tissue, @view embryo[y:z])
end
```

and mutate `tis.values` in `f`. For example, we could have

```julia
function f(du, tissue::Tissue, p, t)
    du .+= randn(3)
end
```

applies normal random numbers to the three values. We could use this to add to the
model the fact that `tissue.values[1:3]` are the tissue's position, and `f` would then be
adding Brownian motion.

Of course, you can keep going and kind of do whatever you want. The power is yours!
