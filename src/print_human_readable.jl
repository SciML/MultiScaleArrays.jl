############### print_human_readable(X::AbstractMultiScaleArray
function add_separator!(first_element, toprint, level)
    if first_element
        toprint[level][end] = "+" * toprint[level][end]
    end
    return false
end

function take_first_char(s::Union{String, SubString{String}}, N::Int64)
    s[1:min(length(s), N)]
end

function toprint_AbstractMultiScaleArray!(toprint, X::AbstractMultiScaleArray, fields,
        levelmax, n_char_per_name, level = 1)
    if length(toprint) < level
        push!(toprint, Vector{String}())
    end
    first_element = true
    if level < levelmax
        for x in X.nodes
            if fieldname(typeof(x), 1) == :nodes
                push!(toprint[level],
                    "|" *
                    take_first_char(split(string(typeof(x)), "{")[1], n_char_per_name))
                first_element = add_separator!(first_element, toprint, level)
                toprint_AbstractMultiScaleArray!(toprint, x, fields, levelmax,
                    n_char_per_name, level + 1)
            else
                if fields == nothing
                    push!(toprint[level],
                        take_first_char(split(string(typeof(x)), "{")[1],
                            n_char_per_name))
                    first_element = add_separator!(first_element, toprint, level)
                else
                    if !(x isa AbstractMultiScaleArrayLeaf)
                        push!(toprint[level],
                            take_first_char(split(string(typeof(x)), "{")[1],
                                n_char_per_name))
                        first_element = add_separator!(first_element, toprint, level)
                    else
                        push!(toprint[level],
                            join(
                                [take_first_char(string(field), n_char_per_name) * ": " *
                                 string(getfield(x, Symbol(field)))
                                 for field in collect(fields)],
                                ", "))
                        first_element = add_separator!(first_element, toprint, level)
                    end
                end
            end
        end
    end
end

"""
```julia
print_human_readable(embryo)
# +|Tissue;                                                 |Tissue
#  +|Popula;           |Popula;           |Popula;          +|Popula;           |Popula;           |Popula
#   +Cell; Cell; Cell; +Cell; Cell; Cell; +Cell; Cell; Cell; +Cell; Cell; Cell; +Cell; Cell; Cell; +Cell; Cell; Cell

print_human_readable(embryo; NcharPerName = 2)
# +|Ti;                                   |Ti
#  +|Po;         |Po;         |Po;        +|Po;         |Po;         |Po
#   +Ce; Ce; Ce; +Ce; Ce; Ce; +Ce; Ce; Ce; +Ce; Ce; Ce; +Ce; Ce; Ce; +Ce; Ce; Ce
```

Here, if the 'AbstractMultiScaleArrayLeaf's contain several fields, you can specify them with fields = [field1,field2,...]

```julia
print_human_readable(embryo; NcharPerName = 2, fields = [:values])
# +|Ti;                                                                                                                                                                             |Ti
#  +|Po;                                                       |Po;                                                       |Po;                                                      +|Po;                                                       |Po;                                                       |Po
#   +va: [1.0, 2.0, 3.0]; va: [3.0, 2.0, 5.0]; va: [4.0, 6.0]; +va: [1.0, 2.0, 3.0]; va: [3.0, 2.0, 5.0]; va: [4.0, 6.0]; +va: [1.0, 2.0, 3.0]; va: [3.0, 2.0, 5.0]; va: [4.0, 6.0]; +va: [1.0, 2.0, 3.0]; va: [3.0, 2.0, 5.0]; va: [4.0, 6.0]; +va: [1.0, 2.0, 3.0]; va: [3.0, 2.0, 5.0]; va: [4.0, 6.0]; +va: [1.0, 2.0, 3.0]; va: [3.0, 2.0, 5.0]; va: [4.0, 6.0]
```

if your screen is small, then print a sub-part of the AbstractMultiScaleArray:

```julia
print_human_readable(embryo.nodes[1].nodes[1]; fields = [:values])
# +values: [1.0, 2.0, 3.0]; values: [3.0, 2.0, 5.0]; values: [4.0, 6.0]
```
"""
function print_human_readable(X::AbstractMultiScaleArray; n_char_per_name = 6,
        fields = nothing, levelmax = Inf, n_item_max_per_levels = Inf) # fields = nothing OR [field1,field2,...]
    #     if X isa AbstractMultiScaleArrayLeaf
    #         println(X)
    #     else
    toprint = Vector{Vector{String}}([[]]) # one vector per row, one string per element
    toprint_AbstractMultiScaleArray!(toprint, X, fields, levelmax, n_char_per_name)
    toprint = map(x -> join(x, "; "), toprint)
    #     toprint[1] = join(split(toprint[1]," "),"; ")
    Done = false
    while !Done
        Done = true # temporary
        for level1 in length(toprint):-1:2
            toprint0 = split(toprint[level1 - 1], "|")
            toprint1 = split(toprint[level1], "+")
            toprint0_, toprint1_ = "", ""
            for (x0, x1) in zip(toprint0, toprint1)
                if length(x0) > length(x1)
                    x1 = join([x1, repeat(" ", length(x0) - length(x1))], "")
                    Done = false
                elseif length(x0) < length(x1)
                    x0 = join([x0, repeat(" ", length(x1) - length(x0))], "")
                    Done = false
                end
                toprint0_ = join([toprint0_, x0], "|")
                toprint1_ = join([toprint1_, x1], "+")
            end
            toprint0_ = toprint0_[2:end]
            toprint1_ = toprint1_[2:end]
            toprint[level1 - 1] = toprint0_
            toprint[level1] = toprint1_
        end
    end
    # handle "ula; +         |Popu" cases
    #              ^^^^^^^^^^^
    for level in length(toprint):-1:1
        while occursin("+ ", toprint[level])
            Move = findfirst("+ ", toprint[level])[1]
            To = findnext(" |", toprint[level], Move[end])[1]
            toprint[level] = join([
                toprint[level][1:(Move - 1)],
                toprint[level][(Move + 1):To],
                toprint[level][Move],
                toprint[level][(To + 1):end]
            ])
            for level_ in 1:(level - 1)
                toprint[level_] = join([
                    toprint[level_][1:(Move - 1)],
                    toprint[level_][(To + 1):(To + To - Move)],
                    toprint[level_][Move:To],
                    toprint[level_][(To + To - Move + 1):end]
                ])
            end
        end
    end
    println(join(toprint, "\n"))
    #     end
end
