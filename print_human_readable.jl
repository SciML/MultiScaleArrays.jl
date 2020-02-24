############### println(X::AbstractMultiScaleArray
function addSeparator!(firstElement,ToPrint,level)
    if firstElement
        ToPrint[level][end] = "+"*ToPrint[level][end]
    end
    return false
end

function Take1stChar(s::Union{String,SubString{String}},N::Int64)     s[1:min(length(s),N)]     end

function ToPrintAbstractMultiScaleArray!(ToPrint,X::AbstractMultiScaleArray,fields,LevelMax,NcharPerName,level = 1)
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
                    if !(typeof(x) <: AbstractMultiScaleArrayLeaf)
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
    println(join(toprint,"\n"))
end
