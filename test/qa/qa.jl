using SciMLTesting, MultiScaleArrays

# Non-public names from Base, Base.Broadcast, ForwardDiff, and RecursiveArrayTools
# are required for broadcast style hooks, generated `similar`, and ODE cache resizing.
run_qa(
    MultiScaleArrays;
    ei_kwargs = (;
        all_qualified_accesses_are_public = (;
            ignore = (
                :AbstractArrayStyle,
                :DefaultArrayStyle,
                :DerivativeConfig,
                :Dual,
                :Tag,
                :_broadcast_getindex_eltype,
                :typename,
            ),
        ),
    ),
)
