using SciMLTesting, MultiScaleArrays

# ExplicitImports only sees an extension module once its trigger package is loaded, so QA
# environments normally add the weakdeps here to bring the extensions into scope.
#
# MultiScaleArraysSparseDiffToolsExt stays unscanned: SparseDiffTools 2 caps SciMLOperators
# at 0.4, while SciMLBase 3.28.1+ and OrdinaryDiffEqDifferentiation 3 (both hard deps)
# require SciMLOperators 1.3+, so SparseDiffTools cannot be installed alongside
# MultiScaleArrays at all and the extension can never load.

run_qa(MultiScaleArrays)
