using Documenter, MultiScaleArrays

include("pages.jl")

makedocs(; sitename = "MultiScaleArrays.jl",
         authors = "Chris Rackauckas",
         modules = [MultiScaleArrays],
         clean = true, doctest = false,
         format = Documenter.HTML(; analytics = "UA-90474609-3",
                                  assets = ["assets/favicon.ico"],
                                  canonical = "https://multiscalearrays.sciml.ai/stable/"),
         pages = pages)

deploydocs(; repo = "github.com/SciML/MultiScaleArrays.jl.git",
           push_preview = true)
