using Documenter, MultiScaleArrays

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

include("pages.jl")

makedocs(sitename = "MultiScaleArrays.jl",
         authors = "Chris Rackauckas",
         modules = [MultiScaleArrays],
         clean = true, doctest = false,
         format = Documenter.HTML(analytics = "UA-90474609-3",
                                  assets = ["assets/favicon.ico"],
                                  canonical = "https://docs.sciml.ai/MultiScaleArrays/stable/"),
         pages = pages)

deploydocs(repo = "github.com/SciML/MultiScaleArrays.jl.git";
           push_preview = true)
