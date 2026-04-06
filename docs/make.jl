using NSTools
using Documenter

DocMeta.setdocmeta!(NSTools, :DocTestSetup, :(using NSTools); recursive=true)

makedocs(;
    modules=[NSTools],
    authors="T. Hashimoto",
    sitename="NSTools.jl",
    format=Documenter.HTML(;
        canonical="https://tkrhsmt.github.io/NSTools.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/tkrhsmt/NSTools.jl",
    devbranch="main",
)
