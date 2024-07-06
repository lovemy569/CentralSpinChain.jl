using CentralSpinChain
using Documenter

DocMeta.setdocmeta!(CentralSpinChain, :DocTestSetup, :(using CentralSpinChain); recursive=true)

makedocs(;
    modules=[CentralSpinChain],
    authors="lovemy569",
    sitename="CentralSpinChain.jl",
    format=Documenter.HTML(;
        canonical="https://lovemy569.github.io/CentralSpinChain.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/lovemy569/CentralSpinChain.jl",
    devbranch="main",
)
