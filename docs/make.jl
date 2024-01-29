using Zerilli
using Documenter

DocMeta.setdocmeta!(Zerilli, :DocTestSetup, :(using Zerilli); recursive=true)

makedocs(;
    modules=[Zerilli],
    authors="Stamatis Vretinaris",
    sitename="Zerilli.jl",
    format=Documenter.HTML(;
        canonical="https://svretina.github.io/Zerilli.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/svretina/Zerilli.jl",
    devbranch="master",
)
