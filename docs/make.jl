using QuantumStateTransfer
using Documenter

DocMeta.setdocmeta!(QuantumStateTransfer, :DocTestSetup, :(using QuantumStateTransfer); recursive=true)

makedocs(;
    modules=[QuantumStateTransfer],
    authors="Luis M. B. Varona <lbvarona@mta.ca>",
    sitename="QuantumStateTransfer.jl",
    format=Documenter.HTML(;
        canonical="https://Luis-Varona.github.io/QuantumStateTransfer.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Luis-Varona/QuantumStateTransfer.jl",
    devbranch="main",
)
