using Documenter
using ConstrainedDynamics

makedocs(;
    modules = [ConstrainedDynamics],
    format = Documenter.HTML(
        canonical = "https://janbruedigam.github.io/ConstrainedDynamics.jl/stable/",
        assets = ["assets/favicon.ico"],
    ),
    pages = [
        "Home" => "index.md",
        "Getting started" => [
            "From URDF" => "gettingstarted/urdf.md",
            "From Code" => "gettingstarted/code.md"
        ],
        "Library" => [
            "Body" => "library/body.md",
            "Shape" => "library/shape.md",
            "Constraint" => "library/equalityconstraint.md",
            "Mechanism" => "library/mechanism.md",
            "Interface" => "library/interface.md",
            "Simulation" => "library/simulation.md",
            "State" => "library/state.md"
        ]
    ],
    sitename = "ConstrainedDynamics.jl"
)

deploydocs(; repo = "github.com/janbruedigam/ConstrainedDynamics.jl")