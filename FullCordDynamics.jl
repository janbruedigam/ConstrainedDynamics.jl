module FullCordDynamics

using LinearAlgebra
using StaticArrays
using Rotations
using Plots

export
    Quaternion,
    Node,
    Link,
    Constraint,
    Robot,

    Axis,
    Socket,
    FixedOrientation,
    FixedPosition,
    Combined,

    box,
    initialPosition,
    setInit!,
    sim!,
    trajSFunc,
    plotTraj


include(joinpath("util", "quaternion.jl"))
include("node.jl")
include("link.jl")

include(joinpath("constraints", "constraint.jl"))
include(joinpath("constraints", "fixedposition.jl"))
include(joinpath("constraints", "fixedorientation.jl"))
include(joinpath("constraints", "socket.jl"))
include(joinpath("constraints", "axis.jl"))
include(joinpath("constraints", "combined.jl"))

include(joinpath("util", "graph.jl"))
include(joinpath("util", "shapes.jl"))

include("robot.jl")
include("sparseldu.jl")

include("newton.jl")
end
