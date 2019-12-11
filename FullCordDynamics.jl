module FullCordDynamics

using LinearAlgebra
using StaticArrays
using Rotations
using Plots
using Base.Threads
using LightXML


export
    Quaternion,
    Node,
    JointNode,
    Link,
    Constraint,
    Joint,
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
    plotTraj,

    parse_urdf


include(joinpath("util", "quaternion.jl"))
include(joinpath("components", "node.jl"))
include(joinpath("components", "link.jl"))

include(joinpath("joints", "joint.jl"))
include(joinpath("joints", "fixedposition.jl"))
include(joinpath("joints", "fixedorientation.jl"))
include(joinpath("joints", "socket.jl"))
include(joinpath("joints", "axis.jl"))
include(joinpath("components", "constraint.jl"))

include(joinpath("util", "graph.jl"))
include(joinpath("util", "shapes.jl"))

include("robot.jl")
include("sparseldu.jl")

include("newton.jl")

include(joinpath("util", "parseurdf.jl"))
end
