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
    FillIn,
    JointNode,
    Link,
    Origin,
    Constraint,
    Joint,
    Robot,

    Axis,
    Socket,
    SocketYZ,
    # FixedOrientation,
    # FixedPosition,
    Combined2,
    Combined3,

    box,
    initialPosition,
    setInit!,
    sim!,
    trajSFunc,
    plotTraj,

    parse_urdf


include(joinpath("util", "quaternion.jl"))
include(joinpath("components", "node.jl"))
include(joinpath("components", "fillin.jl"))
include(joinpath("components", "link.jl"))

include(joinpath("joints", "joint.jl"))
# include(joinpath("joints", "fixedposition.jl"))
# include(joinpath("joints", "fixedorientation.jl"))
include(joinpath("joints", "socket.jl"))
include(joinpath("joints", "socketyz.jl"))
include(joinpath("joints", "axis.jl"))
include(joinpath("components", "constraint.jl"))
include(joinpath("components", "combined2.jl"))
include(joinpath("components", "combined3.jl"))

include(joinpath("util", "util.jl"))
include(joinpath("util", "graph.jl"))
include(joinpath("util", "shapes.jl"))
include(joinpath("util", "storage.jl"))

include("robot.jl")
# include("sparseldu.jl")

include("newton.jl")

end
