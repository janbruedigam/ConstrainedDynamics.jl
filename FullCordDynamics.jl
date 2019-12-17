module FullCordDynamics

using LinearAlgebra
using StaticArrays
using Rotations
using Plots

export
    Quaternion,
    Origin,
    Link,
    Constraint,
    Robot,

    Axis,
    Socket,
    SocketYZ,
    # FixedOrientation,
    # FixedPosition,

    box,
    initialPosition,
    setInit!,
    sim!,
    trajSFunc,
    plotTraj


include(joinpath("util", "quaternion.jl"))
include(joinpath("components", "node.jl"))
include(joinpath("joints", "joint.jl"))
include(joinpath("components", "link.jl"))
include(joinpath("components", "constraint.jl"))

# include(joinpath("joints", "fixedposition.jl"))
# include(joinpath("joints", "fixedorientation.jl"))
include(joinpath("joints", "socket.jl"))
include(joinpath("joints", "socketyz.jl"))
include(joinpath("joints", "axis.jl"))

include(joinpath("util", "util.jl"))
include(joinpath("util", "graph.jl"))
include(joinpath("util", "shapes.jl"))
include(joinpath("util", "storage.jl"))

include("sparseldu2.jl")
include("robot.jl")
include("sparseldu.jl")

include("newton.jl")

end
