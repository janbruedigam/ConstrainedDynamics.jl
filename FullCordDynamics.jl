module FullCordDynamics

# using TimerOutputs
# const to = TimerOutput()

using LinearAlgebra, StaticArrays
using StaticArrays: SUnitRange
using Rotations

using CoordinateTransformations
using GeometryTypes:
    GeometryPrimitive, GeometryTypes, Vec, Point, Rectangle,
    HomogenousMesh, SignedDistanceField, HyperSphere, GLUVMesh, Pyramid
using Blink
using Colors: RGBA, RGB
using FileIO, MeshIO
using MeshCat

using Plots

export
    Box,
    Cylinder,

    Quaternion,
    Origin,
    Link,
    Constraint,
    Robot,

    OriginConnection,
    Axis,
    Line,
    Socket,
    SocketYZ,
    MatchedOrientation,

    setInit!,
    simulate!,
    plotθ,
    plotλ,
    visualize,

    simulate_energy!,
    simulate_drift!,
    simulate_reset!,
    simulate_steptol!


include(joinpath("util", "util.jl"))
include(joinpath("util", "customdict.jl"))
include(joinpath("util", "quaternion.jl"))
include(joinpath("util", "shapes.jl"))
include(joinpath("components", "component.jl"))
include(joinpath("joints", "joint.jl"))
include(joinpath("components", "link.jl"))
include(joinpath("components", "constraint.jl"))

include(joinpath("joints", "matchedorientation.jl"))
include(joinpath("joints", "originconnection.jl"))
include(joinpath("joints", "socket.jl"))
include(joinpath("joints", "socketyz.jl"))
include(joinpath("joints", "axis.jl"))
include(joinpath("joints", "line.jl"))

include(joinpath("util", "graph.jl"))
include(joinpath("util", "storage.jl"))

include(joinpath("solver", "sparseldu.jl"))
include(joinpath("components", "robot.jl"))
include(joinpath("solver", "solverfunctions.jl"))

include(joinpath("solver", "newton.jl"))
include(joinpath("paper_experiments", "experiment_methods.jl"))

include(joinpath("util", "visualize.jl"))

end
