module FullCordDynamics

using LinearAlgebra, StaticArrays
using Rotations

using CoordinateTransformations
using GeometryTypes: # Define geometric shapes
    GeometryPrimitive, GeometryTypes, HyperRectangle, Vec, Point, Rectangle, Cylinder,
    HomogenousMesh, SignedDistanceField, HyperSphere, GLUVMesh, Pyramid
using Blink
using Colors: RGBA, RGB # Handle RGB colors
using FileIO, MeshIO # Load meshes in MeshCat
using MeshCat # Visualize 3D animations

using Plots

export
    Box,

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

    setInit!,
    simulate!,
    plotTraj,
    visualize


include(joinpath("util", "quaternion.jl"))
include(joinpath("util", "shapes.jl"))
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
include(joinpath("util", "storage.jl"))

include(joinpath("solver", "sparseldu2.jl"))
include(joinpath("components", "robot.jl"))
include(joinpath("solver", "sparseldu.jl"))

include(joinpath("solver", "newton.jl"))

include(joinpath("util", "visualize.jl"))

end
