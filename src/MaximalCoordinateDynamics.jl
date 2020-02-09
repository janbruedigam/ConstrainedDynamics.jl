module MaximalCoordinateDynamics

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
    Body,
    EqualityConstraint,
    InequalityConstraint,
    Mechanism,

    OriginConnection,
    Prismatic,
    Spherical,
    Cylindrical,
    Revolute,
    Planar,
    Fixed,
    FixedOrientation,
    CylindricalFree,

    setInit!,
    simulate!,
    simulate_ip!,
    plotθ,
    plotλ,
    visualize


include(joinpath("util", "util.jl"))
include(joinpath("util", "customdict.jl"))
include(joinpath("util", "quaternion.jl"))
include(joinpath("util", "shapes.jl"))
include(joinpath("components", "component.jl"))
include(joinpath("joints", "joint.jl"))
include(joinpath("components", "body.jl"))
include(joinpath("components", "constraint.jl"))
include(joinpath("contacts", "impact.jl"))

include(joinpath("joints", "translationalrotational6.jl"))
include(joinpath("joints", "translational0.jl"))
include(joinpath("joints", "translational1.jl"))
include(joinpath("joints", "translational2.jl"))
include(joinpath("joints", "rotational0.jl"))
include(joinpath("joints", "rotational1.jl"))
include(joinpath("joints", "rotational2.jl"))
include(joinpath("joints", "prototypes.jl"))

include(joinpath("components", "equalityconstraint.jl"))
include(joinpath("components", "inequalityconstraint.jl"))

include(joinpath("util", "graph.jl"))
include(joinpath("util", "storage.jl"))

include(joinpath("solver", "sparseldu.jl"))
include(joinpath("components", "mechanism.jl"))
include(joinpath("solver", "solverfunctions.jl"))

include(joinpath("solver", "newton.jl"))
include(joinpath("solver", "newton_ip.jl"))

include(joinpath("util", "visualize.jl"))

end
