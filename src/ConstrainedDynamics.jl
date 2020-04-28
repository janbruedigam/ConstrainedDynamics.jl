module ConstrainedDynamics

using LinearAlgebra
using StaticArrays
using StaticArrays: SUnitRange
using Rotations

using CoordinateTransformations
using GeometryTypes: GeometryTypes, Vec, Point, GLUVMesh
using Blink
using Colors: RGBA, RGB
using FileIO, MeshIO
using MeshCat
using LightXML

using Plots

export Box,
    Cylinder,
    Sphere,
    Mesh,

    Quaternion,
    Origin,
    Body,
    EqualityConstraint,
    InequalityConstraint,
    Mechanism,
    PID,

    OriginConnection,
    Prismatic,
    Spherical,
    Cylindrical,
    Revolute,
    Planar,
    Fixed,
    FixedOrientation,
    CylindricalFree,

    Impact,
    Friction,

    setPosition!,
    setVelocity!,
    setForce!,
    simulate!,
    plotθ,
    plotλ,
    visualize!,

    disassemble,
    getcomponent,
    getbody,
    geteqconstraint,
    getineqconstraint,
    minimalCoordinates,

    RotX,
    RotY,
    RotZ,
    RGBA


include(joinpath("util", "util.jl"))
include(joinpath("util", "customdict.jl"))
include(joinpath("util", "quaternion.jl"))
include(joinpath("util", "shapes.jl"))
include(joinpath("components", "component.jl"))

include(joinpath("components", "body.jl"))
include(joinpath("joints", "joint.jl"))
include(joinpath("bounds", "bound.jl"))

include(joinpath("joints", "translational.jl"))
include(joinpath("joints", "translational0.jl"))
include(joinpath("joints", "translational1.jl"))
include(joinpath("joints", "translational2.jl"))
include(joinpath("joints", "translational3.jl"))
include(joinpath("joints", "rotational.jl"))
include(joinpath("joints", "rotational0.jl"))
include(joinpath("joints", "rotational1.jl"))
include(joinpath("joints", "rotational2.jl"))
include(joinpath("joints", "rotational3.jl"))
include(joinpath("joints", "prototypes.jl"))

include(joinpath("bounds", "impact.jl"))
include(joinpath("bounds", "friction.jl"))

include(joinpath("components", "equalityconstraint.jl"))
include(joinpath("components", "inequalityconstraint.jl"))

include(joinpath("util", "graph.jl"))
include(joinpath("util", "storage.jl"))

include(joinpath("control", "controller.jl"))

include(joinpath("solver", "sparseldu.jl"))
include(joinpath("components", "mechanism.jl"))
include(joinpath("components", "mechanism_functions.jl"))
include(joinpath("components", "initialize.jl"))
include(joinpath("solver", "solverfunctions.jl"))

include(joinpath("util", "urdf.jl"))

include(joinpath("solver", "newton.jl"))

include(joinpath("util", "visualize.jl"))

end
