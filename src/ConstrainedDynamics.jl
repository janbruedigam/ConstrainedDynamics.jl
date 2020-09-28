module ConstrainedDynamics

using LinearAlgebra
using StaticArrays
using StaticArrays: SUnitRange
using Rotations
using Rotations: RotationError, pure_quaternion, params, lmult, rmult, tmat, vmat, hmat, skew

using Colors: RGBA, RGB
using LightXML

export Box,
    Cylinder,
    Sphere,
    Mesh,

    UnitQuaternion,
    Origin,
    Body,
    EqualityConstraint,
    InequalityConstraint,
    Mechanism,
    LinearMechanism,
    Storage,
    Controller,

    OriginConnection,
    Prismatic,
    Spherical,
    Cylindrical,
    Revolute,
    Planar,
    PlanarFree,
    Fixed,
    FixedOrientation,
    CylindricalFree,

    Impact,
    Friction,

    setPosition!,
    setVelocity!,
    setForce!,
    simulate!,
    initializeConstraints!,

    disassemble,
    getid,
    getcomponent,
    getbody,
    geteqconstraint,
    getineqconstraint,
    minimalCoordinates,
    activate!,
    deactivate!,

    linearsystem,

    RotX,
    RotY,
    RotZ,
    RGBA,

    szeros,
    sones,
    srand


include(joinpath("util", "util.jl"))
include(joinpath("util", "custom_static.jl"))
include(joinpath("util", "customdict.jl"))
include(joinpath("util", "quaternion.jl"))
include(joinpath("util", "shapes.jl"))

include(joinpath("components", "component.jl"))
include(joinpath("components", "state.jl"))
include(joinpath("components", "body.jl"))

include(joinpath("joints", "joint.jl"))
include(joinpath("bounds", "bound.jl"))

include(joinpath("joints", "translational.jl"))
include(joinpath("joints", "rotational.jl"))
include(joinpath("joints", "prototypes.jl"))

include(joinpath("bounds", "contact.jl"))
include(joinpath("bounds", "impact.jl"))
include(joinpath("bounds", "friction.jl"))

include(joinpath("components", "equalityconstraint.jl"))
include(joinpath("components", "inequalityconstraint.jl"))

include(joinpath("util", "graph.jl"))
include(joinpath("util", "storage.jl"))

include(joinpath("components", "controller.jl"))

include(joinpath("solver", "sparseldu.jl"))
include(joinpath("components", "mechanism.jl"))
include(joinpath("components", "linear_mechanism.jl"))
include(joinpath("components", "mechanism_functions.jl"))
include(joinpath("components", "mechanism_ui.jl"))
include(joinpath("components", "simulate.jl"))
include(joinpath("components", "initialize.jl"))
include(joinpath("solver", "solverfunctions.jl"))
include(joinpath("solver", "system.jl"))
include(joinpath("solver", "initconstraints.jl"))

include(joinpath("util", "urdf.jl"))

include(joinpath("solver", "newton.jl"))
include(joinpath("solver", "linesearch.jl"))


include(joinpath("discretization", "SymplecticEuler.jl"))
# include(joinpath("discretization", "ImplicitTrapezoid.jl"))
end
