module ConstrainedDynamics

using LinearAlgebra
using StaticArrays
using ForwardDiff 
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

    Floating,
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

include(joinpath("optional_components", "shapes.jl"))
include(joinpath("optional_components", "storage.jl"))

include(joinpath("main_components", "component.jl"))
include(joinpath("main_components", "state.jl"))
include(joinpath("main_components", "body.jl"))
include(joinpath("main_components", "equalityconstraint.jl"))
include(joinpath("main_components", "inequalityconstraint.jl"))
include(joinpath("main_components", "graph.jl"))
include(joinpath("main_components", "controller.jl"))
include(joinpath("main_components", "sparseldu.jl"))
include(joinpath("main_components", "mechanism_struct.jl"))
include(joinpath("main_components", "mechanism_functions.jl"))

include(joinpath("joints", "abstract_joint.jl"))

include(joinpath("joints", "joint.jl"))
include(joinpath("joints", "translational.jl"))
include(joinpath("joints", "rotational.jl"))
include(joinpath("joints", "genericjoint.jl"))
include(joinpath("joints", "prototypes.jl"))

include(joinpath("bounds", "bound.jl"))
include(joinpath("bounds", "contact.jl"))
include(joinpath("bounds", "impact.jl"))
include(joinpath("bounds", "friction.jl"))

include(joinpath("solver", "solverfunctions.jl"))
include(joinpath("solver", "system.jl"))
include(joinpath("solver", "initconstraints.jl"))
include(joinpath("solver", "newton.jl"))
include(joinpath("solver", "linesearch.jl"))

include(joinpath("discretization", "SymplecticEuler.jl"))
# include(joinpath("discretization", "ImplicitTrapezoid.jl"))

include(joinpath("ui", "mechanism_ui.jl"))
include(joinpath("ui", "simulate.jl"))
include(joinpath("ui", "initialize.jl"))
include(joinpath("ui", "urdf.jl"))


end
