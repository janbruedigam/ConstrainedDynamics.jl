module ConstrainedDynamics

using LinearAlgebra
using StaticArrays
using ForwardDiff 
using StaticArrays: SUnitRange
using Quaternions
using Rotations
using Rotations: RotationError, params, lmult, rmult, tmat, vmat, hmat, skew
using Colors: RGBA, RGB
using LightXML
using GraphBasedSystems
using GraphBasedSystems: Entry

using DocStringExtensions


export Origin,
    Body,
    EqualityConstraint,
    InequalityConstraint,
    Friction,
    Mechanism,
    Controller,
    Storage,
    Quaternion,

    Box,
    Cylinder,
    Sphere,
    Pyramid,
    Mesh,

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
    Quaternion,

    setPosition!,
    setVelocity!,
    setForce!,
    addForce!,
    getid,
    getcomponent,
    getbody,
    geteqconstraint,
    getfriction,
    getineqconstraint,
    simulate!,
    initializeConstraints!,
    disassemble,
    minimalCoordinates,
    minimalVelocities,
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
include(joinpath("main_components", "abstractconstraint.jl"))
include(joinpath("main_components", "equalityconstraint.jl"))
include(joinpath("main_components", "inequalityconstraint.jl"))
include(joinpath("main_components", "friction.jl"))
include(joinpath("main_components", "controller.jl"))
include(joinpath("main_components", "mechanism_struct.jl"))
include(joinpath("main_components", "system.jl"))
include(joinpath("main_components", "mechanism_functions.jl"))

include(joinpath("joints", "abstract_joint.jl"))

include(joinpath("bounds", "bound.jl"))
include(joinpath("bounds", "impact.jl"))
include(joinpath("bounds", "friction_bounds.jl"))

include(joinpath("joints", "joint.jl"))
include(joinpath("joints", "translational.jl"))
include(joinpath("joints", "rotational.jl"))
include(joinpath("joints", "genericjoint.jl"))
include(joinpath("joints", "prototypes.jl"))
# include(joinpath("joints", "friction.jl"))

include(joinpath("solver", "solverfunctions.jl"))
include(joinpath("solver", "initconstraints.jl"))
include(joinpath("solver", "newton.jl"))
include(joinpath("solver", "linesearch.jl"))
include(joinpath("optional_components", "linearization.jl"))

include(joinpath("discretization", "Linear.jl"))
# include(joinpath("discretization", "Quadratic.jl"))

include(joinpath("ui", "mechanism_ui.jl"))
include(joinpath("ui", "simulate.jl"))
include(joinpath("ui", "initialize.jl"))
include(joinpath("ui", "urdf.jl"))


end
