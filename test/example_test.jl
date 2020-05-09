using BenchmarkTools

include("examples/atlas.jl")
simulate!(mech)
@test true

include("examples/chain_in_chain.jl")
simulate!(mech)
@test true

include("examples/dice_nofriction.jl")
simulate!(mech)
@test true

include("examples/dice_tiltedplane.jl")
simulate!(mech)
@test true

include("examples/dice.jl")
simulate!(mech)
@test true

include("examples/disconnected_bodies.jl")
simulate!(mech)
@test true

include("examples/doublependulum_3d.jl")
simulate!(mech)
@test true

include("examples/inverted_pyramid_plane.jl")
simulate!(mech)
@test true

include("examples/joint_force.jl")
simulate!(mech)
@test true

include("examples/joint_torque.jl")
simulate!(mech)
@test true

include("examples/n_fourbars.jl")
simulate!(mech)
@test true

include("examples/n_pendulum.jl")
simulate!(mech)
@test true

include("examples/pendulum_forced.jl")
simulate!(mech)
@test true

include("examples/pendulum.jl")
simulate!(mech)
@test true

include("examples/planar_example.jl")
simulate!(mech)
@test true

include("examples/scissor_lift.jl")
simulate!(mech)
@test true

include("examples/slider_crank.jl")
simulate!(mech)
@test true

include("examples/slider_crank3d.jl")
simulate!(mech)
@test true

include("examples/urdf_doublependulum.jl")
simulate!(mech)
@test true

include("examples/wheel.jl")
simulate!(mech)
@test true

