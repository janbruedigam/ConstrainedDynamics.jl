using BenchmarkTools

files = [
    "atlas"
    "axes"
    "chain_in_chain"
    "dice_nofriction"
    "dice_tiltedplane"
    "dice"
    "disconnected_bodies"
    "doublependulum_3d"
    "doublependulum_disconnection"
    "doublependulum_urdf"
    "elliptic_joint"
    "fourbar_disconnection"
    "inverted_pyramid_plane"
    "joint_force"
    "joint_torque"
    "linear_pendulum"
    "n_fourbars"
    "n_pendulum"
    "nutation"
    "pendulum_forced"
    "pendulum"
    "planar_example"
    "scissor_lift"
    "slider_crank"
    "slider_crank3d"
    "wheel"
]

controlled = [
    "doublependulum_disconnection"
    "joint_force"
    "joint_torque"
    "pendulum_forced"
    "nutation"
    "fourbar_disconnection"
    "slider_crank"
]

for file in files
    include("examples/"*file*".jl")
    if file in controlled
        storage = simulate!(mech, 10., eval(Meta.parse(file*"_control!")), record = true)
    else
        storage = simulate!(mech, 10., record = true)
    end
    n = length(storage.x)
    for i=1:n
        @test !any(ConstrainedDynamics.sisnan.(storage.x[i]))
        @test !any(ConstrainedDynamics.sisnan.(storage.q[i]))
        @test !any(ConstrainedDynamics.sisnan.(storage.v[i]))
        @test !any(ConstrainedDynamics.sisnan.(storage.ω[i]))
    end
end