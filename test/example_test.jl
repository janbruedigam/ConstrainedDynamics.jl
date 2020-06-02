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
    "doublependulum_urdf"
    "football"
    "inverted_pyramid_plane"
    "joint_force"
    "joint_torque"
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

for file in files
    # println(file)
    include("examples/"*file*".jl")
    if file=="joint_force" || file=="pendulum_forced" || file=="nutation" || file=="football"
        storage = simulate!(mech, 10., control!, record = true)
    else
        storage = simulate!(mech, 10., record = true)
    end
    n = length(storage.x)
    for i=1:n
        @test !any(isnan.(storage.x[i]))
        @test !any(isnan.(storage.q[i]))
        @test !any(isnan.(storage.v[i]))
        @test !any(isnan.(storage.ω[i]))
    end
end