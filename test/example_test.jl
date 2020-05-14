using BenchmarkTools

files = [
    "atlas"
    "chain_in_chain"
    "dice_nofriction"
    "dice_tiltedplane"
    "dice"
    "disconnected_bodies"
    "doublependulum_3d"
    "inverted_pyramid_plane"
    "joint_force"
    "joint_torque"
    "n_fourbars"
    "n_pendulum"
    "pendulum_forced"
    "pendulum"
    "planar_example"
    "scissor_lift"
    "slider_crank"
    "slider_crank3d"
    "urdf_doublependulum"
    "wheel"
]

for file in files
    include("examples/"*file*".jl")
    if file=="joint_force" || file=="pendulum_forced"
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


# include("examples/atlas.jl")
# simulate!(mech, save = true)
# n = length(mech.storage.x)
# for i=1:n
#     @test !any(isnan.(mech.storage.x[i]))
#     @test !any(isnan.(mech.storage.q[i]))
# end

# include("examples/chain_in_chain.jl")
# simulate!(mech, save = true)
# n = length(mech.storage.x)
# for i=1:n
#     @test !any(isnan.(mech.storage.x[i]))
#     @test !any(isnan.(mech.storage.q[i]))
# end

# include("examples/dice_nofriction.jl")
# simulate!(mech, save = true)
# n = length(mech.storage.x)
# for i=1:n
#     @test !any(isnan.(mech.storage.x[i]))
#     @test !any(isnan.(mech.storage.q[i]))
# end

# include("examples/dice_tiltedplane.jl")
# simulate!(mech, save = true)
# n = length(mech.storage.x)
# for i=1:n
#     @test !any(isnan.(mech.storage.x[i]))
#     @test !any(isnan.(mech.storage.q[i]))
# end

# include("examples/dice.jl")
# simulate!(mech, save = true)
# n = length(mech.storage.x)
# for i=1:n
#     @test !any(isnan.(mech.storage.x[i]))
#     @test !any(isnan.(mech.storage.q[i]))
# end

# include("examples/disconnected_bodies.jl")
# simulate!(mech, save = true)
# n = length(mech.storage.x)
# for i=1:n
#     @test !any(isnan.(mech.storage.x[i]))
#     @test !any(isnan.(mech.storage.q[i]))
# end

# include("examples/doublependulum_3d.jl")
# simulate!(mech, save = true)
# n = length(mech.storage.x)
# for i=1:n
#     @test !any(isnan.(mech.storage.x[i]))
#     @test !any(isnan.(mech.storage.q[i]))
# end

# include("examples/inverted_pyramid_plane.jl")
# simulate!(mech, save = true)
# n = length(mech.storage.x)
# for i=1:n
#     @test !any(isnan.(mech.storage.x[i]))
#     @test !any(isnan.(mech.storage.q[i]))
# end

# include("examples/joint_force.jl")
# simulate!(mech, save = true)
# n = length(mech.storage.x)
# for i=1:n
#     @test !any(isnan.(mech.storage.x[i]))
#     @test !any(isnan.(mech.storage.q[i]))
# end

# include("examples/joint_torque.jl")
# simulate!(mech, save = true)
# n = length(mech.storage.x)
# for i=1:n
#     @test !any(isnan.(mech.storage.x[i]))
#     @test !any(isnan.(mech.storage.q[i]))
# end

# include("examples/n_fourbars.jl")
# simulate!(mech, save = true)
# n = length(mech.storage.x)
# for i=1:n
#     @test !any(isnan.(mech.storage.x[i]))
#     @test !any(isnan.(mech.storage.q[i]))
# end

# include("examples/n_pendulum.jl")
# simulate!(mech, save = true)
# n = length(mech.storage.x)
# for i=1:n
#     @test !any(isnan.(mech.storage.x[i]))
#     @test !any(isnan.(mech.storage.q[i]))
# end

# include("examples/pendulum_forced.jl")
# simulate!(mech, save = true)
# n = length(mech.storage.x)
# for i=1:n
#     @test !any(isnan.(mech.storage.x[i]))
#     @test !any(isnan.(mech.storage.q[i]))
# end

# include("examples/pendulum.jl")
# simulate!(mech, save = true)
# n = length(mech.storage.x)
# for i=1:n
#     @test !any(isnan.(mech.storage.x[i]))
#     @test !any(isnan.(mech.storage.q[i]))
# end

# include("examples/planar_example.jl")
# simulate!(mech, save = true)
# n = length(mech.storage.x)
# for i=1:n
#     @test !any(isnan.(mech.storage.x[i]))
#     @test !any(isnan.(mech.storage.q[i]))
# end

# include("examples/scissor_lift.jl")
# simulate!(mech, save = true)
# n = length(mech.storage.x)
# for i=1:n
#     @test !any(isnan.(mech.storage.x[i]))
#     @test !any(isnan.(mech.storage.q[i]))
# end

# include("examples/slider_crank.jl")
# simulate!(mech, save = true)
# n = length(mech.storage.x)
# for i=1:n
#     @test !any(isnan.(mech.storage.x[i]))
#     @test !any(isnan.(mech.storage.q[i]))
# end

# include("examples/slider_crank3d.jl")
# simulate!(mech, save = true)
# n = length(mech.storage.x)
# for i=1:n
#     @test !any(isnan.(mech.storage.x[i]))
#     @test !any(isnan.(mech.storage.q[i]))
# end

# include("examples/urdf_doublependulum.jl")
# simulate!(mech, save = true)
# n = length(mech.storage.x)
# for i=1:n
#     @test !any(isnan.(mech.storage.x[i]))
#     @test !any(isnan.(mech.storage.q[i]))
# end

# include("examples/wheel.jl")
# simulate!(mech, save = true)
# n = length(mech.storage.x)
# for i=1:n
#     @test !any(isnan.(mech.storage.x[i]))
#     @test !any(isnan.(mech.storage.q[i]))
# end

