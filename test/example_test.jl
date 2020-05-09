using BenchmarkTools

include("examples/atlas.jl")
simulate!(mech, save = true)
n = length(mech.storage.x)
for i=1:n
    @test !any(isnan.(mech.storage.x[i]))
    @test !any(isnan.(mech.storage.q[i]))
end

include("examples/chain_in_chain.jl")
simulate!(mech, save = true)
n = length(mech.storage.x)
for i=1:n
    @test !any(isnan.(mech.storage.x[i]))
    @test !any(isnan.(mech.storage.q[i]))
end

include("examples/dice_nofriction.jl")
simulate!(mech, save = true)
n = length(mech.storage.x)
for i=1:n
    @test !any(isnan.(mech.storage.x[i]))
    @test !any(isnan.(mech.storage.q[i]))
end

include("examples/dice_tiltedplane.jl")
simulate!(mech, save = true)
n = length(mech.storage.x)
for i=1:n
    @test !any(isnan.(mech.storage.x[i]))
    @test !any(isnan.(mech.storage.q[i]))
end

include("examples/dice.jl")
simulate!(mech, save = true)
n = length(mech.storage.x)
for i=1:n
    @test !any(isnan.(mech.storage.x[i]))
    @test !any(isnan.(mech.storage.q[i]))
end

include("examples/disconnected_bodies.jl")
simulate!(mech, save = true)
n = length(mech.storage.x)
for i=1:n
    @test !any(isnan.(mech.storage.x[i]))
    @test !any(isnan.(mech.storage.q[i]))
end

include("examples/doublependulum_3d.jl")
simulate!(mech, save = true)
n = length(mech.storage.x)
for i=1:n
    @test !any(isnan.(mech.storage.x[i]))
    @test !any(isnan.(mech.storage.q[i]))
end

include("examples/inverted_pyramid_plane.jl")
simulate!(mech, save = true)
n = length(mech.storage.x)
for i=1:n
    @test !any(isnan.(mech.storage.x[i]))
    @test !any(isnan.(mech.storage.q[i]))
end

include("examples/joint_force.jl")
simulate!(mech, save = true)
n = length(mech.storage.x)
for i=1:n
    @test !any(isnan.(mech.storage.x[i]))
    @test !any(isnan.(mech.storage.q[i]))
end

include("examples/joint_torque.jl")
simulate!(mech, save = true)
n = length(mech.storage.x)
for i=1:n
    @test !any(isnan.(mech.storage.x[i]))
    @test !any(isnan.(mech.storage.q[i]))
end

include("examples/n_fourbars.jl")
simulate!(mech, save = true)
n = length(mech.storage.x)
for i=1:n
    @test !any(isnan.(mech.storage.x[i]))
    @test !any(isnan.(mech.storage.q[i]))
end

include("examples/n_pendulum.jl")
simulate!(mech, save = true)
n = length(mech.storage.x)
for i=1:n
    @test !any(isnan.(mech.storage.x[i]))
    @test !any(isnan.(mech.storage.q[i]))
end

include("examples/pendulum_forced.jl")
simulate!(mech, save = true)
n = length(mech.storage.x)
for i=1:n
    @test !any(isnan.(mech.storage.x[i]))
    @test !any(isnan.(mech.storage.q[i]))
end

include("examples/pendulum.jl")
simulate!(mech, save = true)
n = length(mech.storage.x)
for i=1:n
    @test !any(isnan.(mech.storage.x[i]))
    @test !any(isnan.(mech.storage.q[i]))
end

include("examples/planar_example.jl")
simulate!(mech, save = true)
n = length(mech.storage.x)
for i=1:n
    @test !any(isnan.(mech.storage.x[i]))
    @test !any(isnan.(mech.storage.q[i]))
end

include("examples/scissor_lift.jl")
simulate!(mech, save = true)
n = length(mech.storage.x)
for i=1:n
    @test !any(isnan.(mech.storage.x[i]))
    @test !any(isnan.(mech.storage.q[i]))
end

include("examples/slider_crank.jl")
simulate!(mech, save = true)
n = length(mech.storage.x)
for i=1:n
    @test !any(isnan.(mech.storage.x[i]))
    @test !any(isnan.(mech.storage.q[i]))
end

include("examples/slider_crank3d.jl")
simulate!(mech, save = true)
n = length(mech.storage.x)
for i=1:n
    @test !any(isnan.(mech.storage.x[i]))
    @test !any(isnan.(mech.storage.q[i]))
end

include("examples/urdf_doublependulum.jl")
simulate!(mech, save = true)
n = length(mech.storage.x)
for i=1:n
    @test !any(isnan.(mech.storage.x[i]))
    @test !any(isnan.(mech.storage.q[i]))
end

include("examples/wheel.jl")
simulate!(mech, save = true)
n = length(mech.storage.x)
for i=1:n
    @test !any(isnan.(mech.storage.x[i]))
    @test !any(isnan.(mech.storage.q[i]))
end

