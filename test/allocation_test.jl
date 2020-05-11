using BenchmarkTools

include("examples/atlas.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t))
# @test r.time < ...
@test r.memory == 0
display(r.time)

include("examples/chain_in_chain.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
# @test r.time < ...
@test r.memory == 0

include("examples/dice_nofriction.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
# @test r.time < ...
@test r.memory == 0

include("examples/dice_tiltedplane.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
# @test r.time < ...
@test r.memory == 0

include("examples/dice.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
# @test r.time < ...
@test r.memory == 0

include("examples/disconnected_bodies.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
# @test r.time < ...
@test r.memory == 0

include("examples/doublependulum_3d.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
# @test r.time < ...
@test r.memory == 0

include("examples/inverted_pyramid_plane.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
# @test r.time < ...
@test r.memory == 0

include("examples/joint_force.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
# @test r.time < ...
@test r.memory == 0

include("examples/joint_torque.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
# @test r.time < ...
@test r.memory == 0

include("examples/n_fourbars.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
# @test r.time < ...
@test r.memory == 0

include("examples/n_pendulum.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
# @test r.time < ...
@test r.memory == 0

include("examples/pendulum_forced.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
# @test r.time < ...
@test r.memory == 0

include("examples/pendulum.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
# @test r.time < ...
@test r.memory == 0

include("examples/planar_example.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
# @test r.time < ...
@test r.memory == 0

include("examples/scissor_lift.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
# @test r.time < ...
@test r.memory == 0

include("examples/slider_crank.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
# @test r.time < ...
@test r.memory == 0

include("examples/slider_crank3d.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
# @test r.time < ...
@test r.memory == 0

include("examples/urdf_doublependulum.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
# @test r.time < ...
@test r.memory == 0

# include("examples/wheel.jl")
# t = @benchmarkable simulate!($mech)
# r = BenchmarkTools.minimum(run(t,samples=1))
# @test r.time < ...
# @test r.memory == 0

