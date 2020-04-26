using BenchmarkTools

time = zeros(21)
memory = zeros(21)

include("atlas.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
time[1] = r.time
memory[1] = r.memory

include("chain_in_chain.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
time[2] = r.time
memory[2] = r.memory

include("dice_nofriction.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
time[3] = r.time
memory[3] = r.memory

include("dice_tiltedplane.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
time[4] = r.time
memory[4] = r.memory

include("dice.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
time[5] = r.time
memory[5] = r.memory

include("disconnected_bodies.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
time[6] = r.time
memory[6] = r.memory

include("doublependulum_3d.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
time[7] = r.time
memory[7] = r.memory

include("inverted_pyramid_plane.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
time[8] = r.time
memory[8] = r.memory

include("joint_force.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
time[9] = r.time
memory[9] = r.memory

include("joint_torque.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
time[10] = r.time
memory[10] = r.memory

include("n_fourbars.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
time[11] = r.time
memory[11] = r.memory

include("n_pendulum.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
time[12] = r.time
memory[12] = r.memory

include("pendulum_controlled.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
time[13] = r.time
memory[13] = r.memory

include("pendulum_forced.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
time[14] = r.time
memory[14] = r.memory

include("pendulum.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
time[15] = r.time
memory[15] = r.memory

include("planar_example.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
time[16] = r.time
memory[16] = r.memory

include("scissor_lift.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
time[17] = r.time
memory[17] = r.memory

include("slider_crank.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
time[18] = r.time
memory[18] = r.memory

include("slider_crank3d.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
time[19] = r.time
memory[19] = r.memory

include("urdf_doublependulum.jl")
t = @benchmarkable simulate!($mech)
r = BenchmarkTools.minimum(run(t, samples = 1))
time[20] = r.time
memory[20] = r.memory

# include("wheel.jl")
# t = @benchmarkable simulate!($mech)
# r = BenchmarkTools.minimum(run(t,samples=1))
# time[21] = r.time
# memory[21] = r.memory

time /= 1e6
