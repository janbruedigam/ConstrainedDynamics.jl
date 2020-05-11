using BenchmarkTools

time = zeros(20)
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
c1 = maximum(length.(files))+2
c2 = 9

for i=1:19
    include("examples/"*files[i]*".jl")
    t = @benchmarkable simulate!($mech)
    r = BenchmarkTools.minimum(run(t, samples=1))
    @test r.memory == 0
    time[i] = r.time/1e6
end

println(" "^(c1-7)*" Files | Time")
println("-"^(c1+c2+1))
for i=1:20
    l = length(files[i])
    println((" "^(c1-l-1))*files[i]*" | "*string(time[i]))
end

