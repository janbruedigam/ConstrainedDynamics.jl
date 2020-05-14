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
warningtimes = [
    1000.0 # 562.9592
    500.0 # 275.3259
    2000.0 # 873.5784
    2000.0 # 1249.2964
    2000.0 # 900.5005
    50.0 # 18.3757
    100.0 # 30.044299
    500.0 # 227.109
    100.0 # 32.286301
    100.0 # 26.4995
    4000.0 # 2069.365999
    1000.0 # 519.0557
    25.0 # 11.376999
    25.0 # 11.785099
    100.0 # 45.874299
    100.0 # 42.6032
    100.0 # 59.4295
    500.0 # 173.376601
    100.0 # 37.420599
    10.0 # 0.0
]*10
c1 = maximum(length.(files))+2
c2 = 9

for i=1:19
    include("examples/"*files[i]*".jl")
    steps = Base.OneTo(1000)
    storage = Storage{Float64}(steps,length(mech.bodies))
    if files[i]=="pendulum_forced"
        t = @benchmarkable simulate!($mech, $steps, $storage, $control!)
    # elseif files[i]=="joint_force"
        # TODO once zero alloc works
    else
        t = @benchmarkable simulate!($mech, $steps, $storage)
    end
    r = BenchmarkTools.minimum(run(t, samples=1))
    @test r.memory == 0
    time[i] = r.time/1e6
    @test time[i] < warningtimes[i]
end

println(" "^(c1-7)*" Files | Time")
println("-"^(c1+c2+1))
for i=1:20
    l = length(files[i])
    println((" "^(c1-l-1))*files[i]*" | "*string(time[i]))
end

