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

warningtimes = [
    1000.0 # 562.9592
    50.0 # 24.4849
    500.0 # 275.3259
    4000.0 # 873.5784
    5000.0 # 1249.2964
    4000.0 # 900.5005
    50.0 # 18.3757
    100.0 # 30.044299
    100.0 # 37.420599
    4000.0 # 1688.237801
    500.0 # 227.109
    100.0 # 32.286301
    100.0 # 26.4995
    4000.0 # 2069.365999
    1000.0 # 519.0557
    25.0 # 4.226299
    25.0 # 11.376999
    25.0 # 11.785099
    100.0 # 45.874299
    100.0 # 42.6032
    100.0 # 59.4295
    500.0 # 173.376601
    10.0 # 0.0
]*10

time = zeros(length(files))

c1 = maximum(length.(files))+2
c2 = 9

for i=1:length(files)-1
    include("examples/"*files[i]*".jl")
    steps = Base.OneTo(1000)
    storage = Storage{Float64}(steps,length(mech.bodies))
    if false        
    # elseif files[i]=="joint_force" || files[i]=="pendulum_forced" || files[i]=="football" || files[i]=="nutation"
        # TODO once zero alloc works
        # t = @benchmarkable simulate!($mech, $steps, $storage, $control!)
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
for i=1:length(files)
    l = length(files[i])
    println((" "^(c1-l-1))*files[i]*" | "*string(time[i]))
end

