using BenchmarkTools

SUITE = BenchmarkGroup()

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

for i=1:length(files)-1
    include("../test/examples/"*files[i]*".jl")
    steps = Base.OneTo(100)
    storage = Storage{Float64}(steps,length(mech.bodies))
    if files[i]=="joint_force" || files[i]=="pendulum_forced" || files[i]=="football" || files[i]=="nutation"
        SUITE[files[i]] = @benchmarkable simulate!($mech, $steps, $storage, $(eval(Meta.parse(files[i]*"_control!")))) samples=1
    elseif files[i]=="chain_in_chain"
        #
    else
        SUITE[files[i]] = @benchmarkable simulate!($mech, $steps, $storage) samples=1
    end
end
