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
    mech.g = 0.0

    steps = Base.OneTo(100)
    storage = Storage{Float64}(steps,length(mech.bodies))

    if files[i]=="joint_force" 
        SUITE[files[i]] = @benchmarkable simulate!($mech, $steps, $storage, $joint_force_control!) samples=1
    elseif files[i]=="pendulum_forced"
        SUITE[files[i]] = @benchmarkable simulate!($mech, $steps, $storage, $pendulum_forced_control!) samples=1
    elseif files[i]=="football"
        SUITE[files[i]] = @benchmarkable simulate!($mech, $steps, $storage, $football_control!) samples=1
    elseif files[i]=="nutation"
        SUITE[files[i]] = @benchmarkable simulate!($mech, $steps, $storage, $nutation_control!) samples=1
    # elseif files[i]=="chain_in_chain"
    #     #
    else
        SUITE[files[i]] = @benchmarkable simulate!($mech, $steps, $storage) samples=1
    end
end
