files = [
    "atlas"
    "axes"
    "chain_in_chain"
    "dice_nofriction"
    "dice_tiltedplane"
    "dice"
    "disconnected_bodies"
    "doublependulum_3d"
    "doublependulum_disconnection"
    "doublependulum_urdf"
    "elliptic_joint"
    "fourbar_disconnection"
    "inverted_pyramid_plane"
    "joint_force"
    "joint_torque"
    "linear_pendulum"
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

controlled = [
    "doublependulum_disconnection"
    "joint_force"
    "joint_torque"
    "pendulum_forced"
    "nutation"
    "fourbar_disconnection"
    "slider_crank"
]

for i=1:length(files)
    include("../test/examples/"*files[i]*".jl")
    mech.g = 0.0

    steps = Base.OneTo(100)
    storage = Storage{Float64}(steps,length(mech.bodies))

    if files[i] ∈ controlled
        control_function = eval(Meta.parse(files[i]*"_control!"))
        SUITE[files[i]] = @benchmarkable simulate!($mech, $steps, $storage, $control_function) samples=2
    else
        SUITE[files[i]] = @benchmarkable simulate!($mech, $steps, $storage) samples=2
    end
end
