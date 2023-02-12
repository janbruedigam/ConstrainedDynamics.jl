using RigidBodyDynamics
RBD = RigidBodyDynamics
using BenchmarkTools
using LinearAlgebra
using StaticArrays

joint_axis = SVector(1., 0., 0.);
g = -9.81

l = 1.
m = l
I = diagm([0.0845833;0.0845833;0.00125])

# Function to reset joint angles for each run
function reset!(state,jointies,N)
    ϕ = π/4
    set_configuration!(state, jointies[1], ϕ)
    for i = 2:N
        set_configuration!(state, jointies[i], 0.)
    end

    setdirty!(state)
end

function runge_kutta_2(scalar_type::Type{T}) where {T}
    a = zeros(T, 2, 2)
    a[2,1] = 1/2
    b = T[0; 1]
    RBD.ButcherTableau(a, b)
end

function euler_simulate(state0::MechanismState{X}, final_time, control! = RBD.zero_torque!;
        Δt = 1e-4, stabilization_gains=RBD.default_constraint_stabilization_gains(X)) where X
    T = RBD.cache_eltype(state0)
    result = DynamicsResult{T}(state0.mechanism)
    control_torques = similar(velocity(state0))
    closed_loop_dynamics! = let result=result, control_torques=control_torques, stabilization_gains=stabilization_gains # https://github.com/JuliaLang/julia/issues/15276
        function (v̇::AbstractArray, ṡ::AbstractArray, t, state)
            control!(control_torques, t, state)
            dynamics!(result, state, control_torques; stabilization_gains=stabilization_gains)
            copyto!(v̇, result.v̇)
            copyto!(ṡ, result.ṡ)
            nothing
        end
    end
    tableau = runge_kutta_2(T)
    storage = RBD.ExpandingStorage{T}(state0, ceil(Int64, final_time / Δt * 1.001)) # very rough overestimate of number of time steps
    integrator = RBD.MuntheKaasIntegrator(state0, closed_loop_dynamics!, tableau, storage)
    RBD.integrate(integrator, final_time, Δt)
    storage.ts, storage.qs, storage.vs
end

timing1 = zeros(100)
timing2 = zeros(100)

for N = 1:100
    display(N)
    origin = RigidBody{Float64}("origin")
    pendulum = Mechanism(origin; gravity = SVector(0., 0., g))

    joint1 = Joint("joint1", Revolute(joint_axis))
    inertia1 = SpatialInertia(frame_after(joint1),com=SVector(0, 0, -l/2),moment_about_com=I,mass=m)
    link1 = RigidBody(inertia1)
    before_joint1_to_origin = one(Transform3D,frame_before(joint1), default_frame(origin))
    attach!(pendulum, origin, link1, joint1,joint_pose = before_joint1_to_origin)

    jointies = [joint1]
    links = [link1]

    for i = 2:N
        @eval begin
            $(Symbol("joint",i)) = Joint("joint"*string($i), Revolute(joint_axis))
            $(Symbol("inertia",i)) = SpatialInertia(frame_after($(Symbol("joint",i))),com=SVector(0, 0, -l/2),moment_about_com=I,mass=m)
            $(Symbol("link",i)) = RigidBody($(Symbol("inertia",i)))
            $(Symbol("before_joint",i,"_to_after_joint",i-1)) = Transform3D(frame_before($(Symbol("joint",i))), frame_after($jointies[$i-1]), SVector(0, 0., -l))
            attach!($pendulum, $links[$i-1], $(Symbol("link",i)), $(Symbol("joint",i)),joint_pose = $(Symbol("before_joint",i,"_to_after_joint",i-1)))

            push!($links,$(Symbol("link",i)))
            push!($jointies, $(Symbol("joint",i)))
        end
    end

    state = MechanismState(pendulum)

    bn = @benchmarkable euler_simulate($state, 10., Δt = 0.01) setup=(reset!($state,$jointies,$N))
    timing1[N] = BenchmarkTools.minimum(run(bn,samples=100,seconds=100)).time
    bn = @benchmarkable simulate($state, 10., Δt = 0.01) setup=(reset!($state,$jointies,$N))
    timing2[N] = BenchmarkTools.minimum(run(bn,samples=100,seconds=100)).time
end