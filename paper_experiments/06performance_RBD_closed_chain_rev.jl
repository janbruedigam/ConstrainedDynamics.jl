using LinearAlgebra
using RigidBodyDynamics
using StaticArrays
using BenchmarkTools

g = -9.81

l = 1.
m = l
I = 0.0841667

joint_axis = SVector(1., 0., 0.);

function test(state,jointies,N)
    ϕ = π/4
    set_configuration!(state, jointies[1], ϕ)
    for i = 2:floor(Int64,N/2)
        set_configuration!(state, jointies[i], 0.)
    end
    set_configuration!(state, jointies[ceil(Int64,N/2)], π/2-ϕ)
    set_configuration!(state, jointies[ceil(Int64,N/2)+1], π/2+ϕ)
    for i = ceil(Int64,N/2)+2:N-1
        set_configuration!(state, jointies[i], 0.)
    end
    set_configuration!(state, jointies[N], ϕ)

    setdirty!(state)
    simulate(state, 10., Δt = 0.01)
end

val = zeros(101)


for N = 3:2:11
    origin = RigidBody{Float64}("origin")
    chain = Mechanism(origin; gravity = SVector(0., 0., g))

    joint1 = Joint("joint1", Revolute(joint_axis))
    inertia1 = SpatialInertia(frame_after(joint1),com=SVector(0, 0, -l/2),moment_about_com=I*joint_axis*transpose(joint_axis),mass=m)
    link1 = RigidBody(inertia1)
    before_joint1_to_origin = one(Transform3D,frame_before(joint1), default_frame(origin))
    attach!(chain, origin, link1, joint1,joint_pose = before_joint1_to_origin)

    links = [link1]
    jointies = [joint1]

    for i=2:N-1
        @eval begin
            $(Symbol("joint",i)) = Joint("joint"*string($i), Revolute(joint_axis))
            $(Symbol("inertia",i)) = SpatialInertia(frame_after($(Symbol("joint",i))),com=SVector(0, 0, -l/2),moment_about_com=I*joint_axis*transpose(joint_axis),mass=m)
            $(Symbol("link",i)) = RigidBody($(Symbol("inertia",i)))
            $(Symbol("before_joint",i,"_to_after_joint",i-1)) = Transform3D(frame_before($(Symbol("joint",i))), frame_after($jointies[$i-1]), SVector(0, 0., -l))
            attach!($chain, $links[$i-1], $(Symbol("link",i)), $(Symbol("joint",i)),joint_pose = $(Symbol("before_joint",i,"_to_after_joint",i-1)))
            push!($links,$(Symbol("link",i)))
            push!($jointies,$(Symbol("joint",i)))

        end
    end

    @eval begin
        $(Symbol("joint",N)) = Joint("joint"*string($N), Revolute(joint_axis))
        $(Symbol("inertia",N)) = SpatialInertia(frame_after($(Symbol("joint",N))),com=SVector(0, 0, -l/2),moment_about_com=I*joint_axis*transpose(joint_axis),mass=m)
        $(Symbol("link",N)) = RigidBody($(Symbol("inertia",N)))
        $(Symbol("before_joint",N,"_to_origin")) = Transform3D(frame_before($(Symbol("joint",N))), default_frame($origin), SVector(0., 1.,0.))
        attach!($chain, $origin, $(Symbol("link",N)), $(Symbol("joint",N)),joint_pose = $(Symbol("before_joint",N,"_to_origin")))
        push!($links,$(Symbol("link",N)))
        push!($jointies,$(Symbol("joint",N)))

    end

    @eval begin
        $(Symbol("joint",N+1)) = Joint("joint"*string($N+1), Revolute(joint_axis))
        $(Symbol("before_joint",N+1,"_to_joint",N-1)) = Transform3D(frame_before($(Symbol("joint",N+1))), frame_after($jointies[$N-1]), SVector(0, 0., -l))
        $(Symbol("joint",N,"_to_joint",N+1)) = Transform3D(frame_after($jointies[$N]), frame_after($(Symbol("joint",N+1))), SVector(0, 0., l))
        attach!($chain, $links[$N-1], $links[$N], $(Symbol("joint",N+1)),joint_pose = $(Symbol("before_joint",N+1,"_to_joint",N-1)), successor_pose = $(Symbol("joint",N,"_to_joint",N+1)))
        push!($jointies,$(Symbol("joint",N+1)))

    end

    state = MechanismState(chain)

    t = @benchmarkable test($state,$jointies,$N)
    val[N] = BenchmarkTools.minimum(run(t,samples=100,seconds=100)).time
end
