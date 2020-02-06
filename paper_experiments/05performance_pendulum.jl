using Rotations
using BenchmarkTools

!(@isdefined FullCordDynamics) && include(joinpath("..", "FullCordDynamics.jl"))
using Main.FullCordDynamics

# Parameters
joint_axis = [1.;0.;0.]
m = 1. # mass
h = 1. # height
r = .05 # radius
box = Cylinder(r,h,m)

# joint connection points
p1 = [0.;0.;h/2]
p2 = -p1

val = zeros(100)

for N=1:10
    ϕ1 = π/4
    q1 = Quaternion(RotX(ϕ1))
    for i=2:N
        @eval begin
            $(Symbol("q",i)) = Quaternion(RotX($ϕ1))
        end
    end

    # Links
    origin = Origin{Float64}()
    link1 = Link(box)

    xin =[link1.x[1]]
    qin = [q1 for i=1:N]

    links = [link1]

    setInit!(origin,link1,zeros(3),p1,q=q1)
    for i=2:N
        @eval begin
            $(Symbol("link",i)) = Link(box)
            push!($links,$(Symbol("link",i)))
            setInit!($links[$i-1],$links[$i],p2,p1,q=$(Symbol("q",i)))
            push!($xin,$links[$i].x[1])
        end
    end

    # Constraints
    # Remove Axis constraint for spherical joint
    joint01 = Constraint(Socket(origin,link1,zeros(3),p1),Axis(origin,link1,joint_axis))

    constraints = [joint01]

    for i=2:N
        @eval begin
            # Remove Axis constraint for spherical joint
            $(Symbol("joint",i-1,i)) = Constraint(Socket($links[$i-1],$links[$i],p2,p1),Axis($links[$i-1],$links[$i],joint_axis))
            push!($constraints,$(Symbol("joint",i-1,i)))
        end
    end

    # Mechanism
    bot = Robot(origin, links, constraints; tend=10., dt=0.01)

    t = @benchmarkable simulate_reset!($bot,$xin,$qin,tol=1e-10)
    val[N] = BenchmarkTools.minimum(run(t,samples=100,seconds=100)).time
end
