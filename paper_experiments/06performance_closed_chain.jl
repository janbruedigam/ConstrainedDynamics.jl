using Rotations
using BenchmarkTools

!(@isdefined FullCordDynamics) && include(joinpath("..", "FullCordDynamics.jl"))
using Main.FullCordDynamics

# Parameters
joint_axis = [1.;0.;0.]
m = 1. # mass
l = 1.0 # length
x,y = .1,.1 # size of link
box = Box(x,y,l,m)

# joint connection points
p1 = [0.;0.;l/2]
p2 = -p1

val = zeros(101)

for N = 3:2:11
    # Initial orientation
    ϕout, ϕin = π/4, π/2
    qout = Quaternion(RotX(ϕout))
    qout2 = Quaternion(RotX(ϕout+π))
    qin = Quaternion(RotX(ϕin))


    # Links
    origin = Origin{Float64}()

    link1 = Link(box)
    setInit!(origin,link1,zeros(3),p1,q=qout)
    links = [link1]

    xvals = [link1.x[1]]
    qvals = [[qout for i=1:floor(Int64,N/2)];[qin];[qout2 for i=1:floor(Int64,N/2)]]

    for i=2:floor(Int64,N/2)
        @eval begin
            $(Symbol("link",i)) = Link(box)
            setInit!($links[$i-1],$(Symbol("link",i)),p2,p1,q=$qout)
            push!($links,$(Symbol("link",i)))
            push!($xvals,$(Symbol("link",i)).x[1])
        end
    end

    @eval begin
        $(Symbol("link",ceil(Int64,N/2))) = Link(box)
        setInit!($links[ceil(Int64,$N/2)-1],$(Symbol("link",ceil(Int64,N/2))),p2,p1,q=$qin)
        push!($links,$(Symbol("link",ceil(Int64,N/2))))
        push!($xvals,$(Symbol("link",ceil(Int64,N/2))).x[1])
    end

    for i=ceil(Int64,N/2)+1:N-1
        @eval begin
            $(Symbol("link",i)) = Link(box)
            setInit!($links[$i-1],$(Symbol("link",i)),p2,p1,q=$qout2)
            push!($links,$(Symbol("link",i)))
            push!($xvals,$(Symbol("link",i)).x[1])
        end
    end

    @eval begin
        $(Symbol("link",N)) = Link(box)
        setInit!($origin,$(Symbol("link",N)),[0.;1.;0.],p2,q=$qout2)
        push!($links,$(Symbol("link",N)))
        push!($xvals,$(Symbol("link",N)).x[1])
    end

    # Constraints
    # Remove Axis constraint for spherical joint
    @eval begin
        $(Symbol("joint0to1",N)) = Constraint(Socket($origin,$links[1],zeros(3),p1),Axis($origin,$links[1],joint_axis),SocketYZ($origin,$links[$N],[0.;1.;0.],p2))
        joint12 = Constraint(Socket($links[1],$links[2],p2,p1),Axis($links[1],$links[2],joint_axis))
        constraints = [$(Symbol("joint0to1",N));joint12]
    end
    for i=3:N
        @eval begin
            # Remove Axis constraint for spherical joint
            $(Symbol("joint",i-1,i)) = Constraint(Socket($links[$i-1],$links[$i],p2,p1),Axis($links[$i-1],$links[$i],joint_axis))
            push!($constraints,$(Symbol("joint",i-1,i)))
        end
    end

    # mechanism
    bot = Robot(origin,links, constraints,dt=0.01)

    t = @benchmarkable simulate_reset!($bot,$xvals,$qvals)
    val[N] = BenchmarkTools.minimum(run(t,samples=100,seconds=100)).time
end
