using Rotations
using Plots
using BenchmarkTools

!(@isdefined FullCordDynamics) && include(joinpath("..", "FullCordDynamics.jl"))
using Main.FullCordDynamics

# Parameters
ex = [1.;0.;0.]

l1 = 1.0
x,y = .1,.1
b1 = Box(x,y,l1,l1,color=RGBA(1.,1.,0.))

vert11 = [0.;0.;l1/2]
vert12 = -vert11

# val = zeros(101)

for N = 3:2:11
    # Initial orientation
    phiout, phiin = pi/4, pi/2
    qout = Quaternion(RotX(phiout))
    qout2 = Quaternion(RotX(phiout+pi))
    qin = Quaternion(RotX(phiin))


    # Links
    origin = Origin{Float64}()

    link1 = Link(b1)
    setInit!(origin,link1,zeros(3),vert11,q=qout)
    links = [link1]

    xvals = [link1.x[1]]
    qvals = [[qout for i=1:floor(Int64,N/2)];[qin];[qout2 for i=1:floor(Int64,N/2)]]

    for i=2:floor(Int64,N/2)
        @eval begin
            $(Symbol("link",i)) = Link(b1)
            setInit!($links[$i-1],$(Symbol("link",i)),vert12,vert11,q=$qout)
            push!($links,$(Symbol("link",i)))
            push!($xvals,$(Symbol("link",i)).x[1])
        end
    end

    @eval begin
        $(Symbol("link",ceil(Int64,N/2))) = Link(b1)
        setInit!($links[ceil(Int64,$N/2)-1],$(Symbol("link",ceil(Int64,N/2))),vert12,vert11,q=$qin)
        push!($links,$(Symbol("link",ceil(Int64,N/2))))
        push!($xvals,$(Symbol("link",ceil(Int64,N/2))).x[1])
    end

    for i=ceil(Int64,N/2)+1:N-1
        @eval begin
            $(Symbol("link",i)) = Link(b1)
            setInit!($links[$i-1],$(Symbol("link",i)),vert12,vert11,q=$qout2)
            push!($links,$(Symbol("link",i)))
            push!($xvals,$(Symbol("link",i)).x[1])
        end
    end

    @eval begin
        $(Symbol("link",N)) = Link(b1)
        setInit!($origin,$(Symbol("link",N)),[0.;1.;0.],vert12,q=$qout2)
        push!($links,$(Symbol("link",N)))
        push!($xvals,$(Symbol("link",N)).x[1])
    end

    # Constraints
    @eval begin
        $(Symbol("joint0to1",N)) = Constraint(Socket($origin,$links[1],zeros(3),vert11),SocketYZ($origin,$links[$N],[0.;1.;0.],vert12))
        joint12 = Constraint(Socket($links[1],$links[2],vert12,vert11))
        constraints = [$(Symbol("joint0to1",N));joint12]
    end
    for i=3:N
        @eval begin
            $(Symbol("joint",i-1,i)) = Constraint(Socket($links[$i-1],$links[$i],vert12,vert11))
            push!($constraints,$(Symbol("joint",i-1,i)))
        end
    end

    bot = Robot(origin,links, constraints,dt=0.01)

    t = @benchmarkable simulate!($bot,$xvals,$qvals)
    val[N] = BenchmarkTools.minimum(run(t,samples=100,seconds=100)).time
end
