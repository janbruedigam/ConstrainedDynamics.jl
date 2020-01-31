using Rotations
using Plots

!(@isdefined FullCordDynamics) && include(joinpath("..", "FullCordDynamics.jl"))
using Main.FullCordDynamics

# Parameters
ex = [1.;0.;0.]

m=1.
l1 = 1.
l2 = l1#sqrt(2)/2
x = .001
b1 = Box(x,x,x,m,color=RGBA(1.,1.,0.))
b2 = Box(x,x,x,m,color=RGBA(1.,1.,0.))

vert11 = [0.;0.;l1]
vert12 = zeros(3)

vert21 = [0.;0.;l2]

# Initial orientation
phi1, phi2 = pi/2, 0.
q1, q2 = Quaternion(RotX(phi1)), Quaternion(RotX(phi2))

# Links
origin = Origin{Float64}()

link1 = Link(b1)
setInit!(origin,link1,zeros(3),vert11,q=q1)

link2 = Link(b2)
setInit!(link1,link2,vert12,vert21,q=q2)

# Constraints
joint0to1 = Constraint(Socket(origin,link1,zeros(3),vert11),Axis(origin,link1,ex))
joint1to2 = Constraint(Socket(link1,link2,vert12,vert21),Axis(link1,link2,ex))

links = [link1;link2]
constraints = [joint0to1;joint1to2]
shapes = [b1,b2]


bot = Robot(origin,links, constraints,tend=3600.,dt=.01)

E1 = simulate!(bot,save=false,debug=false)
E2 = E1[:,1]+E1[:,2]

function filt(ein,w)
    e = copy(ein)
    for i=1:length(e)
        ind1 = maximum([i-floor(Int64,w/2);1])
        ind2 = minimum([i+ceil(Int64,w/2);length(e)])
        e[i] = median(e[ind1:ind2])
    end
    return e
end

E3 = filt(E2,100)
