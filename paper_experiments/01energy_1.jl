using Rotations
using Statistics

!(@isdefined FullCordDynamics) && include(joinpath("..", "FullCordDynamics.jl"))
using Main.FullCordDynamics

# Parameters
joint_axis = [1.;0.;0.]
m = 1. # mass
l = 1. # length
d = .1 # length of one side of a cube
box = Box(d,d,d,m)

p1 = [0.;0.;l] # joint connection point


# Initial orientation
ϕ1, ϕ2 = π/2, 0.
q1, q2 = Quaternion(RotX(ϕ1)), Quaternion(RotX(ϕ2))

# Links
origin = Origin{Float64}()

link1 = Link(box)
setInit!(origin,link1,zeros(3),p1,q=q1)

link2 = Link(box)
setInit!(link1,link2,zeros(3),p1,q=q2)

# Constraints
joint0to1 = Constraint(Socket(origin,link1,zeros(3),p1),Axis(origin,link1,joint_axis))
joint1to2 = Constraint(Socket(link1,link2,zeros(3),p1),Axis(link1,link2,joint_axis))

links = [link1;link2]
constraints = [joint0to1;joint1to2]
shapes = [box]

# Mechanism
bot = Robot(origin, links, constraints, tend=3600., dt=.01)

E = simulate_energy!(bot)
E = E[:,1]+E[:,2] # Kinetic and Potential energies

# Filter noise due to missmatched velcity and position.
# Actual and mean meassured energy are correct.
function filt(ein,w)
    e = copy(ein)
    for i=1:length(e)
        ind1 = maximum([i-floor(Int64,w/2);1])
        ind2 = minimum([i+ceil(Int64,w/2);length(e)])
        e[i] = median(e[ind1:ind2])
    end
    return e
end

E = filt(E,100)
