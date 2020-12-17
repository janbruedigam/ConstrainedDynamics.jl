using ConstrainedDynamics


# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 0.5
width, depth = 0.5, 0.5
box1 = Box(width, depth, length1, 1., color = RGBA(1., 1., 0.))
b1 = Box(0.1, 0.1, .1, .1, color = RGBA(1., 0., 0.))

# Corner vectors
corners = [
    [[length1 / 2;length1 / 2;-length1 / 2]]
    [[length1 / 2;-length1 / 2;-length1 / 2]]
    [[-length1 / 2;length1 / 2;-length1 / 2]]
    [[-length1 / 2;-length1 / 2;-length1 / 2]]
    [[length1 / 2;length1 / 2;length1 / 2]]
    [[length1 / 2;-length1 / 2;length1 / 2]]
    [[-length1 / 2;length1 / 2;length1 / 2]]
    [[-length1 / 2;-length1 / 2;length1 / 2]]
]

# Links
origin = Origin{Float64}()
link1 = Body(box1)

# Constraints
impacts = [InequalityConstraint(Friction(link1,[0;-0.1;1.0], 0.2; p = corners[i])) for i=1:8]

joint0to1 = EqualityConstraint(Floating(origin, link1))

links = [link1]
eqcs = [joint0to1]
ineqcs = impacts
shapes = [box1]


mech = Mechanism(origin, links, eqcs, ineqcs, shapes = shapes)
setPosition!(link1,x = [0.;-2;1.5])

ωtemp = [0.1;0.1;0.1]
setVelocity!(link1,v = [0;3;7.],ω = ωtemp)
