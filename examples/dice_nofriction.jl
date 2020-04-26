using ConstrainedDynamics


# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 0.5
width, depth = 0.5, 0.5
box1 = Box(width, depth, length1, 1., color = RGBA(1., 1., 0.))
b1 = Box(0.1, 0.1, .1, .1, color = RGBA(1., 0., 0.))

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
links = [link1;[Body(b1) for i = 1:8]]

# Constraints
ineqcs = [InequalityConstraint(Impact(links[i + 1], [0;0;1.0])) for i = 1:8]

joint0to1 = EqualityConstraint(OriginConnection(origin, link1))
eqcs = [joint0to1;[EqualityConstraint(Fixed(link1, links[i + 1], corners[i], zeros(3))) for i = 1:8]]

shapes = [box1;b1]


mech = Mechanism(origin, links, eqcs, ineqcs, shapes = shapes)
setPosition!(mech,link1,x = [0.;-2;1.5])
for i = 1:8
    setPosition!(mech, link1, links[i + 1], p1 = corners[i])
end

setVelocity!(mech,link1,v = [0;3;7.],ω = (rand(3) .- 0.5) * 100)
for i = 1:8
    setVelocity!(mech, link1, links[i + 1], p1 = corners[i])
end


simulate!(mech,save = true)
visualize!(mech)
