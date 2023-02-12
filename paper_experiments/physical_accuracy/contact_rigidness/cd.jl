using ConstrainedDynamics


# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 0.5
width, depth = 0.5, 0.5

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
link1 = Box(width, depth, length1, 1., color=RGBA(0,0.3,0.8))

# Constraints
fricsandineqs = [Friction(link1, [0;0;1.0], 0.3; p = corners[i]) for i=1:8]
frics = getindex.(fricsandineqs,1)
ineqcs = vcat(getindex.(fricsandineqs,2)...)

joint0to1 = EqualityConstraint(Floating(origin, link1))

links = [link1]
eqcs = [joint0to1]


mech = Mechanism(origin, links, eqcs, ineqcs, frics)

depth = zeros(100)

function controller!(mechanism, k)
    depth[k] = link1.state.xc[3]
    # display(depth[k])
end

for d in [0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25]
    setPosition!(link1,x = [0.;0;d])
    setVelocity!(link1)
    
    storage = simulate!(mech, 1, controller!; ε = 1e-3, record=true)
    display(minimum(depth))
    # visualize(mech, storage)
end

