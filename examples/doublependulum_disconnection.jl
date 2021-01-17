using ConstrainedDynamics
using ConstrainedDynamicsVis


# Parameters
ex = [1.;0.;0.]

h = 1.
r = .05

vert11 = [0.;0.;h / 2]
vert12 = -vert11

# Initial orientation
phi = pi / 4
q1 = UnitQuaternion(RotX(phi))

# Links
origin = Origin{Float64}()
links = [Cylinder(r, h, h, color = RGBA(1., 0., 0.)) for i = 1:2]

# Constraints
jointb1 = EqualityConstraint(Revolute(origin, links[1], ex; p2=vert11))
joint12 = EqualityConstraint(Revolute(links[1], links[2], ex; p1=vert12, p2=vert11))
constraints = [jointb1;joint12]

mech = Mechanism(origin, links, constraints;Δt = 0.01)
setPosition!(origin, links[1], p2 = vert11,Δq = q1)
setPosition!(links[1], links[2], p1 = vert12, p2 = vert11)

function control!(mechanism,k)
    if k==75
        deactivate!(mechanism,constraints[2].id)
    end
end

storage = simulate!(mech, 10., control!,record = true)
visualize(mech, storage)