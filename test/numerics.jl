using ConstrainedDynamics
using LinearAlgebra

# Parameters
ex = [1.;0.;0.] # joint axis for revolute joint

depth = 1.
width = 1.
height = 1.
mass = 1.
b1 = Box(depth,width,height,mass, color = RGBA(1., 0., 0.)) # link type
# b1.J = [1.0 2 3;4 5 6;7 8 9]; b1.m = 10.0   to change inertial properties for testing

# joint connection points
vert11 = [0.;0.;height / 2]
vert12 = -vert11

# Initial orientation
phi = pi / 4
q1 = UnitQuaternion(RotX(phi))

# Links
N = 100 # Number of links (and joints)

origin = Origin{Float64}()
links = [Body(b1) for i = 1:N]

# Constraints
jointb1 = EqualityConstraint(Revolute(origin, links[1], zeros(3), vert11, ex))
# jointb1 = EqualityConstraint(Spherical(origin, links[1], zeros(3), vert11))
if N>1
    constraints = [jointb1;[EqualityConstraint(Revolute(links[i - 1], links[i], vert12, vert11, ex)) for i = 2:N]]
    # constraints = [jointb1;[EqualityConstraint(Spherical(links[i - 1], links[i], vert12, vert11)) for i = 2:N]]
else
    constraints = [jointb1]
end

shapes = [b1]

mech = Mechanism(origin, links, constraints;shapes = shapes)

# set initial configuration
setPosition!(origin,links[1],p2 = vert11,Δq = q1)
previd = links[1].id
for body in Iterators.drop(mech.bodies, 1)
    global previd
    setPosition!(ConstrainedDynamics.getbody(mech, previd), body, p1 = vert12, p2 = vert11)
    previd = body.id
end

# simulate and visualize
# simulate!(mech,save = true)
# visualize(mech)

# numerics evaluation
ConstrainedDynamics.setentries!(mech) # writes the entries into the A matrix and b vector for the current configuration
# ConstrainedDynamics.factor!(mech.graph,mech.ldu) # inplace factorization of the A matrix
# ConstrainedDynamics.solve!(mech) # solve system (stores result in Δs)
A = ConstrainedDynamics.formAMatrix(mech) # returns the A matrix
b = ConstrainedDynamics.formbVector(mech) # returns the b vector
# Δs = ConstrainedDynamics.formΔsVector(mech) # returns the Δs vector
cond(A)
