using ConstrainedDynamics
using ConstrainedDynamicsVis
using ConstrainedDynamics: my_constraint
using StaticArrays


# Parameters
length1 = 1.0
width, depth = 0.1, 0.1
box = Box(width, depth, length1, length1)

# Links
origin = Origin{Float64}()
link1 = Body(box)

@inline function g(joint::my_constraint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    a=10
    b=15
    G= SA[xb[1]-xa[1]; ((xb[2]-xa[2])^2/a^2)+((xb[3]-xa[3])^2/b^2)-1]
    return G
end

@inline function g(joint::my_constraint, xb::AbstractVector, qb::UnitQuaternion)
    a=10
    b=15
    G= SA[xb[1]; (xb[2]^2/a^2)+(xb[3]^2/b^2)-1]
    return G  
end

# Constraints
joint_between_origin_and_link1 = EqualityConstraint(my_constraint{Float64,2}(origin,link1,g))


links = [link1]
constraints = [joint_between_origin_and_link1]
shapes = [box]


mech = Mechanism(origin, links, constraints, shapes = shapes)
setPosition!(origin,link1,p2=[5;6;7],Δq = UnitQuaternion(RotX(0.1)))

#solve_Eqc(mech,1e-10,30) 
initializeConstraints!(mech)

storage = simulate!(mech, 10., record = true)
visualize(mech, storage, shapes)