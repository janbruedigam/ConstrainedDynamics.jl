using ConstrainedDynamics
using ConstrainedDynamicsVis
using ConstrainedDynamics: vrotate,orthogonalrows,Vmat
using StaticArrays

# Parameters
joint_axis = [1.0;0.0;0.0]
length1 = 1.0
width, depth = 0.1, 0.1
box = Box(width, depth, length1, length1)
p2 = [0.0;0.0;length1 / 2] # joint connection point

# Initial orientation
ϕ1 = π / 2
q1 = UnitQuaternion(RotX(ϕ1))

# Links
origin = Origin{Float64}()
link1 = Body(box)

@inline function g(joint::my_constraint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    G= [vrotate(xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)), inv(qa));
        V12*Vmat(qa \ qb / _qoffset)]

    return G
    
end
@inline function g(joint::my_constraint, xb::AbstractVector, qb::UnitQuaternion)
    G= [ xb + vrotate(vertices[2], qb) - vertices[1];
        V12*Vmat(qb / _qoffset)]
    
    return G
    
end
@inline function my_Joint(body1::AbstractBody{T}, body2, axis; p1 = szeros(T, 3), p2 = szeros(T, 3), qoffset = one(UnitQuaternion{T})) where T  
    global vertices = (p1, p2)
    global V1, V2, V3 = orthogonalrows(axis)
    global V12 = [V1;V2]
    global _qoffset=qoffset
    tr=my_constraint{T,5}(body1, body2,g)
    return tr
end

# Constraints
joint_between_origin_and_link1 = EqualityConstraint(my_Joint(origin, link1, joint_axis; p2=p2))


links = [link1]
constraints = [joint_between_origin_and_link1]
shapes = [box]


mech = Mechanism(origin, links, constraints, shapes = shapes)
setPosition!(origin,link1,p2=[5;6;7],Δq = UnitQuaternion(RotX(0.1)))

#solve_Eqc(mech,1e-10,30) 
initializeConstraints!(mech)

storage = simulate!(mech, 10., record = true)
visualize(mech, storage, shapes)
