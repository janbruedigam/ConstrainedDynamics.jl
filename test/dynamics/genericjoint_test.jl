using ConstrainedDynamics
using ConstrainedDynamicsVis
using ConstrainedDynamics: vrotate, Vmat
using Test


joint_axis = [1.0;0.0;0.0]
l = 1.0
w, d = .1, .1
b1 = Box(w, d, l, l)
b2 = Box(w, d, l, l)

vert11 = [0.;0.;l / 2]
vert12 = -vert11
vert21 = [0.;0.;l / 2]

origin = Origin{Float64}()
origin_abstract=Base.deepcopy(origin)
link1 = Body(b1)
link1_abstract=Base.deepcopy(link1)
link2 = Body(b2)
link2_abstract=Base.deepcopy(link2)


@inline function abstract_Fixed(body1::AbstractBody{T}, body2; p1 = szeros(T, 3), p2 = szeros(T, 3), qoffset = one(UnitQuaternion{T})) where T  
     vertices = (p1, p2)
     @inline function g1(joint::my_constraint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)

        G= [vrotate(xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)), inv(qa));
            Vmat(qa \ qb / qoffset)]
        return G
        
    end
    @inline function g1(joint::my_constraint, xb::AbstractVector, qb::UnitQuaternion)
        G= [ xb + vrotate(vertices[2], qb) - vertices[1];
            Vmat(qb / qoffset)]
        return G
        
    end
    return my_constraint{T,6}(body1, body2,g1)
     
end
@inline function abstract_Revolute(body1::AbstractBody{T}, body2, axis; p1 = szeros(T, 3), p2 = szeros(T, 3), qoffset = one(UnitQuaternion{T})) where T  
     vertices = (p1, p2)
     V1, V2, V3 = orthogonalrows(axis)  
     V12 = [V1;V2]
     
     @inline function g2(joint::my_constraint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
       
        G= [vrotate(xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)), inv(qa));
        V12*Vmat(qa \ qb / qoffset)]
        return G
        
    end
    @inline function g2(joint::my_constraint, xb::AbstractVector, qb::UnitQuaternion)
        
        G= [ xb + vrotate(vertices[2], qb) - vertices[1];
        V12*Vmat(qb / qoffset)]
        return G
        
    end
    return my_constraint{T,5}(body1, body2,g2)
    
end
joint0to1 = EqualityConstraint(Fixed(origin, link1; p2=vert11))
joint1to2 = EqualityConstraint(Revolute(link1, link2,joint_axis; p1=vert12, p2=vert21))
joint0to1_abstract = EqualityConstraint(abstract_Fixed(origin_abstract, link1_abstract; p2=vert11))
joint1to2_abstract = EqualityConstraint(abstract_Revolute(link1_abstract, link2_abstract,joint_axis; p1=vert12, p2=vert21))

links = [link1;link2]
links_abstract = [link1_abstract;link2_abstract]
constraints = [joint0to1;joint1to2]
constraints_abstract = [joint0to1_abstract;joint1to2_abstract]
shapes = [b1,b2]

mech = Mechanism(origin, links, constraints, shapes=shapes)
mech_abstract= Mechanism(origin_abstract, links_abstract, constraints_abstract, shapes=shapes)

setPosition!(origin,link1,p2 = vert11,Δq = UnitQuaternion(RotX(pi/4)))
setPosition!(origin_abstract,link1_abstract,p2 = vert11,Δq = UnitQuaternion(RotX(pi/4)))

setPosition!(link1,link2,p1 = vert12,p2 = vert21,Δq = inv(UnitQuaternion(RotX(pi/4)))*UnitQuaternion(RotY(pi/4)))
setPosition!(link1_abstract,link2_abstract,p1 = vert12,p2 = vert21,Δq = inv(UnitQuaternion(RotX(pi/4)))*UnitQuaternion(RotY(pi/4)))


initializeConstraints!(mech) 
initializeConstraints!(mech_abstract)


storage = simulate!(mech, 10., record = true)
visualize(mech, storage, shapes)
storage_abstract = simulate!(mech_abstract, 10., record = true)
visualize(mech_abstract, storage_abstract, shapes)

@test storage==storage_abstract