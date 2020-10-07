using ConstrainedDynamics
using ConstrainedDynamicsVis
using ConstrainedDynamics: my_constraint,Vmat,UnitQuaternion,params
using StaticArrays


# Parameters
length1 = 1.0
width, depth = 0.1, 0.1
box = Box(width, depth, length1, length1)

# Links
origin = Origin{Float64}()
link1 = Body(box)

@inline function g(joint::my_constraint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    a=15
    b=10
    function f()
        if (xb[2]==xa[2]) 
            θ=0
        else
            m1=((xb[3]-xa[3])*a^2)/((xb[2]-xa[2])*b^2)                                 #Slope of the normal
            α=atan(m1)                                                                 #Angle between normal and the horizontal axis 
            (xb[2]-xa[2])*(xb[3]-xa[3])<0 ? θ=(pi/2)+α : θ=-(pi/2)+α                   #Angle between normal and the vertical axis 
        end
        return sin(θ/2)
    end
      
    vb=params(qb)[2:4]                                   
    eqc1=SA[vb[1]-f();vb[2];vb[3]]
    eqc2=SA[xb[1]-xa[1]; ((xb[2]-xa[2])^2/a^2)+((xb[3]-xa[3])^2/b^2)-1]
    G= [eqc1 ; eqc2]
    return G
end

@inline function g(joint::my_constraint, xb::AbstractVector, qb::UnitQuaternion)
    a=15
    b=10
    function f()
        if (xb[2]==0) 
            θ=0
        else
            m1=(xb[3]*a^2)/(xb[2]*b^2)                                 #Slope of the normal
            α=atan(m1)                                                 #Angle between normal and the horizontal axis 
            xb[2]*xb[3]<0 ? θ=(pi/2)+α : θ=-(pi/2)+α                   #Angle between normal and the vertical axis 
        end
        return sin(θ/2)
    end
    vb=params(qb)[2:4]
    eqc1=SA[vb[1]-f();vb[2];vb[3]]
    eqc2=SA[xb[1]; (xb[2]^2/a^2)+(xb[3]^2/b^2)-1]
    G= [eqc1;eqc2]
    return G
end

# Constraints
joint_between_origin_and_link1 = EqualityConstraint(my_constraint{Float64,5}(origin,link1,g))


links = [link1]
constraints = [joint_between_origin_and_link1]
shapes = [box]


mech = Mechanism(origin, links, constraints, shapes = shapes)
setPosition!(origin,link1,p2=[5;6;7],Δq = UnitQuaternion(RotX(0.1)))

initializeConstraints!(mech)

storage = simulate!(mech, 10., record = true)
visualize(mech, storage, shapes)