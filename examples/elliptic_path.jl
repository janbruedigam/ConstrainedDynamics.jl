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
    a=5
    b=2

    if xb[3]==xa[3]
        xb[2] > xa[2] ? α = -pi/2 : α = pi/2
    elseif xb[3] > xa[3]
        α = -atan(b^2/a^2*(xb[2]-xa[2])/(xb[3]-xa[3]))
    else
        α = -atan(b^2/a^2*(xb[2]-xa[2])/(xb[3]-xa[3])) - pi
    end

    qα = UnitQuaternion(cos(α/2),sin(α/2),0,0,false)

    eqc1=Vmat(qb\qα)  #Orientation relativ to qa ??                                  
    eqc2=SA[xb[1]-xa[1]; ((xb[2]-xa[2])^2/a^2)+((xb[3]-xa[3])^2/b^2)-1]
    G= [eqc1 ; eqc2]
    return G
end

@inline function g(joint::my_constraint, xb::AbstractVector, qb::UnitQuaternion)
    a=5
    b=2
    # function f()
    #     if (xb[2]==0) 
    #         xb[3]>0 ? θ = pi : θ = 0
    #     else
    #         m1=(xb[3]*a^2)/(xb[2]*b^2)                                 #Slope of the normal
    #         α=atan(m1)                                                 #Angle between normal and the horizontal axis 
    #         xb[2]>0 ? θ=(pi/2)+α : θ=-(pi/2)+α                   #Angle between normal and the vertical axis 
    #     end
    #     return sin(θ/2)
    # end

    if xb[3]==0
        xb[2] > 0 ? α = -pi/2 : α = pi/2
    elseif xb[3] > 0
        α = -atan(b^2/a^2*xb[2]/xb[3])
    else
        α = -atan(b^2/a^2*xb[2]/xb[3]) - pi
    end

    qα = UnitQuaternion(cos(α/2),sin(α/2),0,0,false)

    # vb=params(qb)[2:4]
    eqc1=Vmat(qb\qα)
    eqc2=SA[xb[1]; (xb[2]^2/a^2)+(xb[3]^2/b^2)-1]
    G= [eqc1;eqc2]
    return G
end

# Constraints
joint_between_origin_and_link1 = EqualityConstraint(my_constraint{Float64,5}(origin,link1,g))


links = [link1]
constraints = [joint_between_origin_and_link1]
shapes = [box]


mech = Mechanism(origin, links, constraints, shapes = shapes, g=-10)
setPosition!(link1,x=[0;0.01;2],q = UnitQuaternion(RotX(0)))
setVelocity!(link1,v=[0;2.0;0])

initializeConstraints!(mech,newtonIter = 400)

steps = Base.OneTo(1000)
storage = Storage{Float64}(steps,1)

simulate!(mech, storage, record = true)
visualize(mech, storage, shapes)

