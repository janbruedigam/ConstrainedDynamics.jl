using LinearAlgebra
using ForwardDiff

function gc(xa,qa,xb,qb,a::Vector)
    return [xa+qa+xb+qb;2*xa]
end

function g(x)
    return  gc(x[1],x[2],x[3],x[4],joint)
end

function gall(xa,qa,xb,qb,a::Vector)
    x=[xa;qa;xb;qb]
    global joint=a
    function der(x,a::Vector)
        return ForwardDiff.jacobian(g,x)[:,1:2]
    end
    return der(x,a) 
end
#println(gall(1,2,3,4,[1;2])) 
#println(ForwardDiff.jacobian(x->[2*x[1]],[5]))

#=length1 = 1.0
width, depth = 0.1, 0.1
box = Box(width, depth, length1, length1)
origin = Origin{Float64}()
link1 = Body(box)
constraint1=myconstraint{Float64}(origin,link1)
xa=[0;0;0]
qa=UnitQuaternion([1;0;0;0]...,false)
qb=UnitQuaternion(RotX(0.1))
xb=-vrotate([0;0;0],qb) 

 v=g(constraint1[1], xa, qa, xb, qb)
println(v)

p1=[0.0;0.0;0.0]
    p2=[0.0;0.0;0.5]
    vertices = (p1, p2)
    return vrotate(xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)), inv(qa))=#
