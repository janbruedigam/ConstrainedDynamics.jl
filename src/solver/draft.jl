using ConstrainedDynamics
using ConstrainedDynamics: Translational, g, params,Rotational,Joint,Translational3,Rotational2
using LinearAlgebra


function Liniensuche(xold,qold,sx,sq,j)
    
    sqtemp = sq/(2^(j-1)) 
    sxtemp = sx/(2^(j-1))
    w = sqrt(1-norm(sqtemp)^2) 
    sqexpanded = [w;sq]
    qnew=Lmat(qold)*sqexpanded 
    xnew=xold+sxtemp 
    return xnew,qnew    

end
function VLᵀmat(q)
    [
        -q[2]  q[1]  q[4] -q[3];
        -q[3] -q[4]  q[1]  q[2];
        -q[4]  q[3] -q[2]  q[1];
    ]
end
function Lmat(q)
    [
        q[1]  -q[2]  -q[3] -q[4];
        q[2]   q[1]  -q[4]  q[3];
        q[3]   q[4]   q[1] -q[2];
        q[4]  -q[3]   q[2]  q[1]
    ]
end

# Constraint functions
function gfunc(x,joint) 
    xc = x[1:3]
    qc = UnitQuaternion(x[4:7]...,false)
    g(joint,xc,qc)
end

function my_g(x,joint::Joint{T,1}) where T
    joint.V3 * gfunc(x,joint)
end

function my_g(x,joint::Joint{T,2}) where T
    joint.V12 * gfunc(x,joint)
end

function my_g(x,joint::Joint{T,3}) where T
    gfunc(x,joint)
end

function my_gTR(x,joint)  
    return [my_g(x,joint[1]) ; my_g(x,joint[2])]
end


function generatederivative(x,joint)
    xc = x[1:3]
    qc = UnitQuaternion(x[4:7]...,false)
    return ConstrainedDynamics.∂g∂posb(joint,xc,qc)
end

function myderivative(x,joint::Joint{T,1}) where T 
    X, Q = generatederivative(x,joint)
    return joint.V3 * [X Q]
end

function myderivative(x,joint::Joint{T,2}) where T 
    X, Q = generatederivative(x,joint)
    return joint.V12 * [X Q]
    
end

function myderivative(x,joint::Joint{T,3}) where T
    X, Q = generatederivative(x,joint)
    return [X Q]
end

function myderivativeTR(x,joint) 
    X1 = myderivative(x,joint[1])[:,1:3]
    Q1 = myderivative(x,joint[1])[:,4:7]
    X2 = myderivative(x,joint[2])[:,1:3]
    Q2 = myderivative(x,joint[2])[:,4:7]
    return [X1 Q1;X2 Q2]
    #return [myderivative(x,joint[1]);myderivative(x,joint[2])]
end

function newton(f,link,joint,epsilon,Nmax)

    x0 = [link.state.xc;params(link.state.qc)]
    println("x0",x0)

    xold=x0[1:3]
    qold=x0[4:7]
    xnew=x0[1:3]
    qnew=x0[4:7]
    conv=false

    for k=1:Nmax 
        if (norm(f([xnew;qnew],joint)) <= epsilon)
            conv=true;
            break
        end

        #step calculation
        M=[Matrix{Float64}(I, 3, 3) zeros(3,4); zeros(3,3) VLᵀmat(qold)] #6*7
        #println("M :",M)
        invderiv=pinv(myderivativeTR([xold;qold],joint)) #7*n
        #println("invderiv :",invderiv)
        Δs=-M*invderiv*f([xold;qold],joint) #6*1
        #println("Δs",Δs)
        sx=Δs[1:3]
        sq=Δs[4:6]
        
        for j=1:10
            sqtemp = sq/(2^(j-1)) 
            sxtemp = sx/(2^(j-1))
            w = sqrt(1-norm(sqtemp)^2) 
            sqexpanded = [w;sq]
            qnew=Lmat(qold)*sqexpanded 
            xnew=xold+sxtemp 
            if norm(f([xnew;qnew],joint)) <= norm(f([xold;qold],joint)) 
                break
            end
        end

        println("xnew,qnew",[xnew;qnew])
        xold = xnew
        qold = qnew

    end
    println("xnew,qnew",[xnew;qnew])
    return conv,[xnew;qnew]
end

length1 = 1.0
width, depth = 0.1, 0.1
box = Box(width, depth, length1, length1)

# Links
origin = Origin{Float64}()
link1 = Body(box)
link2 = Body(box)

t,_,_, =Translational3{Float64}(origin, link1, p2 =[0.0;0.0;length1 / 2])
r,_,_, =Rotational2{Float64}(origin, link1, axis = [1.0;0.0;0.0])

setPosition!(origin,link1,p2 = [0;0;length1 / 2],Δq = UnitQuaternion(RotX(π/2)))
(conv,x)= newton(my_gTR,link1,[t r] ,1e-10,10) 



