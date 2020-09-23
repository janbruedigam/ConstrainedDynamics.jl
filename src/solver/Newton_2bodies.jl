using ConstrainedDynamics
using ConstrainedDynamics: Translational, g, params,Rotational,Joint
using LinearAlgebra


function Liniensuche(xold,qold,sx,sq,j)
    
    sqtemp = sq/(2^(j-1)) 
    sxtemp = sx/(2^(j-1))
    w = sqrt(1-norm(sqtemp)^2) 
    sqexpanded = [w;sqtemp]
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

function gfunc2(va,vb,joint) 
    xa = va[1:3]
    qa = UnitQuaternion(va[4:7]...,false)
    xb = vb[1:3]
    qb = UnitQuaternion(vb[4:7]...,false)
    g(joint,xa,qa,xb,qb)
end

function my_g2(va,vb,joint::Joint{T,1}) where T
    joint.V3 * gfunc2(va,vb,joint)
end

function my_g2(va,vb,joint::Joint{T,2}) where T
    joint.V12 * gfunc2(va,vb,joint)
end

function my_g2(va,vb,joint::Joint{T,3}) where T
    gfunc2(va,vb,joint)
end

function my_g2(va,vb,joint)  
    return [my_g2(va,vb,joint[1]) ; my_g2(va,vb,joint[2])]
end

# Derivatives

function generatederivativeb(va,vb,joint)
    xa = va[1:3]
    qa = UnitQuaternion(va[4:7]...,false)
    xb = vb[1:3]
    qb = UnitQuaternion(vb[4:7]...,false)
    return ConstrainedDynamics.∂g∂posb(joint,xa,qa,xb,qb)
end

function generatederivativea(va,vb,joint)
    xa = va[1:3]
    qa = UnitQuaternion(va[4:7]...,false)
    xb = vb[1:3]
    qb = UnitQuaternion(vb[4:7]...,false)
    return ConstrainedDynamics.∂g∂posa(joint,xa,qa,xb,qb)
end

function myderivativea(va,vb,joint::Joint{T,1}) where T 
    X, Q = generatederivativea(va,vb,joint)
    return joint.V3 * [X Q]
end

function myderivativeb(va,vb,joint::Joint{T,1}) where T 
    X, Q = generatederivativeb(va,vb,joint)
    return joint.V3 * [X Q]
end

function myderivativea(va,vb,joint::Joint{T,2}) where T 
    X, Q = generatederivativea(va,vb,joint)
    return joint.V12 * [X Q] 
end

function myderivativeb(va,vb,joint::Joint{T,2}) where T 
    X, Q = generatederivativeb(va,vb,joint)
    return joint.V12 * [X Q] 
end

function myderivativea(va,vb,joint::Joint{T,3}) where T
    X, Q = generatederivativea(va,vb,joint)
    return [X Q]
end

function myderivativeb(va,vb,joint::Joint{T,3}) where T
    X, Q = generatederivativeb(va,vb,joint)
    return [X Q]
end

function myderivative2(va,vb,joint) 
    A=[myderivativea(va,vb,joint[1]);myderivativea(va,vb,joint[2])]
    B=[myderivativeb(va,vb,joint[1]);myderivativeb(va,vb,joint[2])]
    return [A B]
end

function myderivative2(va,vb,joint::Joint) #to test pure translational or rotational
    A=myderivativea(va,vb,joint)
    B=myderivativeb(va,vb,joint)
    return [A B]
end


#Newton für 2 Kôrper 
function newton(f,link1,link2,joint,epsilon,Nmax)

    x0a = [link1.state.xc;params(link1.state.qc)]
    x0b=  [link2.state.xc;params(link2.state.qc)]
    println("initiale Koordinaten a: ",x0a)
    println("initiale Koordinaten b: ",x0b)
    xolda=x0a[1:3]
    qolda=x0a[4:7]
    xnewa=x0a[1:3]
    qnewa=x0a[4:7]
    xoldb=x0b[1:3]
    qoldb=x0b[4:7]
    xnewb=x0b[1:3]
    qnewb=x0b[4:7]
    conv=false

    for k=1:Nmax 
        if (norm(f([xnewa;qnewa],[xnewb;qnewb],joint)) <= epsilon)
            conv=true;
            break
        end

        #step calculation 
        Ma=[Matrix{Float64}(I, 3, 3) zeros(3,4); zeros(3,3) VLᵀmat(qolda)] #6*7
        Mb=[Matrix{Float64}(I, 3, 3) zeros(3,4); zeros(3,3) VLᵀmat(qoldb)] #6*7
        M2=[Ma zeros(6,7);zeros(6,7) Mb] #12*14
        invderiv=pinv(myderivative2([xolda;qolda],[xnewb;qnewb],joint)) #14*n
        Δs=-M2*invderiv*f([xnewa;qnewa],[xnewb;qnewb],joint) #12*1
        sxa=Δs[1:3]
        sqa=Δs[4:6]
        sxb=Δs[7:9]
        sqb=Δs[10:12]

        for j=1:Nmax
            xnewa,qnewa = Liniensuche(xolda,qolda,sxa,sqa,j)
            xnewb,qnewb = Liniensuche(xoldb,qoldb,sxb,sqb,j)
            if norm(f([xnewa;qnewa],[xnewb;qnewb],joint)) < norm(f([xolda;qolda],[xoldb;qoldb],joint)) 
                break
            end
        end

        println("a  ", [xnewa;qnewa])
        println("b  ", [xnewb;qnewb])
        xolda = xnewa
        qolda = qnewa
        xoldb = xnewb
        qoldb = qnewb

    end
    
    return conv,[xnewa;qnewa],[xnewb;qnewb]
end

length1 = 1.0
width, depth = 0.1, 0.1
box = Box(width, depth, length1, length1)

# Links
origin = Origin{Float64}()
link1 = Body(box)
link2 = Body(box)


t,_,_, = Translational{Float64,3}(link1,link2,axis=[1;0;0],p2=[0;0;0.5])
r,_,_, = Rotational{Float64,2}(link1,link2,axis=[1;0;0])


setPosition!(link1,x=[0;0;0]) 
setPosition!(link2,x=[5;6;7]) 
(conv,va,vb)= newton(my_g2,link1,link2,[t r] ,1e-10,10) 