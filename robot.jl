using StaticArrays

struct Robot{T,Nl,Nc}
    dt::T
    g::T

    nLl::Int64 # Parameters
    nCl::Int64 # Number of constraint equation lines (One constraint typically has several equations)

    links::SVector{Nl,Link{T}}
    constraints::SVector{Nc,Constraint{T}}
end

function Base.show(io::IO, robot::Robot{T,Nl,Nc}) where {T,Nl,Nc}
    heading = string("Robot{",T,",",Nl,",",Nc,"} with ",Nl," links and ",Nc," constraints:")
    links = string("\n Links (links): ",robot.links)
    constraints = string("\n Constraints (constraints): ",robot.constraints)

    print(io,heading,links,constraints)
end


function Robot(linksIn::AbstractVector{<:Link{T}}, constraintsIn::AbstractVector{<:Constraint}; dt::T=.01, g::T=9.81) where T
    Nl = length(linksIn)
    nLl = 6*Nl
    Nc = 0
    nCl = 0

    for (i,link) in enumerate(linksIn)
        link.id[1] = i
        link.dt[1] = dt
        link.g[1] = g
    end

    for constraint in constraintsIn
        Nc += 1
        nCl += getNc(constraint)
    end

    links = convert(SVector{Nl,Link{T}},linksIn)
    constraints = convert(SVector{Nc,Constraint{T}},constraintsIn)


    Robot{T,Nl,Nc}(dt,g,nLl,nCl,links,constraints)
end

Robot(links::AbstractVector{<:Link{T}}; dt=.01, g=9.81) where T = Robot(links, Array{Constraint{T},1}(undef,0), dt=dt, g=g)


function getState(robot::Robot{T,Nl}) where {T,Nl}
    k = 2
    nP = 3

    x = @MVector zeros(T,Nl*nP)
    φ = @MVector zeros(T,Nl)
    for (i,link) in enumerate(robot.links)
        x[sind3(i)] = link.x[k]
        φ[sind1(i)] = angleAxis(link.q[k])[1]*sign(angleAxis(link.q[k])[2][1])
    end

    return x, φ
end

function getDynamics(robot::Robot{T,Nl}) where {T,Nl}
    nP = 3
    fv = @MVector zeros(T,Nl*nP)
    fω = @MVector zeros(T,Nl*nP)
    for (i,link) in enumerate(robot.links)
        fv[sind3(i)] = link.dynT()
        fω[sind3(i)] = link.dynR()
    end

    return [fv-Gxλ(robot);fω-Gqλ(robot)]
end

function Gxλ(robot::Robot{T,Nl}) where {T,Nl}
    nP = 3
    Gλ = @MVector zeros(T,Nl*nP)

    for (i, constraint) in enumerate(robot.constraints)
            id1, id2 = constraint.linkids()
            Gλ[sind3(id1)] += constraint.∂g∂xa()'*constraint.λ
            Gλ[sind3(id2)] += constraint.∂g∂xb()'*constraint.λ
    end

    return Gλ
end

function Gqλ(robot::Robot{T,Nl}) where {T,Nl}
    nP = 3
    Gλ = @MVector zeros(T,Nl*nP)

    for (i, constraint) in enumerate(robot.constraints)
        id1, id2 = constraint.linkids()
        Gλ[sind3(id1)] += constraint.∂g∂qa()'*constraint.λ
        Gλ[sind3(id2)] += constraint.∂g∂qb()'*constraint.λ
    end

    return Gλ
end

function getConstraints(robot::Robot{T}) where T
    nCl = robot.nCl

    g = @MVector zeros(T,nCl)

    startpos = 1
    endpos = 0
    for (i,constraint) in enumerate(robot.constraints)
        endpos += getNc(constraint)
        g[sind(startpos,endpos)] = constraint.g()
        startpos = endpos+1
    end

    return g
end

function initialGuess(robot::Robot{T,Nl}) where {T,Nl}
    nLl = robot.nLl
    nCl = robot.nCl
    nP = 3

    vvec = @MVector zeros(T,Nl*nP)
    ωvec = @MVector zeros(T,Nl*nP)
    λvec = @MVector zeros(T,nCl)

    for (i,link) in enumerate(robot.links)
        vvec[sind3(i)] = link.vnew
        ωvec[sind3(i)] = link.ωnew
    end
    startpos = 1
    endpos = 0
    for (i,constraint) in enumerate(robot.constraints)
        endpos += getNc(constraint)
        λvec[sind(startpos,endpos)] = constraint.λ
        startpos = endpos+1
    end

    return [vvec;ωvec;λvec]
end

function setvars!(robot::Robot{T,Nl},x0::AbstractArray{T}) where {T,Nl}
    nLl = robot.nLl
    nCl = robot.nCl
    nP = 3

    vvec = x0[sind(1,Nl*nP)]
    ωvec = x0[sind(Nl*nP+1,nLl)]
    λvec = x0[sind(nLl+1,nLl+nCl)]

    for (i,link) in enumerate(robot.links)
        link.vnew[sind3(1)] = vvec[sind3(i)]
        link.ωnew[sind3(1)] = ωvec[sind3(i)]
    end
    startpos = 1
    endpos = 0
    for (i,constraint) in enumerate(robot.constraints)
        Nc = getNc(constraint)
        endpos += Nc
        constraint.λ = λvec[sind(startpos,endpos)]
        startpos = endpos+1
    end
end

function updateRobot!(robot::Robot{T,Nl,Nc}, solution::AbstractVector{T}) where {T,Nl,Nc}
    nLl = robot.nLl
    nP = 3
    nCl = robot.nCl

    vnew = convert(SVector{Nl*nP,T},solution[1:Nl*nP])
    ωnew = convert(SVector{Nl*nP,T},solution[Nl*3+1:nLl])
    λnew = convert(SVector{nCl,T},solution[nLl+1:nLl+nCl])

    updateState!(robot,vnew,ωnew)
    updateλ!(robot,λnew)
end

function updateState!(robot::Robot{T}, v2::AbstractVector{T}, ω2::AbstractVector{T}) where T
    nP = 3

    for (i,link) in enumerate(robot.links)
        link.x[1] = link.x[2]
        link.q[1] = link.q[2]

        v2T = v2[sind3(i)]
        ω2T = ω2[sind3(i)]

        link.vnew[sind(1,nP)] = v2T
        link.ωnew[sind(1,nP)] = ω2T

        link.x[2] += v2T.*link.dt[1]
        link.q[2] = Quaternion(link.dt[1]/2 .*Lmat(link.q[2])*ωbar(ω2T,link.dt[1]))
    end
end

function updateλ!(robot::Robot{T},λ::AbstractVector{T}) where T
    startpos=1
    endpos=0
    for (i,constraint) in enumerate(robot.constraints)
        endpos += getNc(constraint)
        constraint.λ = λ[sind(startpos,endpos)]
        startpos = endpos+1
    end
end

function gradient(robot::Robot{T,Nl}) where {T,Nl}
    nLl = robot.nLl
    nCl = robot.nCl
    nP = 3

    Z1 = @SMatrix zeros(T,Nl*nP,Nl*nP)
    Z2 = @SMatrix zeros(T,nCl,nCl)
    ∂dyn∂v, ∂dyn∂ω = dynDerivative(robot)
    ∂g∂v, ∂g∂ω, Gx, Gq = constrDerivative(robot)

    return [∂dyn∂v Z1 -Gx';Z1 ∂dyn∂ω -Gq';∂g∂v ∂g∂ω Z2]
end

function dynDerivative(robot::Robot{T,Nl}) where {T,Nl}
    nLl = robot.nLl
    nP = 3

    ∂dyn∂v = @MMatrix zeros(T,Nl*nP,Nl*nP)
    ∂dyn∂ω = @MMatrix zeros(T,Nl*nP,Nl*nP)

    for (i,link) in enumerate(robot.links)
        ∂dyn∂v[sind3(i),sind3(i)] = link.∂dynT∂vnew()
        ∂dyn∂ω[sind3(i),sind3(i)] = link.∂dynR∂ωnew()
    end

    return ∂dyn∂v, ∂dyn∂ω
end

function constrDerivative(robot::Robot{T,Nl}) where {T,Nl}
    nLl = robot.nLl
    nCl = robot.nCl
    nP = 3

    ∂g∂v = @MMatrix zeros(T,nCl,Nl*nP)
    ∂g∂ω = @MMatrix zeros(T,nCl,Nl*nP)

    Gx = @MMatrix zeros(T,robot.nCl,Nl*nP)
    Gq = @MMatrix zeros(T,robot.nCl,Nl*nP)

    startpos=1
    endpos=0
    for constraint in robot.constraints
        id1, id2 = constraint.linkids()
        endpos += getNc(constraint)

        ∂g∂v[sind(startpos,endpos),sind3(id1)] += constraint.∂g∂va()
        ∂g∂v[sind(startpos,endpos),sind3(id2)] += constraint.∂g∂vb()
        ∂g∂ω[sind(startpos,endpos),sind3(id1)] += constraint.∂g∂ωa()
        ∂g∂ω[sind(startpos,endpos),sind3(id2)] += constraint.∂g∂ωb()
        Gx[sind(startpos,endpos),sind3(id1)] += constraint.∂g∂xa()
        Gx[sind(startpos,endpos),sind3(id2)] += constraint.∂g∂xb()
        Gq[sind(startpos,endpos),sind3(id1)] += constraint.∂g∂qa()
        Gq[sind(startpos,endpos),sind3(id2)] += constraint.∂g∂qb()

        startpos=endpos+1
    end

    return ∂g∂v, ∂g∂ω, Gx, Gq
end


function invDynGrad(robot::Robot{T,Nl}) where {T,Nl}
    nLl = robot.nLl
    nP = 3

    Ginv = @MMatrix zeros(T,nLl,nLl)

    for (i,link) in enumerate(robot.links)
        Ginv[sind3(i),sind3(i)] = inv(link.∂dynT∂vnew())
        Ginv[sind3(i+Nl),sind3(i+Nl)] = inv(link.∂dynR∂ωnew())
    end

    return Ginv
end

sind1(i::T) where T = SVector{1,T}(i)
sindn(i::T,n::T) where T = SVector{n,T}(collect((i-1)*n+1:i*n)...)
sind3(i::T) where T =  SVector{3,T}((i-1)*3+1,(i-1)*3+2,i*3)
sind(startpos::T,endpos::T) where T = SVector{endpos-startpos+1,T}(collect(startpos:endpos)...)
sind(::Constraint{T,Nc}) where {T,Nc} = SVector{Nc,Int64}(collect(1:Nc)...)

getNc(::Constraint{T,Nc}) where {T,Nc} = Nc
