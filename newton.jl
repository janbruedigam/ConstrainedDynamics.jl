using ForwardDiff
using LinearAlgebra

include("ldu.jl")


function lduinv(A,x)
    F = ldu(A)
    return F\x, F
end

function luinv(A,x)
    F = lu(A)
    return F\x, F
end

function newton(robot::Robot{T,Nl}; ε=1e-10, μ=1e-5, newtonIter=100, refinementIter=20, lineIter=20, warning::Bool=false, refinement::Bool=true) where {T,Nl}
    nLl = robot.nLl # Number of link equation lines (6*Nl = 6*NumberOfLinks)
    nCl = robot.nCl # Number of constraint equation  lines
    x0 = initialGuess(robot) # Current State
    function f(x) # Stores variables in Robot struct and calculates dynamics
        setvars!(robot,x)
        [getDynamics(robot);getConstraints(robot)]
    end

    x1 = x0
    l = length(x0)
    # E = μ*[Matrix(I, nLl, nLl) zeros(nLl,nCl);zeros(nCl,nLl) -Matrix(I, nCl, nCl)]
    E = μ*Matrix(I, l, l)

    n = 1
    for outer n=1:newtonIter
        # Gradient
        f0 = f(x0) # !!! Sets variables in Robot struct. Change calculation order in program with care
        nf0 = norm(f0)

        G = gradient(robot)

        Gμ = G + E

        # Update calculation (x̂ = x0 - x1)
        x̂1, LDUGμ = lduinv(Gμ,f0) # Own LDU seemed to be faster than Julia's lu. Can be replaced with lduinv(Gμ,f0)
        refinement ? x̂1 = iterativeRefinement(x̂1,f0,G,LDUGμ,ε/10;iter=refinementIter, lineIter=lineIter, warning=warning) : nothing

        # Newton update with line search
        x1 = x0 - x̂1
        nf1 = norm(f(x1))
        nf1>nf0 ? x1 = lineSearch(f,x0,nf0,-x̂1;iter=lineIter, warning=warning) : nothing

        if nf1 < ε || norm(x1-x0)<ε
            return x1, n
        else
            x0 = x1
        end
    end

    if warning
        display(string("WARNING:  newton did not converge. n = ",newtonIter,", tol = ",norm(x1-x0),"."))
    end

    return x1, -1
end

function lineSearch(f,x0,nf0,d;iter=20, warning::Bool=false)
    α = 1
    x1 = x0 + 1/(2^α)*d

    for n=1:iter
        α += 1
        norm(f(x1)) >= nf0 ? x1 = x0 + 1/(2^α)*d : (return x1)
    end

    if warning
        display(string("WARNING:  lineSearch did not converge. n = ",iter,"."))
    end

    return x1
end

function iterativeRefinement(x̃0, f0, G, LUGμ, ε; iter=20, lineIter=20, warning::Bool=false)
    f = (x) -> f0 - G*x

    x̃1 = x̃0
    r0 = f(x̃0)
    r1 = r0
    nr0 = norm(r0)

    for n=1:iter
        d = LUGμ\r0
        x̃1 = x̃0 + d
        r1 = f(x̃1)
        nr1 = norm(r1)

        if nr1>nr0
            x̃1 = lineSearch(f,x̃0,nr0,d;iter=lineIter)
            r1 = f(x̃1)
            nr1 = norm(r1)
            if nr1 > nr0
                x̃1 = lineSearch(f,x̃0,nr0,-d;iter=lineIter)
                r1 = f(x̃1)
                nr1 = norm(r1)
            end
        end
        if warning && nr1 > nr0
            display("Warning:  residual increased in iterativeRefinement.")
        end

        if norm(r1 - r0) < ε
            return x̃1
        else
            r0 = r1
            nr0 = nr1
            x̃0 = x̃1
        end
    end

    if warning
        display(string("WARNING:  iterativeRefinement did not converge. n = ",iter,", tol = ",norm(r1 - r0),"."))
    end
    return x̃1
end
