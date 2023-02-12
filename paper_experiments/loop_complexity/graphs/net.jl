using GraphBasedSystems
using LinearAlgebra
GBS = GraphBasedSystems
using BenchmarkTools

function init1!(system)
    for entry in system.matrix_entries.nzval
        GBS.initialize!(entry)
    end
    for entry in system.vector_entries
        GBS.initialize!(entry)
    end
end

function init2!(A,b)
    A[:,:] = rand(size(A)...)
    b[:] = rand(size(b)[1])
end

timing1 = zeros(18)
timing2 = zeros(18)

for Nlevels = 1:18
    display(Nlevels)

    A = zeros(Int,350,350)
    for i=1:350
        if i+1<=350
            A[i,i+1] = 1
        end
        if i+10<=350
            A[i,i+10] = 1
        end
    end

    nodes = vcat([collect(10*i+1:10*i+Nlevels) for i=0:Nlevels-1]...)
    A = A[nodes,nodes]

    A += A'

    system = System{Float64}(A, ones(Int,Nlevels^2)*6)

    # init1!(system)

    F = full_matrix(system)
    f = full_vector(system)
    f = [f...]
    # ldu_solve!(system)

    # norm(F\f - full_vector(system))

    bm1 = @benchmarkable ldu_solve!($system) setup=(init1!($system))
    timing1[Nlevels] = BenchmarkTools.minimum(run(bm1,samples=100,seconds=100)).time
    bm2 = @benchmarkable ($F\$f) setup=(init2!($F,$f))
    timing2[Nlevels] = BenchmarkTools.minimum(run(bm2,samples=100,seconds=100)).time
end

