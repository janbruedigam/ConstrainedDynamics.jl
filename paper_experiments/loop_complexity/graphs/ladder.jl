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

timing1 = zeros(151)
timing2 = zeros(151)

for Nlevels = 1:151
    display(Nlevels)

    A0 = [
        0 1 1 0
        0 0 0 1
        0 0 0 1
        0 0 0 0
    ]

    A = zeros(Int,Nlevels*2,Nlevels*2)
    
    if Nlevels == 1
        A[1:2,1:2] = [0 1;0 0]
    else
        A[1:4,1:4]=A0
    end

    for i=3:Nlevels
        A[(i-1)*2-1:(i-1)*2+1,(i-1)*2+1] = [1;0;0]
        A[(i-1)*2-0:(i-1)*2+2,(i-1)*2+2] = [1;1;0]
    end

    A += A'

    system = System{Float64}(A, ones(Int,Nlevels*2)*6)

    # init(system)

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

