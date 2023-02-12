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

timing1 = zeros(6)
timing2 = zeros(6)

for Nlevels = 1:6
    display(Nlevels)

    # net
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
    A0 = A[nodes,nodes]

    # crystal
    N2 = Nlevels^2
    A = [A0 zeros(N2,N2*(Nlevels-1));zeros(N2*(Nlevels-1),N2) zeros(N2*(Nlevels-1),N2*(Nlevels-1))]
    for i=2:Nlevels
        A[(i-1)*N2+1:i*N2,(i-1)*N2+1:i*N2] = A0
        for j=1:N2
            A[N2*(i-1)+j-N2,N2*(i-1)+j] = 1
        end
    end

    A += A'

    system = System{Float64}(A, ones(Int,Nlevels^3)*6)

    init1!(system)

    F = full_matrix(system)
    f = full_vector(system)
    f = [f...]
    ldu_solve!(system)

    norm(F\f - full_vector(system))

    bm1 = @benchmarkable ldu_solve!($system) setup=(init1!($system))
    timing1[Nlevels] = BenchmarkTools.minimum(run(bm1,samples=100,seconds=100)).time
    bm2 = @benchmarkable ($F\$f) setup=(init2!($F,$f))
    timing2[Nlevels] = BenchmarkTools.minimum(run(bm2,samples=100,seconds=100)).time
end