# using BenchmarkTools
#
# N=10
#
# for i=1:N
#     @eval begin
#         $(Symbol("A",i)) = rand($i*9,$i*9)
#         $(Symbol("b",i)) = rand($i*9)
#     end
# end
#
# function test(A,b)
#     A\b
# end
#
# meanvec=zeros(N)
#
# for i = 1:N
#     @eval begin
#         $(Symbol("t",i)) = @benchmarkable test($(Symbol("A",i)),$(Symbol("b",i)))
#         tune!($(Symbol("t",i)))
#         meanvec[$i] = BenchmarkTools.minimum(run($(Symbol("t",i)),samples=100)).time
#     end
# end

using BenchmarkTools

suite = BenchmarkGroup()

suite["matinv"] = BenchmarkGroup(["time", "inv"])

for n=1:10
    suite["matinv"][string(n)] = @benchmarkable \($(rand(n*9,n*9)),$(rand(n*9)))
end

results = run(suite)
minvec = [results.data["matinv"][string(i)].times[1] for i=1:10]
plot(minvec)
