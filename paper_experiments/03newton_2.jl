using BenchmarkTools

suite = BenchmarkGroup()

suite["matinv"] = BenchmarkGroup(["time", "inv"])

for n=1:10
    # Inverting a dense matrix of the same size as the system in the newton step
    suite["matinv"][string(n)] = @benchmarkable \($(rand(n*9,n*9)),$(rand(n*9)))
end

results = run(suite)
minimum_time = [results.data["matinv"][string(i)].times[1] for i=1:10]
