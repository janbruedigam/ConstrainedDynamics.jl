using BenchmarkTools

SUITE = BenchmarkGroup()
SUITE["sum"] = @benchmarkable sum($(randn(10_000)))
SUITE["sum2"] = @benchmarkable sum($(randn(10_000)))