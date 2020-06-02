using BenchmarkTools

SUITE = BenchmarkGroup()
SUITE["sum"] = @benchmarkable sum($(randn(10_000)))
SUITE["sum2"] = @benchmarkable sum($(randn(10_000)))

include("examples/atlas.jl")
steps = Base.OneTo(1000)
storage = Storage{Float64}(steps,length(mech.bodies))
SUITE["atlastest"] = @benchmarkable simulate!($mech, $steps, $storage) samples=1