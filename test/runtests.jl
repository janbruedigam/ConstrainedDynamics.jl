using Test
using SafeTestsets


@safetestset "Joint Diff Tests" begin
    include("jointdiff_test.jl")
end

@safetestset "Pendulum Period Test" begin
    include("pendulum_test.jl")
end