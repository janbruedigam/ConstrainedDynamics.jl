using Test
using SafeTestsets


@safetestset "Dynamics Diff Tests" begin
    include("dyndiff_test.jl")
end

@safetestset "Joint Diff Tests" begin
    include("jointdiff_test.jl")
end

@safetestset "Pendulum Period Tests" begin
    include("pendulum_test.jl")
end

@safetestset "Initialization Tests" begin
    include("initialize_test.jl")
end

@safetestset "Example Tests" begin
    include("example_test.jl")
end

@safetestset "Allocation Tests" begin
    include("allocation_test.jl")
end