using Test
using SafeTestsets


@safetestset "Joint Diff Tests" begin
    include("jointdiff_test.jl")
end

@safetestset "Pendulum Period Tests" begin
    include("pendulum_test.jl")
end

@safetestset "Example Tests" begin
    include("example_test.jl")
end

# @safetestset "Allocation Tests" begin
#     include("allocation_test.jl")
# end