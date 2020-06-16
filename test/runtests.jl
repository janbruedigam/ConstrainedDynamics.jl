using Test
using SafeTestsets


include("diff/diff_test.jl")

@safetestset "Factorization Test" begin
    include("factorization_test.jl")
end

@safetestset "Initialization Tests" begin
    include("initialization/initialization_test.jl")
end

@safetestset "Shape Tests" begin
    include("shape_test.jl")
end

@safetestset "Pendulum Period Tests" begin
    include("pendulum_test.jl")
end

@safetestset "Nutation Behavior Test" begin
    include("nutation_test.jl")
end

@safetestset "Example Tests" begin
    include("example_test.jl")
end