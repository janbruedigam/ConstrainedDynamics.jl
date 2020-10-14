using Test
using SafeTestsets


# include("diff/diff_test.jl")

# @safetestset "Factorization Test" begin
#     include("factorization_test.jl")
# end

# include("initialization/initialization_test.jl")

# @safetestset "Shape Tests" begin
#     include("shape_test.jl")
# end

@safetestset "Dynamics Tests" begin
    include("dynamics/dynamics_test.jl")
end

# @safetestset "Example Tests" begin
#     include("example_test.jl")
# end