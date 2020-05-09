using Test
using SafeTestsets


@safetestset "Joint Diff Tests" begin
    include("jointdiff_test.jl")
end
