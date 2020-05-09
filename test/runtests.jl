using ConstrainedDynamics
using StaticArrays
using Test


@testset "ConstrainedDynamics.jl" begin
    @test ConstrainedDynamics.skew([1;2;3]) == SMatrix{3,3,Int64,9}(0, 3, -2, -3, 0, 1, 2, -1, 0)
    @test ConstrainedDynamics.skew([1;2;4]) == SMatrix{3,3,Int64,9}(0, 4, -2, -4, 0, 1, 2, -1, 0)
end
