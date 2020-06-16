@safetestset "Body Init Tests" begin
    include("body_initialization_test.jl")
end

@safetestset "Fixed Joint Tests" begin
    include("fixed_joint_test.jl")
end

@safetestset "Prismatic Joint Tests" begin
    include("prismatic_joint_test.jl")
end

@safetestset "Revolute Joint Tests" begin
    include("revolute_joint_test.jl")
end
