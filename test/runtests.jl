using JuDyn
using Test

@testset "JuDyn.jl" begin
    include("NL_beam_tests/test_beam_static.jl")
    sol = test_beam_static()
    val = 2.158401e+00
    prec = 1.e-5
    @test (sol[3]- val)/val <= prec

end
