using JuDyn
using Test

@testset "JuDyn.jl" begin

    prec = 1.e-5

    global model_container, node_container, element_container = create_model()

    include("NL_beam_tests/test_beam_dyn.jl")
    res_x, res_z,res_energy = test_beam_dyn()
    val_x = [-5.242273e-01, -2.041743e+00]
    val_z = [1.358509e+00, 3.602518e+00]
    val_energy = 2.390450e+03
    @test  (((abs(sum(val_x - res_x))+ abs(sum(val_z - res_z))) <= prec && abs(res_energy-val_energy)/val_energy) <= prec)

    model_container, node_container, element_container = create_model()

    include("shell_tests/test_plate_static.jl")
    sol = test_plate_static()
    val = [2.701850e-01, 2.701850e-01, 2.702028e-01]
    @test  sum(broadcast(abs,(sol - val))./val) <= prec

    model_container, node_container, element_container = create_model()

    include("NL_beam_tests/test_beam_static.jl")
    sol = test_beam_static()
    val = 2.158401e+00
    @test (sol[3]- val)/val <= prec

    model_container, node_container, element_container = create_model()
    include("top_tests/test_top_with_frame_link.jl")
    min_pot_energy, max_pot_energy, min_kin_energy, max_kin_energy = test_top_with_frame_link()
    vmin_pot_energy = 1.377719e+01
    vmax_pot_energy = 5.992618e+01
    vmin_kin_energy = 1.537849e+03
    vmax_kin_energy = 1.584051e+03

    @test (abs((vmin_pot_energy-min_pot_energy)/vmin_pot_energy) <= prec &&
        abs((vmin_kin_energy-min_kin_energy)/vmin_kin_energy) <= prec && 
        abs((vmax_pot_energy-max_pot_energy)/vmax_pot_energy) <= prec && 
        abs((vmax_kin_energy-max_kin_energy)/vmax_kin_energy) <= prec )
end
