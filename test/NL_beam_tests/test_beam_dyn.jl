function test_beam_dyn()

    model_container, node_container, element_container = create_model()
    println(" ")
    println("nonlinear cantilever beam under dynamic load with NL beams")

    nbr_el =  20

    name = "beam_dyn"

    gravity = Vec3()
    update = true

    set_general_data(name,gravity,update)


    L = 5.0
    stiffness_properties = [1.008000e+09, 3.231000e+08, 3.231000e+08, 1.437700e+07, 9.345000e+06, 9.345000e+06]
    mass_properties= [3.298000e+01, 6.116000e-01, 3.058000e-01, 3.058000e-01]

    input_dir = "test/NL_beam_tests/"

    x = Array{Vec3,1}(undef,nbr_el+1)
    x[1]=Vec3()
    set_node(1,x[1],"frame")
    for i=1:nbr_el
        x[i+1]=Vec3((L*i)/nbr_el, 0., 0.)
        set_node(i+1,x[i+1],"frame")
    end

    end_nodes()

    for i = 1:nbr_el
        connected_nodes = [i,i+1]
        set_beam(i, connected_nodes,stiffness_properties, mass_properties)
    end

    t1 = 0.
    amp = 600000.
    ncomp = 3

    set_node_force(nbr_el+1,nbr_el+1,ncomp,"step",[amp, t1])

    end_elements()

    set_node_BC(1,"clamped")
    end_BC()

    sol_type = "dynamic"

    assemble(sol_type)

    print_model(model_container)

    initial_conditions()
    end_initial_conditions()

    JSON_file = input_dir*"test_beam_dyn.json"

    h5_file = solve(JSON_file,sol_type)

    global times, y, ydot, p, str_energy, pot_energy, kin_energy, ext_work, nitmax, node_positions,model_data = read_results_dynamic(h5_file)

    nodes = [floor(Int,nbr_el/2)+1; nbr_el+1]

    global idx, idy, idz  = select_results_XYZ_nodes(nodes)

    solx = y[:,idx]
    L = solx[1,2]
    solx[:,1] = solx[:,1] .- 0.5*L
    solx[:,2] = solx[:,2] .- L
    soly = y[:,idy]
    solz = y[:,idz]
    println("max. x displacement at mid-length",  findmin(solx[:,1]))
    println("max. y displacement at mid-length",  findmax(soly[:,1]))
    println("max. z displacement at mid-length",  findmax(solz[:,1]))
    println("max. x displacement at tip",  findmin(solx[:,2]))
    println("max. y displacement at tip",  findmax(soly[:,2]))
    println("max. z displacement at tip",  findmax(solz[:,2]))

    res_x = [findmin(solx[:,1])[1], findmin(solx[:,2])[1]]
    res_z = [findmin(solz[:,1])[1], findmin(solz[:,2])[1]]
    res_energy = findmax(-kin_energy-str_energy+ext_work)[1]

    rm(h5_file)

    return res_x, res_z, res_energy
end

