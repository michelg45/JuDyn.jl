function test_plate_static()

    model_container, node_container, element_container = create_model()


    println(" ")
    println("flat plate with quadrangular shells")


    name = "thin_plate_alu"

    name_1 = name*"_static_2x12_deg1"

    gravity = Vec3(0.0,0.0,0.0)

    rotation = false

    matrix_update = true

    set_general_data(name,gravity,matrix_update)


    #
    # shell_properties
    #

    thickness =  1.0E-02
    stiffness_properties = [7.3E+10, 0.3]
    mass_properties =  [2.7E03]


    input_dir = "test/shell_tests/"

    fname = input_dir*"plate_2x12.msh"

    global num_nodes, x_nodes = read_mesh_nodes(fname)

    Nnodes = size(x_nodes,2)

    for inode =1:Nnodes
        set_node(num_nodes[inode],Vec3(x_nodes[1:3,inode]), "frame")
    end
    end_nodes()

    NumElements, TypeElements, Connected_nodes = read_mesh_elements(fname)

    el_max = findmax(NumElements)[2]

    iquad =  findall(x -> x == 3, TypeElements)

    numElements = NumElements[iquad]

    connected_nodes = Connected_nodes[1:4,iquad]

    nbr_elements= size(iquad,1)

    or = RV3()

    ngauss_points = 4

    for iel = 1:nbr_elements
        set_shell(numElements[iel], connected_nodes[1:4,iel],thickness, 
        stiffness_properties, mass_properties, ngauss_points)
    end


    amp = 1000.0
    axe = Vec3(0.0,0.0,1.0)

    T = 1.0
    set_node_force(el_max+1,3,axe,"linear",[T, amp])
    set_node_force(el_max+2,4,axe,"linear",[T, amp])
    set_node_force(el_max+3,17,axe,"linear",[T, 2*amp])

    # =====================================================================
    # end of element input
    # =====================================================================
    end_elements()
    # =====================================================================
    # set boundary conditions
    # ======================================================================

    i_clamped = num_nodes[findall(x -> x ==0.0, x_nodes[2,:])]

    n_clamped = size(i_clamped,1)

    for i = 1:n_clamped
        set_node_BC(i_clamped[i],"clamped")
    end

    end_BC()

    # =====================================================================
    # model assembly
    # =====================================================================

    sol_type = "static"

    assemble(sol_type,rotation)

    print_model(model_container)

    initial_conditions()
    end_initial_conditions()


    JSON_file2 = input_dir*"plate_static.json"

    # =======================================================================
    h5_file = solve(JSON_file2,sol_type)
    # =======================================================================

    global times, y, ydot, p, str_energy, pot_energy,  ext_work, nitmax,node_positions = read_results_static(h5_file)

    nodes = [3,4,17]

    global ipos, irot  = select_results_nodes(nodes)

    println(name, " : node positions")

    Nnodes_out = size(nodes,1)

    Npas = size(times,1)
    println("Nbr of increments ", Npas-1)

    for i = 1:Nnodes_out
        println("node ", nodes[i], " component u_Z ", y[Npas, ipos[i][3]])
        println("node ", nodes[i], " component phi_X ", y[Npas, irot[i][1]])
    end

    locs = [ipos[i][3] for i in 1:3]
    sol =  y[Npas, locs]

    rm(h5_file)

    return sol 
end