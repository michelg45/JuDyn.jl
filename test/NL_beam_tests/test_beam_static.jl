function test_beam_static()
        model_container, node_container, element_container = create_model()
        println(" ")
        println("nonlinear beam under static load with NL beams")
        nbr_elements =  40
        beam = "Cardona"
        name = beam*"_"*string(nbr_elements)*"_NL_elements"

        gravity = Vec3(0.0,0.0,0.0)
        rotation = true
        rotation_speed = Vec3(0.0, 0.0,0.0)
        matrix_update = true
        set_general_data(name,gravity,rotation,rotation_speed,matrix_update)
        #
        # Beam properties
        #
        L = 5.0
        stiffness_properties = [1.008000e+09, 3.231000e+08, 3.231000e+08, 1.437700e+07, 9.345000e+06, 9.345000e+06]
        mass_properties= [3.298000e+01, 6.116000e-01, 3.058000e-01, 3.058000e-01]
        #
        # nodes
        #
        x = Array{Vec3,1}(undef,nbr_elements+1)
        x[1]=Vec3()
        set_node(1,x[1],"frame")
        for i=1:nbr_elements
            x[i+1]=Vec3((L*i)/nbr_elements, 0., 0.)
            set_node(i+1,x[i+1],"frame")
        end
        end_nodes()
        #
        # elements
        #
        or = RV3()
        for i = 1:nbr_elements
            connected_nodes = [i,i+1]
            set_beam(i, connected_nodes,or, stiffness_properties, mass_properties)
        end
        amp = 600000.
        axe = 3
        T = 1.0
        set_node_force(nbr_elements+1,nbr_elements+1,axe,"linear",[T, amp])
        end_elements()
        # 
        # boundary conditions
        # 
        set_node_BC(1,"clamped")
        end_BC()
        # 
        # model assembly
        # 
        sol_type = "static"
        assemble(sol_type,rotation)
        #
        # initial conditions
        #
        initial_conditions()
        end_initial_conditions()
        #
        # static solution
        #
        JSON_file2 = "test/NL_beam_tests/test_beam_static.json"
        h5_file = solve(JSON_file2,sol_type)
        #
        # results selection
        #
        global times, y, p, str_energy, pot_energy,  ext_work, nitmax, model_data = read_results_static(h5_file)
        nodes = [floor(Int,nbr_elements/2)+1; nbr_elements+1]
        global ipos, irot  = select_results_nodes(nodes)
        println(beam, " : node positions")
        Npas = size(times,1)
        println("Nbr of increments ", Npas-1)
        println("node ", nodes[2], " component u_X ", y[Npas, ipos[2][1]])
        println("node ", nodes[2], " component u_Z ", y[Npas, ipos[2][3]])
        println("node ", nodes[2], " component phi_Y ", y[Npas, irot[2][2]])

        global model_container = nothing
        global node_container = nothing
        global element_container = nothing
        global SetElements.beam_container = nothing
        rm("test_NL_beam.h5")

        return y[Npas, ipos[2][1:3]]
    end
    