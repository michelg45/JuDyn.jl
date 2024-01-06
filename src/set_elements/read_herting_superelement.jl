function read_herting_superelement(SE_file)

        file = jldopen(SE_file)

        herting_model = read(file,"herting_model")

        SE_name = SE_file[1:5]
        SE_type = "herting"

        Gyro = herting_model.Gyro
        Nrig = herting_model.Nrig
        N_B = herting_model.nbr_boundary_dofs
        N_I = herting_model.nbr_vibration_modes
        N_bnodes = herting_model.nbr_boundary_nodes

        X_B = copy(transpose(herting_model.boundary_node_coordinates))
        
        connected_components = herting_model.boundary_node_components
        Mrig = herting_model.Mrig
        K = herting_model.K
        M = herting_model.M
        S_rot = herting_model.S_rot
        A_rot = herting_model.A_rot
        
        close(file)

        return Gyro,SE_name,SE_type, N_bnodes, Nrig, N_B, N_I, connected_components,
         X_B, Mrig, K, M, S_rot, A_rot

    end
