"""
        read_SE_herting
"""
function read_SE_herting(SE_file)

        file = jldopen(SE_file)
        SE_name = SE_file[1:5]
        SE_type = "herting"

        Nrig = read(file,"Nrig")
        N_B = read(file,"nbr_bdofs")
        N_I = read(file,"nbr_modes")
        N_bnodes = read(file,"nbr_bnodes")
        X_B = copy(transpose(read(file,"node_coordinates")))
        connected_components = read(file,"node_components")
        Mrig = read(file,"Mrig")
        K =  read(file,"K")
        M = read(file,"M")
        S_rot = read(file,"S_rot")
        
        close(file)

        return SE_name,SE_type, N_bnodes, Nrig, N_B, N_I, connected_components,
         X_B, Mrig, K, M, S_rot

    end
