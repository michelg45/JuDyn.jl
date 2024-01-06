function read_superelement(SE_file)

        file = matopen(SE_file)

        SE_type = read(file,"SE_type")
        N_bnodes = convert(Int,read(file,"N_boundary_nodes"))
        N_I = convert(Int,read(file,"N_internal_modes"))
        x = read(file,"Connected_components")
        x = convert(Array{Int,2},x)
        connected_components = Int[x[i] for i in 1:N_bnodes]
        X_B = read(file,"X_B")
        Mrig = read(file,"Mrig")
        Nrig=size(Mrig,1)
        K_BB = read(file,"K_BB")
        K_II = read(file,"K_II")
        N_B=size(K_BB,1)
        U = read(file,"U")
        if SE_type != "CB"
                K_BI = read(file,"K_BI")
                M_BI = read(file,"M_BI")
                M_II = read(file,"M_II")
                M_BB = read(file,"M_BB")
                A = zeros(1,1)
                MU = zeros(1,1)
                MV = zeros(1,1)
                M_vv = zeros(1,1)
        else
                A = read(file,"A")
                MV = read(file,"MV")
                MU = read(file,"MU")
                M_vv = read(file,"M_vv")
                K_BI = zeros(1,1)
                M_BB = zeros(1,1)
                M_II = zeros(1,1)
                M_BI = zeros(1,1)

        end




        close(file)

        SE_name=replace(SE_file,".mat" => "")

        return SE_name,SE_type, N_bnodes, Nrig, N_B, N_I, connected_components, X_B, Mrig, K_BB,M_BB,K_BI,M_BI,K_II,M_II,M_vv,MU,MV, A, U

    end
