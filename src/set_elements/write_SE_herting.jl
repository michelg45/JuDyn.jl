"""
        write_SE_herting
"""
function write_SE_herting(SE_file::String,SE_name::String, N_bnodes::Int, Nrig::Int, N_B::Int, N_I::Int, connected_components::Vector{Int64}, X_B::Matrix{Float64}, Mrig::Matrix{Float64}, K::Matrix{Float64}, M::Matrix{Float64}, S_rot::Vector{Matrix{Float64}})

        file = jldopen(SE_file,"r+")
        "Nrig" in keys(file) && (delete!(file,"Nrig"))
        "nbr_bdofs" in keys(file) && (delete!(file,"nbr_bdofs"))
        "nbr_modes" in keys(file) && (delete!(file,"nbr_modes"))
        "nbr_bnodes" in keys(file) && (delete!(file,"nbr_bnodes"))
        "node_coordinates" in keys(file) && (delete!(file,"node_coordinates"))
        "node_components" in keys(file) && (delete!(file,"node_components"))
        "Mrig" in keys(file) && (delete!(file,"Mrig"))
        "K" in keys(file) && (delete!(file,"K"))
        "M" in keys(file) && (delete!(file,"Main"))
        "S_rot" in keys(file) && (delete!(file,"S_rot"))
     
        write(file,"SE_name",SE_name)
        write(file,"Nrig",Nrig)
        write(file,"nbr_bdofs",N_B)
        write(file,"nbr_modes",N_I)
        write(file,"nbr_bnodes",N_bnodes)
        write(file,"node_coordinates",X_B)
        write(file,"node_components",connected_components)
        write(file,"Mrig",Mrig)
         write(file,"K",K)
        write(file,"M",M)
        write(file,"S_rot",S_rot)
        close(file)

    end