function read_shape(jld_file::String)

    jldopen(jld_file, "r") do file 

    global node_positions = read(file, "node_positions")
    global triangles = read(file, "triangles")
    end

    
    return node_positions,triangles
    
    end