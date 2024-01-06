function read_shapes(jld_file::String)

  
    jldopen(jld_file, "r") do file
        global node_positions = read(file,"node_positions")
        global connected_nodes = read(file,"connected_nodes")
        global triangles= read(file,"triangles")
    end
  
    println("model shapes read from ", jld_file)
  
    return node_positions, connected_nodes, triangles
  
  end