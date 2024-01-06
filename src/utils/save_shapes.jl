function save_shapes(jld_file::String,node_positions::Array,connected_nodes::Array,triangles::Array)

  
    jldopen(jld_file,"w") do file
        write(file,"node_positions", node_positions)
        write(file,"connected_nodes", connected_nodes)
        write(file,"triangles", triangles)
    end
  
    println("model shapes saved on ", jld_file)
  
    return 
  
  end