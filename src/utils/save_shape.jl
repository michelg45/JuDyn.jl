function save_shape(label::String,node_positions,triangles)

    jld_file = label*".jld"
  
    jldopen(jld_file,"w")  do file 

    file["node_positions"] = node_positions
    file["triangles"] = triangles
      
    println("shape saved on ", jld_file)

    end
  
    return jld_file 
  
  end

  function save_shape(label::String,node_positions,triangles,node_coordinates)

    jld_file = label*".jld"
  
    jldopen(jld_file,"w")  do file 

    file["node_positions"] = node_positions
    file["triangles"] = triangles
    file["node_coordinates"] = node_coordinates
      
    println("shape saved on ", jld_file)

    end
  
    return jld_file 
  
  end