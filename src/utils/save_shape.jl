function save_shape(label::String,node_positions,triangles)

    jld_file = label*".jld"
  
    jldopen(jld_file,"w")  do file 

    file["node_positions"] = node_positions
    file["triangles"] = triangles
      
    println("shape saved on ", jld_file)

    end
  
    return jld_file 
  
  end