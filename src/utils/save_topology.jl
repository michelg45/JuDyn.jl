"""
   save_topology

      function to save on o JLD file the three global containers of a model:
      model_container, element_container and node_container.  
      
      calling sequence: save_topology(file_name)

"""
function save_topology(label::String)
  
  jld_file = label*".jld"

  jldopen(jld_file,"w") do file
      write(file,"model_container", Main.model_container)
      write(file,"node_container", Main.node_container)
      write(file,"element_container", Main.element_container)
  end

  println("topology saved on ", jld_file)

  return jld_file 

end
