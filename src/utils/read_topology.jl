function read_topology(jld_file::String)

jldopen(jld_file, "r") do file
      global  node_container = read(file,"node_container")
      global  element_container = read(file,"element_container")
      global  model_container = read(file,"model_container")
end

return node_container, element_container, model_container

end
