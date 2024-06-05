function read_shape(jld_file::String)

    jldopen(jld_file, "r") do file 

        global node_positions = read(file, "node_positions")
        global triangles = read(file, "triangles")
        global nc = haskey(file,"node_coordinates") 
        nc == true && ( global node_coordinates = read(file, "node_positions"))  

    end
    
    nc == false ? (return node_positions,triangles) : (return node_positions,triangles, node_coordinates) 
    
end