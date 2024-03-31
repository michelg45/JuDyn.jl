function read_beam_container(jld_file::String)

    jldopen(jld_file, "r") do file
          global  beam_container = read(file,"beam_container")
    end
    
    return beam_container
    
    end