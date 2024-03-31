function read_shell_container(jld_file::String)

    jldopen(jld_file, "r") do file
          global  shell_container = read(file,"shell_container")
    end
    
    return shell_container
    
    end