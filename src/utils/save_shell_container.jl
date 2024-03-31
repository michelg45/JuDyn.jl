"""
    save_shell_container

    function saving the "SetElements.shell_container"  on a JLD file named "name_shell.jld".

        calling sequence: save_shell_container(shell_container,"name.jld")

"""
function save_shell_container(shell_container,file_name::String)
    jldopen(file_name, "w") do file
    # addrequire(file,ShellArray)
    write(file,"shell_container",shell_container)
    end
    return
end