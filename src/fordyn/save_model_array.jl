"""
    save_model_array

    function saving the "JuDyn.FORDYN.model.container"  on a JLD file named "name.jld".

        calling sequence: save_model_array(model_container,"name.jld")

"""
function save_model_array(model_container,file_name::String)
    jldopen(file_name, "w") do file
    addrequire(file,FORDYN)
    write(file,"model",model_container)
    end
    return
end