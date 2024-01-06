"""
    load_model_array

    function recuperating the "JuDyn.FORDYN.model.container" from a JLD file named "name.jld".

        calling sequence: model_container = load_model_array("name.jld")
        
"""
function load_model_array(file_name::String)
    jldopen(file_name, "r") do file
    addrequire(file,FORDYN)
    model_container=read(file,"model")
    end
    return model_container
end