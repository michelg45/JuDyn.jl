"""
    save_beam_container

    function saving the "SetElements.beam_container"  on a JLD file named "name_beam.jld".

        calling sequence: save_beam_container(beam_container,"name.jld")

"""
function save_beam_container(beam_container,file_name::String)
    jldopen(file_name, "w") do file
    # addrequire(file,"BeamArray.jl")
    write(file,"beam_container",beam_container)
    end
    return
end