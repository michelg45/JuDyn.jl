function read_element_initial_conditions(label::String)

    mc = Main.model_container
    

    if mc.Shells > 0
        sc = Main.SetElements.shell_container
        container_name = label*"_shell_container"
        shell_container_old = read_shell_container(container_name*".jld")
        sc.strains_g = shell_container_old.strains_g
        sc.visco_strains_g = shell_container_old.visco_strains_g
    end

    if mc.Beams > 0
        sc = Main.SetElements.beam_container
        container_name = label*"_beam_container"
        beam_container_old = read_beam_container(container_name*".jld")
        sc.strains_g = beam_container_old.strains_g
        sc.strains = beam_container_old.strains
    end

end