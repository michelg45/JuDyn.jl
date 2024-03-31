function save_element_containers(label::String)

    mc = Main.model_container
    se = Main.SetElements

    if mc.visco == true
        if mc.Shells > 0
        container_name = label*"_shell_container.jld"
        save_shell_container(se.shell_container,container_name)
        end
        if mc.Beams > 0
        container_name = label*"_beam_container.jld"
        save_beam_container(se.beam_container,container_name)
        end
    end

end