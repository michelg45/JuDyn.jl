
"""
    get_initial_configuration
"""
function get_initial_configuration()

    nc=Main.node_container
    mc=Main.model_container
    Ndofs=copy(mc.Ndofs)
    N_nodes=copy(mc.Nodes)

    x_0=Vector{Float64}(zeros(Ndofs))
    for in = 1:N_nodes
        inv_loc=copy(nc.inv_locs[in])
        x_init=copy(nc.init_positions[in].v)

        rot_init=copy(nc.init_orientations[in].v)
        for ic = 1:3
            if inv_loc[ic] > 0
                x_0[inv_loc[ic]]=x_init[ic]
            end
            if inv_loc[ic+3] > 0
                x_0[inv_loc[ic+3]]=rot_init[ic]
            end
        end
    end
    println("initial configuration assembed from node input")

    return x_0
end