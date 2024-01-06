"""
    print_initial_configuration
"""
function print_initial_configuration()

    nc=Main.node_container
    mc=Main.model_container
    Ndofs=copy(mc.Ndofs)
    N_nodes=copy(mc.Nodes)
    ec=Main.element_container
    Nel=mc.Elements

    println("Initial conditions")
    println("==================")
    println()
    for in = 1:N_nodes
        node=copy(nc.node_numbers[in])
        inv_loc=copy(nc.inv_locs[in])
        ncomp=size(inv_loc)[1]
        pos=zeros(ncomp)
        vit=zeros(ncomp)
        for ic = 1:ncomp
            if inv_loc[ic] > 0
                pos[ic]=x_0[inv_loc[ic]]
                vit[ic]=xdot_0[inv_loc[ic]]
            end
        end

        println("Node ", node)
        println("Position / orientation:    ", pos)
        println("Linear / angular velocity: ", vit)
        println()
    end

    for iel=1:Nel
        inv_loc_v=copy(ec.inv_loc_v[iel])
        element=copy(ec.element_numbers[iel])
        ncomp=size(inv_loc_v)[1]
        if ncomp > 0
                vit=zeros(ncomp)
                for ic = 1:ncomp
                    if inv_loc_v[ic] > 0
                    vit[ic]=x_0[inv_loc_v[ic]]
                    end
                end
                println("Element ", element)
                println("velocities: ", vit)
                println()
        end
    end

end