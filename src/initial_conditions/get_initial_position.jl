"""
get_initial_position

    function get_initial_position()

    Initial node positions are collected from node input and assembled in vector x_0
    Initial velocities can be introduced through the function "set_initial_velocity".

"""
function get_initial_position()

  
        nc=Main.node_container
        mc=Main.model_container
        Ndofs=copy(mc.Ndofs)
        N_nodes=copy(mc.Nodes)
    
    
        for in = 1:N_nodes
            inv_loc=copy(nc.inv_locs[in])
            X=copy(nc.init_positions[in])
    
            if nc.types[in] == "linked"
          
                n2 = findfirst(x -> x == nc.parent_nodes[in], nc.node_numbers)[1]
                RV=copy(nc.init_orientations[n2])
                x=rot(RV,X)+nc.init_positions[n2]
    
            else
    
                x=X
                RV=copy(nc.init_orientations[in])
    
            end
            for ic = 1:3
                if inv_loc[ic] > 0
                    x_0[inv_loc[ic]]=x[ic]
                end
                if inv_loc[ic+3] > 0
                    x_0[inv_loc[ic+3]]=RV[ic]
                end
            end
        end
        println("initial position assembled from node input")
    
        return
    end