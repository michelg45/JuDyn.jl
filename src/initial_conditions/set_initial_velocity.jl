    """
        set_initial_velocity
    
    
        input:
            node:   node number at which initial velocity is defined
            v:      node asolute translation velocity
            Omega:  node material angular velocity
    
    """
    function set_initial_velocity(node,vt::Vec3,Omega::Vec3)


    
    
    
        nc=Main.node_container
        mc=Main.model_container
        ec=Main.element_container
        Ndofs=copy(mc.Ndofs)
    
        N_nodes=mc.Nodes
        Nel=mc.Elements
    
    
        vit=[vt.v;Omega.v]
    
        inode=findfirst(x -> x == node, nc.node_numbers)
        rv=nc.init_orientations[inode]
        inv_loc=copy(nc.inv_locs[inode])
     #    xdot_0[inv_loc]=vit
    
        for iel=1:Nel
            if ec.element_types[iel] == "rigid_body"
                if ec.element_nodes[iel][1] == inode
                    inv_loc_v= ec.inv_loc_v[iel]
                    x_0[inv_loc_v .+ mc.Ndofs_q]=vit
                end
            end
            if ec.element_types[iel] == "rigid_mass"
                if ec.element_nodes[iel][1] == inode
                    inv_loc_v= ec.inv_loc_v[iel]
                    x_0[inv_loc_v .+ mc.Ndofs_q]=vit[1:3]
                end
            end
        end
    
        n_slaves=size(nc.slave_nodes[inode],1)
    
        for isl = 1:n_slaves
            slave_node = nc.slave_nodes[inode][isl]
            nsl=findfirst(x -> x == nc.slave_nodes[inode][isl], nc.node_numbers)
            inv_loc=copy(nc.inv_locs[nsl][1:3])
            X=nc.init_positions[nsl]
            v_slave=vt+rot(rv,crossp(Omega,X))
            min=findmin(inv_loc)
            if min[1] > 0
     #        xdot_0[inv_loc]=v_slave.v
            nel=findfirst(x -> x == nc.slave_nodes[inode][isl],ec.element_nodes)
            end
            for iel=1:Nel
                if ec.element_types[iel] == "rigid_body"
                    if ec.element_nodes[iel][1] == slave_node
                    inv_loc_vt= copy(ec.inv_loc_v[iel][1:3])
                    inv_loc_vr= copy(ec.inv_loc_v[iel][4:6])
                        x_0[inv_loc_vt]=v_slave.v
                        x_0[inv_loc_vr]=Omega.v
                    end
                end
                        if ec.element_types[iel] == "rigid_mass"
                    if ec.element_nodes[iel][1] == slave_node
                        inv_loc_v= copy(ec.inv_loc_v[iel])
                        x_0[inv_loc_v]=v_slave.v
                    end
                end
            end
    
        end
    
    end