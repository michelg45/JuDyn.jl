
"""
    init_frames

        Function creating the "current_frames" array and initiating the node reference configuration at first time step.

"""
function init_frames(y::Vector{Float64}, ydot::Vector{Float64},initial_shape::Bool)

    nc = Main.node_container

    Nnodes = Main.model_container.Nodes

    global current_frames=Vector{CurrentFrame}(undef,  Nnodes)
    for i = 1:Nnodes
        current_frames[i] = CurrentFrame()
    end

    cf = current_frames

    for inode=1:Nnodes
        cf[inode].nbr = copy(nc.node_numbers[inode])
        cf[inode].xn = copy(nc.init_positions[inode])
        cf[inode].x = copy(nc.init_positions[inode])
        
        inv_loc = copy(nc.inv_locs[inode])
        cf[inode].inv_loc = copy(inv_loc)
        cf[inode].type = nc.types[inode]
        inv_loc_x = inv_loc[1:3]
        if findmin(inv_loc_x)[1] > 0
            cf[inode].xdot = Vec3(copy(ydot[inv_loc_x]))
            initial_shape == true && (cf[inode].x =  Vec3(copy(y[inv_loc_x])))
            initial_shape == false && (y[inv_loc_x] = (cf[inode].x).v)
        else
            for ic = 1:3
                inv_loc_x[ic] > 0 && (initial_shape == true && (cf[inode].x[ic] = y[inv_loc_x[ic]]); 
                      initial_shape == false && (y[inv_loc_x[ic]] = cf[inode].x[ic]))
            end
        end

        if cf[inode].type == "frame"
            cf[inode].psin = copy(nc.init_orientations[inode])
            cf[inode].psi = copy(nc.init_orientations[inode])
            inv_loc_r = inv_loc[4:6]
                if findmin(inv_loc_r)[1] > 0
                    cf[inode].psidot = Vec3(copy(ydot[inv_loc_r]))
                    initial_shape == true && (cf[inode].psi =  RV3(y[inv_loc_r]))
                    initial_shape == false && (y[inv_loc_r] = (cf[inode].psi).v)
                else
                    for ic=1:3
                        if inv_loc_r[ic] > 0  
                            cf[inode].psidot[ic]=copy(ydot[inv_loc_r[ic]]) 
                            initial_shape == true && (cf[inode].psi[ic] =  y[inv_loc_r[ic]])
                            initial_shape == false  && ( y[inv_loc_r[ic]] = cf[inode].psin[ic] )
                            
                        end 
                    end
                end
                cf[inode].T = tang(cf[inode].psi)
                cf[inode].R = rot(cf[inode].psi)
        end
    end

end  # end initial_shapes
