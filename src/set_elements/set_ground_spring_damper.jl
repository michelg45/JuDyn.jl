"""
    set_ground_spring_damper

        function defining a spring damper at a ground node

        calling sequence: 

        set_ground_spring_damper(nbr::Int, node::Int, params::Vector{Float64})
    
        Input:
    
            nbr::Int = element number
            node::Int = node at which force is applied
            params::Vector{Float64} = set of stiffness and damping coefficients 
            [k; c] (isotropic) or [kx; ky; kz; cx; cy; cz] (anisotropic)
            if uniform rotation ==  true, the spring-dammper element must be isotropic and its 
            properties cancellled in the direction of the rotation
    
        Output: entry in SetElements.ground_spring_damper_container::GroundSpringDamperArray
    
"""
function set_ground_spring_damper(nbr::Int, node::Int, params::Vector{Float64})

    
        mc = Main.model_container
        nc = Main.node_container
    
        uniform_rotation = mc.uniform_rotation
      
        if mc.Ground_spring_dampers == 0
            global ground_spring_damper_container = GroundSpringDamperArray()
        else
            ground_spring_damper_container = SetElements.ground_spring_damper_container
        end
    
             gsdc = ground_spring_damper_container
    
    size_par = size(params,1)
    
    size_par != 6 && size_par != 2 &&  (error("el. nbr ", nbr, " spring-damper element - wrong number of parameters"))
    (uniform_rotation == true && size_par != 2)  &&  (error("el. nbr ", nbr, " spring-damper element   in rotation equires only 2 parameters"))
    size_par  == 2 && (params = [params[1]; params[1]; params[1]; params[2]; params[2]; params[2]])
    
    inode = findfirst(x -> x== node, nc.node_numbers)
    typeof(inode) == Nothing && error("element ", nbr, " node ",  "does not exist")
    x_0=copy(nc.init_positions[inode])
    append!(gsdc.number,nbr)
    append!(gsdc.node_order,inode)
    push!(gsdc.position,x_0)
    push!(gsdc.params,params)
    
    # ush!(gsdc.uniform_rotation,uniform_rotation)
    # push!(gsdc.rotation_speed, rotation_speed)
    
    mc.Ground_spring_dampers += 1
    
    loc_x = nc.locs[inode][1:3]
    loc_int = Int[]
    loc_v = Int[]
    loc_mult = Int[]
    append_element(nbr,"ground_spring_damper",[inode],loc_x,loc_int,loc_v,loc_mult)
    
    return
    
    end
    
