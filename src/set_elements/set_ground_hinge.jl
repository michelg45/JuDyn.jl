"""
    set_ground_hinge
    
        Defines the topology of a ground hinge connection at  node 'node'
        Takes into account the initial position  and orientation of  the nodes and of the
        relative orientation of axes axis_1 and axis_2
        Defines a set of 3 Lagrange multipliers
        constructs the array hinge_container[iel]

        calling sequence: 

        set_ground_hinge(nbr::Int,node::Int,,axe::Vec3,mode::String)
    
"""
    function set_ground_hinge(nbr::Int,node::Int,axe::Vec3,scale_factor::Float64,mode::String,func::String,params::Vector{Float64})


    
        mc = Main.model_container
        nc= Main.node_container
    
        if mc.Ground_hinges == 0
            global ground_hinge_container=GroundHingeArray()
        else
             ground_hinge_container = SetElements.ground_hinge_container
        end
    
        ghc = ground_hinge_container
    
        append!(ghc.numbers,nbr)
        append!(ghc.scale_factor,scale_factor)
    
        inode =findfirst(x -> x==node, nc.node_numbers)[1]
    
        push!(ghc.node_orders,inode)
        X = copy(nc.init_positions[inode])
        RV = nc.init_orientations[inode]
        axis=RV3(1/norm(axe.v)*axe.v)
    
    
        push!(ghc.position,X)
        push!(ghc.orientation,RV)
        push!(ghc.axis,axis)
        push!(ghc.mode,mode)
        push!(ghc.time_function,func)
        push!(ghc.params,params)
    
    
        mc.Ground_hinges  += 1
    
    
        # nc.locs[inode][1:3] .= 0
        # nc.inv_locs[inode][1:3] .= 0
        loc_x = copy(nc.locs[inode][1:6])
        loc_v = Int[]
    
        if mode == "driven"
            loc_int= Int[]
        else
            loc_int= [(mc.max_int+1)]
            mc.max_int += 1
        end
        loc_mult= collect(mc.max_mult +i for i=1:6)
        mc.max_mult += 6
    
    
        append_element(nbr,"ground_hinge",[inode],loc_x,loc_int,loc_v,loc_mult)
    
    end
    
    function set_ground_hinge(nbr::Int,node::Int,axe::Vec3,mode::String,func::String,params::Vector{Float64})
        scale_factor = 1.0
        set_ground_hinge(nbr,node,axe,scale_factor,mode,func,params)
    end
    
    function set_ground_hinge(nbr::Int,node::Int,axe::Vec3)    
        mode = "free"
        func = " "
        params = [0.0]
        scale_factor = 1.0
        set_ground_hinge(nbr,node,axe,scale_factor,mode,func,params)
    end
    
    function set_ground_hinge(nbr::Int,node::Int,axe::Vec3,mode::String,params::Vector{Float64})
        scale_factor = 1.0
        func = " "
        set_ground_hinge(nbr,node,axe,scale_factor,mode,func,params)
    end