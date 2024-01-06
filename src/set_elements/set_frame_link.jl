"""
    set_frame_link


    function which 

    * defines the topology of a rigid connection between nodes node1 and node 2
    * Takes into account the relative position between the nodes
    * defines a set of 6 Lagrange multipliers
    * constructs the array JuDyn.SetElements.frame_link_container[iel] 

    Calling squence: 

        set_frame_link(nbr,node1,node2)

    Input:
        
        nbr: element number
        node1, node2: nodes connected by the frame link element.

"""
function set_frame_link(nbr::Int,node1::Int,node2::Int,scale_factor::Float64)

  
        mc = Main.model_container
        nc = Main.node_container
    
    
    
    
        if mc.Frame_links == 0
            global frame_link_container=FrameLinkArray()
        else
             frame_link_container = SetElements.frame_link_container
        end
    
        flc = frame_link_container
    
        append!(flc.numbers,nbr)
        append!(flc.scale_factor,scale_factor)
    
        n1=findfirst(x -> x==node1, nc.node_numbers)[1]
    
        n2=findfirst(x -> x==node2, nc.node_numbers)[1]
        push!(flc.node_orders,[n1,n2])
        x_1=copy(nc.init_positions[n1])
        x_2=copy(nc.init_positions[n2])
        RV_1 = nc.init_orientations[n1]
        RV_2 = nc.init_orientations[n2]
        RV_0 = RV3(-RV_1,RV_2)
        X_0 =rot(-RV_1,x_2-x_1)
    
    
        push!(flc.relative_position,X_0)
        push!(flc.relative_orientation,RV_0)
    
        mc.Frame_links  += 1
    
        loc_x = [copy(nc.locs[n1]);copy(nc.locs[n2])]
        loc_mult= collect(mc.max_mult +i for i=1:6)
        mc.max_mult +=6
        loc_int = Int[]
        loc_v = Int[]
        append_element(nbr,"frame_link",[n1,n2],loc_x,loc_int,loc_v,loc_mult)
    
    end
    
function set_frame_link(nbr::Int,node1::Int,node2::Int)
        
        scale_factor = 1.0
        set_frame_link(nbr,node1,node2,scale_factor)
    
    end