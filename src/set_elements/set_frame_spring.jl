"""
    set_frame_spring
 
        Function defining the topology of a spring connection between nodes node1 and node 2
            * Takes into account the relative position between the nodes
            * Constructs the array frame_spring_container[iel]

"""
function set_frame_spring(nbr::Int,node1::Int,node2::Int,k_x::Vec3,k_psi::Vec3)

    Main.model_container=Main.model_container
    node_container=node_container
#     frame_spring_container=frame_spring_container

     if Main.model_container.Frame_springs == 0
         global frame_spring_container=FrameSpringArray()
     else
         frame_spring_container=SetElements.frame_spring_container
    end

    append!(frame_spring_container.numbers,nbr)
    n1=findfirst(x -> x==node1, node_container.node_numbers)[1]
    n2=findfirst(x -> x==node2, node_container.node_numbers)[1]
    push!(frame_spring_container.node_orders,[n1,n2])
    x_1=copy(node_container.init_positions[n1])
    x_2=copy(node_container.init_positions[n2])
    RV_1 = node_container.init_orientations[n1]
    RV_2 = node_container.init_orientations[n2]
    RV_0 = RV3(-RV_1,RV_2)
    X_0 =rot(-RV_1,x_2-x_1)

    push!(frame_spring_container.extension_stiffness,k_x)
    push!(frame_spring_container.rotation_stiffness,k_psi)

    push!(frame_spring_container.initial_displacements,X_0)
    push!(frame_spring_container.initial_rotations,RV_0)

    Main.model_container.Frame_springs +=1

    loc_x = [copy(node_container.locs[n1]);copy(node_container.locs[n2])]
    loc_v = Int[]
    loc_int = Int[]
    loc_mult = Int[]
    append_element(nbr,"frame_spring",[node1,node2],loc_x,loc_int,loc_v,loc_mult)

end