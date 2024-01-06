"""
    set_node_link

    Function establishing a rigid link between two nodes. Node 2 must be of "linked" type.

    Input:

        * nbr::Int element number
        * node1::Int  "frame" node. 
        * node2::Int  "linked" node.
        * scale_factor::Float64  constraint scaling factor (equal to 1 if omitted).

    Calling sequences:

        set_node_link(nbr, node1,node2,scale_factor)
        set_node_link(nbr, node1,node2)

"""
function set_node_link(nbr::Int,node1::Int,node2::Int,scale_factor::Float64)

    mc = Main.model_container
    nc = Main.node_container


    if mc.Node_links == 0
        global node_link_container=NodeLinkArray()
    else
         node_link_container = SetElements.node_link_container
    end

    nlc = node_link_container

    append!(nlc.numbers,nbr)
    append!(nlc.scale_factor,scale_factor)
    n1=findfirst(x -> x==node1, nc.node_numbers)[1]
    n2=findfirst(x -> x==node2, nc.node_numbers)[1]
    push!(nlc.node_orders,[n1,n2])
    x_1=copy(nc.init_positions[n1])
    x_2=copy(nc.init_positions[n2])
    RV=nc.init_orientations[n2]
    X=rot(-RV,x_1-x_2)
    # nc.init_positions[n1]=copy(x)
    push!(nlc.relative_position,X)
    if nc.types[n1] != "linked"
        error("element node_link ", nbr, " node 1 must be 'linked' type ")
    end

    loc_x = [copy(nc.locs[n1][1:3]);copy(nc.locs[n2])]
    loc_mult= collect(mc.max_mult +i for i=1:3)
    mc.max_mult +=3
    mc.Node_links +=1
    loc_int = Int[]
    loc_v = Int[]
    append_element(nbr,"node_link",[n1,n2],loc_x,loc_int,loc_v,loc_mult)
end

function set_node_link(nbr::Int,node1::Int,node2::Int)
    scale_factor = 1.0
    set_node_link(nbr,node1,node2,scale_factor)
end