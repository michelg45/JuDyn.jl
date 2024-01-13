
"""
    set_super_element
    
        Function defining the topology of a superelement.
        The set of nodes consists of a reference node  'ref_node' at the floating center of mass, itnitially located at the CM location of the rigid elements, and of NB boundary nodes of coordinates XB. The set of dofs at the boundary nodes can be 3 (translation only) or 6 (translation and rotation).  NI is the number of internal (elastic) modes complementing the model. The superelement matrices will be read from an array 'SEMatrixSet' type containing at set of superelements used in the model.


"""
function set_super_element(nbr::Int,node_ref::Int,N_B::Int,N_I::Int,boundary_nodes::Vector{Int},connected_components::Vector{Int}, X_B::Array{Float64,2},SE_name::String,matrix_set::Int,SE_type::String,nl_corr::Bool)


    mc = Main.model_container
    nbc= Main.node_container


    if mc.SuperElements == 0
        global superelement_container=SuperElementArray()
        sec = superelement_container
    else
        sec = Main.SuperElements.superelement_container
    end

    append!(sec.numbers,nbr)
    n_ref=findfirst(x -> x == node_ref, nbc.node_numbers)[1]
    loc_x = copy(nbc.locs[n_ref])
    loc_v= collect(mc.max_v +i for i=1:6)
    mc.max_v += 6
    N_b_nodes = size(boundary_nodes,1)
    
    append!(sec.ref_node_orders,n_ref)
    append!(sec.internal_mode_numbers,N_I)
    append!(sec.matrix_sets,matrix_set)
    boundary_node_orders=Vector{Int}(undef,N_b_nodes)
    pos_ref = nbc.init_positions[n_ref]
    psi_ref = nbc.init_orientations[n_ref]

    pos_rel_B = Vec3[]
    psi_rel_B = RV3[]

    
    append!(sec.boundary_node_numbers,N_b_nodes)
    
    


    for ino = 1:N_b_nodes
        bnode = boundary_nodes[ino]
        ib =findfirst(x -> x == bnode, nbc.node_numbers)[1]
        boundary_node_orders[ino] = ib
        nc = connected_components[ino]
        loc_x = [loc_x;copy(nbc.locs[ib][1:nc])]
        loc_v= [loc_v;collect(mc.max_v +i for i=1:nc)]
        mc.max_v += nc

        pos_rel_P = rot(-psi_ref, (nbc.init_positions[ib] - pos_ref))

        PREC = 1.0E-04
        if norm(pos_rel_P.v) > PREC && norm(pos_rel_P.v - X_B[ino,:])/norm(pos_rel_P.v) > PREC
            println("rel. norm ", norm(pos_rel_P.v - X_B[ino,:])/norm(pos_rel_P.v))
            error("SE nbr. ",nbr, " geometric incompatibility ", " X_P ", X_B[ino,:], " pos_rel_P ", pos_rel_P,
            "rel. norm ", norm(pos_rel_P.v - X_B[ino,:])/norm(pos_rel_P.v))
        end
        nc > 3 ? psi_rel_P = RV3(-psi_ref, nbc.init_orientations[ib]) : psi_rel_P = RV3(zeros(3))
        push!(pos_rel_B,pos_rel_P)
        push!(psi_rel_B,psi_rel_P)
    end


    push!(sec.superelement_types,SE_type)
    push!(sec.superelement_names,SE_name)
    push!(sec.boundary_node_orders,boundary_node_orders)
    push!(sec.boundary_node_components,connected_components)
    push!(sec.local_node_coordinates,pos_rel_B)
    push!(sec.local_node_orientations,psi_rel_B)
    push!(sec.nl_correction,nl_corr)

    loc_int = collect(mc.max_int +i for i=1:N_I)
    mc.max_int += N_I
    loc_v = [loc_v; collect(mc.max_v +i for i=1:N_I)]
    mc.max_v += N_I
    if SE_type == "CB"
        N_v = size(loc_v,1) - 6
        loc_v = loc_v[1:N_v]
        println("N_v, loc_v ", N_v, loc_v)
        mc.max_v -= 6
        loc_mult = collect(mc.max_mult +i for i=1:6)
        mc.max_mult += 6
    else
        loc_mult = Int[]
    end

    mc.SuperElements +=1

    el_nodes = [n_ref;boundary_node_orders]

    println(el_nodes,loc_x,loc_int,loc_v,loc_mult)

    append_element(nbr,"superelement",el_nodes,loc_x,loc_int,loc_v,loc_mult)

end

function set_super_element(nbr::Int,node_ref::Int,N_B::Int,N_I::Int,boundary_nodes::Vector{Int},connected_components::Vector{Int}, X_B::Array{Float64,2},SE_name::String,matrix_set::Int,SE_type::String)

    nl_corr = true

    """
    function set_super_element : nonlinear correction not applied il nl_correction is not specified

    """

    set_super_element(nbr,node_ref,N_B,N_I,boundary_nodes,connected_components, X_B,SE_name,matrix_set,SE_type,nl_corr)

end
