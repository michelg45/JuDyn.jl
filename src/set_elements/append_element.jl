"""
    append_element

            function called fy all the set_element functions (e.g 'set_beam') to update the content of the 
            'Main.element_container' with the pre-processed data of the current element.

            calling sequence: append_element(nbr,el_type,el_nodes,loc_x,loc_int,loc_v,loc_mult)

            with
                nbr::Int                number of the element.
                el_type::String         type of the element. 
                el_nodes::Vector{Int}   set of nodes of the element.
                loc_x::Vector{Int}      localization vector of "displacement" dofs.
                loc_int::Vector{Int}      localization vector of internal dofs.
                loc_int::Vector{Int}      localization vector of internal dofs.
                loc_mult::Vector{Int}      localization vector of Lagrange multipliers.

"""
function append_element(el_nbr::Int,el_type::String,el_nodes::Vector{Int},loc_x::Vector{Int},loc_int::Vector{Int},loc_v::Vector{Int},loc_mult::Vector{Int})

    ec=Main.element_container
    mc=Main.model_container
    nc=Main.node_container
    val = findfirst(x -> x== el_nbr, ec.element_numbers)
    val !== nothing  &&  error("element ", el_nbr, " already defined")
    mc.end_of_nodes == false && error("element input not allowed before closing node input")

    append!(ec.element_numbers,el_nbr)
    n_nodes=size(el_nodes,1)
    for i=1:n_nodes
        inode = el_nodes[i]
        push!(nc.connected_elements[inode],el_nbr)
    end
    push!(ec.element_types,el_type)
    push!(ec.element_nodes,el_nodes)
    n_x=size(loc_x,1)
    size(loc_v,1)  == 0 ? n_v = 0 : n_v = size(loc_v,1)
    size(loc_mult,1)  == 0 ? n_mult = 0 : n_mult=size(loc_mult,1)
    size(loc_int,1)   == 0 ? n_int = 0 : n_int=size(loc_int,1)
    append!(ec.n_x,n_x)
    append!(ec.n_int,n_int)
    append!(ec.n_v,n_v)
    append!(ec.n_mult,n_mult)
    push!(ec.loc_x,loc_x)
    push!(ec.loc_int,loc_int)
    push!(ec.loc_v,loc_v)
    push!(ec.loc_mult,loc_mult)
    push!(ec.inv_loc_x,copy(loc_x))
    push!(ec.inv_loc_int,copy(loc_int))
    push!(ec.inv_loc_v,copy(loc_v))
    push!(ec.inv_loc_mult,copy(loc_mult))
    mc.Elements+=1


end
