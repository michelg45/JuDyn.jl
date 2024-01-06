"""
    set_node_BC

    allows setting the boundary condtions at a node.

    option 1: application of fixations at node according to mask [x,x,x,x,x,x]
            0 = fixed
            1 = free
        The node boundary conditions are also applied to elements whenever appropriate.
        Input:
            node::Int
            mask::Array{Int,1}
        Calling sequence: 
            set_node_BC(node,mask)

        Option 2: application of BC of type "clamped" or "pinned" at a node.
        Input: 
            node::Int
            BC_type::String = "clamped" or "inned".
        Calling sequence: 
            set_node_BC(node,BC_type) 
                    
"""
function set_node_BC(node::Int,mask::Array{Int,1})

    ec=Main.element_container
    nc=Main.node_container
    mc=Main.model_container

    if size(mask,1) != 6
        error("error of BC mask definition  at node ", node)
    end
    for i =1:6
        if mask[i] != 0 &&  mask[i] != 1
            error("BC error of BC mask definition at node ", node)
        end
    end
    nx=findfirst(x -> x==node, nc.node_numbers)
    nc.inv_locs[nx]=nc.inv_locs[nx].*mask
    unmask = -(mask .- 1)
    fixed_dofs = nc.locs[nx].*unmask
    nfixed, list_fixed_dofs = find_positive_integers(fixed_dofs)
    # println("fixations at node ",nx, " fixed components ", fixed_dofs, "list of fixed dofs ", list_fixed_dofs)
    nel = mc.Elements
    for i=1:nfixed
        dof = list_fixed_dofs[i]
        for iel = 1:nel
            #  inode = findfirst(x -> x == node, ec.element_nodes[iel])
            # if typeof(inode) != Nothing
            #     println(" ec.loc_x[iel]", iel, " ", ec.loc_x[iel])
                ncomp=findfirst(x -> x == dof, ec.loc_x[iel])
                if typeof(ncomp) != Nothing
                    ec.inv_loc_x[iel][ncomp] = 0
                end
            # end
        end
    end

end

function set_node_BC(node::Int,BC_type::String)

"""
function set_node_BC(node::Int,BC_type::String)

    application of BC of type "clamped" or "pinned" at node
        0 = fixed
        1 = free

"""

    if BC_type == "clamped"
        mask = [0,0,0,0,0,0]
        set_node_BC(node,mask)
    elseif BC_type == "pinned"
        mask = [0,0,0,1,1,1]
        set_node_BC(node,mask)
    else
        error("BC type error at node ", node, " ", BC_type)
    end
end
