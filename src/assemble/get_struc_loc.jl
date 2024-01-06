"""
    get_struc_loc
"""
function get_struc_loc(sol_type::String,uniform_rotation::Bool)

    nc=Main.node_container
    mc=Main.model_container
    ec=Main.element_container

    struc_loc_v=Vector{Int}[]
    struc_loc_q=Vector{Int}[]
    idx=0

    if sol_type != "static" &&  sol_type != "dynamic" &&  sol_type != "static_constrained"
        error(" get_struc_loc : incorrect solution type ", sol_type)
    end

    """
    Scanning of node loc vectors
    """

    mc.Ndofs_x=0
    for inode = 1:mc.Nodes
            for ic=1:6
            comp=copy(nc.inv_locs[inode][ic])
            if comp > 0 && mc.Ndofs_x == 0
                struc_loc_q=[struc_loc_q;comp]
                idx+=1
                nc.inv_locs[inode][ic] = idx
                mc.Ndofs_x+=1
            elseif comp > 0
                val = findfirst(x -> x == comp, struc_loc_q)
                if typeof(val) == Nothing
                    struc_loc_q=[struc_loc_q;comp]
                    idx+=1
                    nc.inv_locs[inode][ic] = idx
                    mc.Ndofs_x+=1
                else
                    nc.inv_locs[inode][ic]=val
                end
            end
        end

    end

    """
    Scanning of modified locs and application of loc modification to element loc_x vectors
    """
    nbr_modified = size(mc.modified_locs)[1]

    for i = 1:nbr_modified
        loc_a = mc.modified_locs[i][1]
        loc_b = mc.modified_locs[i][2]
        # println("loc_a, loc_b ", loc_a, loc_b)

        for inel = 1:mc.Elements
            n = findfirst(x -> x==loc_a, ec.loc_x[inel])
            if typeof(n) != Nothing
                ec.loc_x[inel][n] = loc_b
                # println("n, in , loc_b ", n, inel , loc_b)
            end
        end
    end

    """
    Scanning of element loc_x vectors
    """

    for in = 1:mc.Elements
        nc=size(ec.inv_loc_x[in],1)
        for ic=1:nc
            comp=copy(ec.inv_loc_x[in][ic])
            if comp > 0 && mc.Ndofs_x == 0
                struc_loc_q=[struc_loc_q;comp]
                idx+=1
                ec.inv_loc_x[in][ic] = idx
                mc.Ndofs_x+=1
            elseif comp > 0
                val = findfirst(x -> x == comp, struc_loc_q)
                if typeof(val) == Nothing
                    idx+=1
                    struc_loc_q=[struc_loc_q;comp]
                    ec.inv_loc_x[in][ic] = idx
                    mc.Ndofs_x+=1
                else
                    ec.inv_loc_x[in][ic] = val
                end
            end
        end

    end

    mc.Ndofs=copy(mc.Ndofs_x)

    """
    Scanning of element loc_int vectors
    """

        for in = 1:mc.Elements
            if ec.n_int[in] > 0
                nc=size(ec.inv_loc_int[in],1)
                for ic=1:nc
                    comp=copy(ec.inv_loc_int[in][ic])
                    if comp == 0
                        ec.inv_locs_int[in][ic] = 0
                    elseif comp > 0 && mc.Ndofs_int == 0
                        idx+=1
                        ec.inv_loc_int[in][ic]=idx
                        struc_loc_q=[struc_loc_q;comp]
                        mc.Ndofs_int+=1
                    elseif comp > 0
                        # val = findfirst(x -> x == comp, struc_loc)
                        # if typeof(val) == Nothing
                            idx+=1
                            ec.inv_loc_int[in][ic]=idx
                            struc_loc_q=[struc_loc_q;comp]
                            mc.Ndofs_int+=1
                        # end
                    end
                end
            end
        end

        mc.Ndofs+=mc.Ndofs_int

    """
    Scanning of element loc_mult vectors
    """

    for in = 1:mc.Elements
        if ec.n_mult[in] > 0
            nc=size(ec.inv_loc_mult[in],1)
            for ic=1:nc
                comp=copy(ec.inv_loc_mult[in][ic])
                if comp > 0 && mc.Ndofs_mult == 0
                    idx+=1
                    ec.inv_loc_mult[in][ic]=idx
                    struc_loc_q=[struc_loc_q;comp]
                    mc.Ndofs_mult+=1
                elseif comp > 0
                    # val = findfirst(x -> x == comp, struc_loc)
                    # if typeof(val) == Nothing
                        idx += 1
                        ec.inv_loc_mult[in][ic]=idx
                        struc_loc_q=[struc_loc_q;comp]
                        mc.Ndofs_mult += 1
                    # end
                end
            end
        end
    end
    mc.Ndofs += mc.Ndofs_mult

    mc.Ndofs_q = mc.Ndofs




struc_loc_v=Vector{Int}[]
idv =0

println("uniform_rotation ", uniform_rotation)

    for in = 1:mc.Elements
        if ec.n_v[in] > 0
            nc=size(ec.inv_loc_v[in],1)
            maskx  = [ec.inv_loc_x[in]; ec.inv_loc_int[in]]
            for ic=1:nc
                comp=copy(ec.inv_loc_v[in][ic])
                if comp > 0
                    if sol_type == "dynamic" || uniform_rotation == true
                        idv += 1
                        ec.inv_loc_v[in][ic]=idv
                        struc_loc_v=[struc_loc_v;comp]
                        mc.Ndofs_v+=1
                    else
                        ec.inv_loc_v[in][ic] = 0
                    end
                end
"""

    The sequence below would fix the velocity degrees of freedom
    corresponding to a fixed displacement.
    It has beeen temporarily desactivated.


               if comp == 0 || maskx[ic] == 0
                    ec.inv_loc_v[in][ic] = 0
                elseif comp > 0 && mc.Ndofs_v == 0 && maskx[ic] > 0
                    idx+=1
                    ec.inv_loc_v[in][ic]=idx
                    struc_loc=[struc_loc;comp]
                    mc.Ndofs_v+=1
                elseif comp > 0 && maskx[ic] > 0
                    # val = findfirst(x -> x == comp, struc_loc)
                    # if typeof(val) == Nothing
                        idx+=1
                        ec.inv_loc_v[in][ic]=idx
                        struc_loc=[struc_loc;comp]
                        mc.Ndofs_v+=1
                    # end
                end
"""
            end
        end
    end

    # println("element ", in, " inv_loc_v ", ec.inv_loc_v)



    mc.Ndofs += mc.Ndofs_v



    return struc_loc_q, struc_loc_v
end

