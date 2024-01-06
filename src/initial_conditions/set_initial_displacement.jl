"""
    set_initial_displacement
"""
function set_initial_displacement(node,x::Vec3,phi::RV3)

    nc=Main.node_container
    mc=Main.model_container
    Ndofs=copy(mc.Ndofs)
    N_nodes=mc.Nodes

    init_disp=[x.v;phi.v]
    inode=findfirst(x -> x == node, nc.node_numbers)
    inv_loc=copy(nc.inv_locs[inode])
    Dx_0[inv_loc]=init_disp

end