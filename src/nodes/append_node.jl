"""
    append_node
        function called by the  different instances of the"set_nodes" function to append 
        the nodes data in the  "Main.node_container" array collecting nodes data.       
"""
function append_node(nbr::Int,npar::Int,loc::Vector{Int},type::String,init_pos::Vec3,init_rot::RV3)
    nc = Main.node_container
    mc = Main.model_container
    push!(nc.node_numbers,nbr)
    push!(nc.parent_nodes,npar)
    push!(nc.slave_nodes,Vector{Int}[])
    push!(nc.locs,loc)
    push!(nc.inv_locs,loc)
    push!(nc.connected_elements,Vector{Int}[])
    push!(nc.types,type)
    push!(nc.init_positions,init_pos)
    push!(nc.init_orientations,init_rot)
    mc.Nodes+=1
    if npar > 0
        ipar = findfirst(x -> x == npar,nc.node_numbers)
        append!(nc.slave_nodes[ipar],nbr)
    end
end