"""
    NodeArray
        mutable structure of the Main.node_container hosting the initial data
         at the nodes of the model. It gathers the following information: 
            node_numbers::Array{Int,1}
            parent_nodes::Array{Int,1}
            slave_nodes::Array{Vector{Int},1}
            connected_elements::Array{Vector{Int},1}
            locs::Array{Vector{Int},1}
            inv_locs::Array{Vector{Int},1}
            init_positions::Array{Vec3,1}
            init_orientations::Array{RV3,1}
            types::Array{String,1}
"""
mutable struct  NodeArray

    node_numbers::Array{Int,1}
    parent_nodes::Array{Int,1}
    slave_nodes::Array{Vector{Int},1}
    connected_elements::Array{Vector{Int},1}
    locs::Array{Vector{Int},1}
    inv_locs::Array{Vector{Int},1}
    init_positions::Array{Vec3,1}
    init_orientations::Array{RV3,1}
    types::Array{String,1}

end


"""
    NodeArray()

        function to create the Main.node_container::NodeArray (called by the 'create_model' functionà.)

        NodeArray() = NodeArray(node_numbers,parent_nodes,slave_nodes,connected_elements,
                      locs,inv_locs, init_positions,init_orientations,types)    
"""
NodeArray()=NodeArray(Array{Int,1}[],Array{Int,1}[],Array{Vector{Int},1}[],Array{Vector{Int},1}[],Array{Vector{Int},1}[],Array{Vector{Int},1}[],Array{Vec3,1}[],Array{RV3,1}[],Array{String,1}[])
