"""
    NodeArray

Data structure of the Main.node_container hosting the initial data
at the nodes of the model. It gathers the following information: 

|  |  |
|:------------------------------------------|----------------------------------|
| node_numbers::Array{Int,1} | node numbering |
| connected_elements::Array{Vector{Int},1} | elements connected to specified node ]
| locs::Array{Vector{Int},1} | node localization vectors | 
| inv_locs::Array{Vector{Int},1} | inverses of node localization vectors 
| init_positions::Array{Vec3,1} |initial position of nodes | 
| init_orientations::Array{RV3,1} | initial node rotation frames |
| types::Array{String,1} | node types ("frame",  "linked", "point_3D") |

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

Function creating an array of NodeArray type.

Calling sequence (once by the create_model function) : 
  
        global node_container = NodeArray()    
"""
NodeArray()=NodeArray(Array{Int,1}[],Array{Int,1}[],Array{Vector{Int},1}[],Array{Vector{Int},1}[],Array{Vector{Int},1}[],Array{Vector{Int},1}[],Array{Vec3,1}[],Array{RV3,1}[],Array{String,1}[])
