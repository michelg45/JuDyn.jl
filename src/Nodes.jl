__precompile__()

"""
    Nodes

    The module "Nodes"  allows the definition and management of nodal frame data.
    Nodal information is stored in the "Main.node_container" (NodeArray structure).
    Nodes are introduced using the "set_node" function.
"""
module Nodes



using LinearAlgebra
using ..MyAlgebra

dir = "nodes/"


include(dir*"nodeframe_loc.jl")

include(dir*"append_node.jl")

include(dir*"set_node.jl")

include(dir*"set_node_connection.jl")

include(dir*"struct_loc.jl")

include(dir*"print_node.jl")

include(dir*"select_results_XYZ_nodes.jl")

include(dir*"select_results_nodes.jl")

include(dir*"find_node_components.jl")

include(dir*"end_nodes.jl")



# ===========================================================================

# ==============================================

export nodeframe_loc
export append_node
export set_node, end_nodes
export set_node_connection
export struct_loc
export print_node
export find_node_components
export select_results_XYZ_nodes,select_results_nodes


# ========================================================================================

end # module
