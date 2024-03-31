"""
    SphericalJointArray
   
Data structure for `spherical_joint` element
        
|                                   |                                              |    
|:--------------------------------------|:---------------------------------------------| 
| numbers::Vector{Int} | numbering of spherical_joint elements  |
| node_orders::Vector{Vector{Int}} | position of nodes in the node_container array. |  
| relative_positions::Vector{Vector{Vec3}}| relative position of the connecting node inode frames 1 and 2|
| | scale_factor::Vector{Float64} | scaling factors of the constraints |

"""
mutable struct SphericalJointArray

    numbers::Vector{Int}
    node_orders::Vector{Vector{Int}}
    relative_positions::Vector{Vector{Vec3}}
    scale_factor::Vector{Float64}

    
        function SphericalJointArray()
    
            numbers = []
            node_orders = []
            relative_positions = []
            scale_factor = []
    
            return new(numbers, node_orders, relative_positions, scale_factor)
            
        end
    
    end