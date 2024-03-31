"""
    FrameLinkArray

   
Data structure for the  `frame_link_container` data set.  It contains the following data for each element of the data set:
        
|                                   |                                              |    
|:----------------------------------|:---------------------------------------------| 
| numbers::Vector{Int} | numbering of frame_link elements  |
| node_orders::Vector{Vector{Int}} | position of nodes in the node_container array. |  
| relative_position::Vector{Vec3} | relative position of the nodes|
| relative_orientation::Vector{RV3} | relative orientation of the nodes |
| scale_factor::Vector{Float64} | scaling factors of the constraints |


Creation sequence: 

````{verbatin}
        global frame_link_container = FrameLinkArray()
````
"""
mutable struct FrameLinkArray

    numbers::Vector{Int}
    node_orders::Array{Vector{Int},1}
    relative_position::Array{Vec3,1}
    relative_orientation::Array{RV3,1}
    scale_factor::Vector{Float64}

    
        function FrameLinkArray()
    
            numbers = []
            node_orders = []
            relative_position = []
            relative_orientation = []
            scale_factor = []
    
            return new(numbers, node_orders, relative_position, relative_orientation,scale_factor)
            
        end
    
    end