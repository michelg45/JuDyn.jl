"""
    FrameLinkArray
    
        Data structure for frame_link_element
            * elemnt number
            * position of nodes in "node_container"
            * relative position
            * relative orientation
    
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
    
    # FrameLinkArray()=FrameLinkArray(Vector{Int}[],Array{Vector{Int},1}[],Array{Vec3,1}[],Array{RV3,1}[])
    
