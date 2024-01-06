"""
    NodeLinkArray

    Data structure for the "Main.SetElements.node_link_container" array. Contains the following data for each node link element of the element set.

        * numbers::Vector{Int}  element number.
        * node_orders::Array{Vector{Int},1} linked nodes. 
        * relative_position::Array{Vec3,1} relative position of node 2  to node 1.
        * scale_factor::Vector{Float64} constraint scaling factor

"""
mutable struct NodeLinkArray

    numbers::Vector{Int}
    node_orders::Array{Vector{Int},1}
    relative_position::Array{Vec3,1}
    scale_factor::Vector{Float64}
end


NodeLinkArray()=NodeLinkArray(Vector{Int}[],Array{Vector{Int},1}[],Array{Vec3,1}[],Vector{Float64}[])