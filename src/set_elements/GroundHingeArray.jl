
"""
    GroundHingeArray

    Data structure for ground hinge element
        * elemnt number
        * position of node in "node_container"
        * position of node
        * orientation of node
        * orientations of hinge axis in node frame

"""
mutable struct GroundHingeArray

    
        numbers::Vector{Int}
        node_orders::Vector{Int}
        position::Array{Vec3,1}
        orientation::Array{RV3,1}
        axis::Array{RV3,1}
        scale_factor::Vector{Float64}
        mode::Array{String,1}
        time_function::Array{String,1}
        params::Vector{Vector{Float64}}
    
        function GroundHingeArray()
            
            numbers = []
            node_orders = []
            positions = []
            orientations = []
            axis = []
            scale_factor = []
            mode = []
            time_function = []
            params = []
    
            return new(numbers,node_orders,positions,orientations,axis,scale_factor,mode,time_function,params)
    
        end
    
    end
    
    # GroundHingeArray()=GroundHingeArray(Vector{Int}[],Vector{Int}[],Array{Vec3,1}[],Array{RV3,1}[],Array{RV3,1}[],Array{String,1}[],Array{String,1}[],Vector{Vector{Float64}}[])