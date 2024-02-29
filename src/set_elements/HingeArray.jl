"""
    HingeArray

Data structure for  hinge element

* elemnt number
* position of node 1 and 2 in "node_container"
* hinge position relatively to nodes 1 and 2
* hinge orientation relatively to nodes 1 and 2
* orientation of hinge axis in relative frame
* driving function label

"""
mutable struct HingeArray


    
        numbers::Vector{Int}
        node_orders::Array{Vector{Int},1}
        positions::Array{Vector{Vec3},1}
        orientations::Array{Vector{RV3},1}
        axis::Array{RV3,1}
        scale_factor::Vector{Float64}
        mode::Array{String,1}
        time_function::Array{String,1}
        params::Vector{Vector{Float64}}
    
        function HingeArray()
            
            numbers = []
            node_orders = []
            positions = []
            orientations = []
            axis = []
            mode = []
            time_function = []
            params = []
            scale_factor = []
    
            return new(numbers,node_orders,positions,orientations,axis,scale_factor,mode,time_function,params)
    
        end
    end
    
    # HingeArray()=HingeArray(Vector{Int}[],Array{Vector{Int},1}[],Array{Vector{Vec3},1}[],Array{Vector{RV3},1}[],Array{RV3,1}[],Array{String,1}[],Array{String,1}[],Vector{Vector{Float64}}[])