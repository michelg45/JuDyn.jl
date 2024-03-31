
"""
    GroundHingeArray

Data structure for the  `ground_hinge_container` data set.  It contains the following data for each element of the data set:

|                                   |                                              |    
|:--------------------------------------|:---------------------------------------------| 
| numbers::Vector{Int} | elemnt numbers |
| node_orders::Array{Vector{Int},1} | positions of node frame in "node_container" | 
| position::Array{Vector{Vec3},1} | hinge position relatively to node frame| 
| orientation::Array{Vector{RV3},1} | hinge orientations relatively to node frame |
| axis::Array{RV3,1} | orientations of hinge axis in relative frame |
| scale_factor::Vector{Float64} | constraint scaling factors |
| mode::Array{String,1} | hinge actuation modes |
| time_function::Array{String,1} | time description of actuation in "driven" or "force" modes.
 params::Vector{Vector{Float64}} | hinge properties ("spring" mode) or time function parameters "driven" or "force" modes). |
            
Creation sequence:
            
````{verbatin}
    global ground_hinge_container = GroundHingeArray()
````
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