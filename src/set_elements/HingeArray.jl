"""
    HingeArray

Data structure for the  `hinge_container` data set.  It contains the following data for each element of the data set:

|                                   |                                              |    
|:--------------------------------------|:---------------------------------------------| 
| numbers::Vector{Int} | elemnt numbers |
| node_orders::Array{Vector{Int},1} | positions of nodes 1 and 2 in "node_container" | 
| positions::Array{Vector{Vec3},1} | hinge positions relatively to node frames 1 and 2 | 
| orientations::Array{Vector{RV3},1} | hinge orientations relatively to node frames 1 and 2 |
| axis::Array{RV3,1} | orientations of hinge axis in relative frame |
| scale_factor::Vector{Float64} | constraint scaling factors |
| mode::Array{String,1} | hinge actuation modes |
| time_function::Array{String,1} | time description of actuation in "driven" or "force" modes.
| params::Vector{Vector{Float64}} | hinge properties ("spring" mode) or time function parameters "driven" or "force" modes). |

Creation sequence:

````{verbatin}
    global hinge_container = HingeArray()
````
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