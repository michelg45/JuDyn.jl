"""
    LinConstrArray

Data structure for the  `lin_constr_container` data set.  It contains the following data for each element of the data set:

|                                   |                                              |    
|:--------------------------------------|:---------------------------------------------| 
| numbers::Vector{Int} | elemnt numbers |
| node_orders::Array{Vector{Int},1} | positions of node frames in "node_container" .| 
| components::Vector{Vector{Int}}| components involved in the linear constraint.| 
| dofs::Vector{Vector{Int}} | degrees of freedom of the constraints. |
| coefs::Vector{Vector{Float64}} | coefficients of the constraint. |
| scale_factor::Vector{Float64} | constraint scaling factors |
                
Creation sequence:
                
````{verbatin}
    global lin_constr_container =  LinConstrArray()
````

"""
mutable struct LinConstrArray

    numbers::Vector{Int}
    node_orders::Vector{Vector{Int}}
    components::Vector{Vector{Int}}
    dofs::Vector{Vector{Int}}
    coefs::Vector{Vector{Float64}}
    val::Vector{Float64}
    scale_factor::Vector{Float64}


    function LinConstrArray()

        numbers = []
        dofs = []
        node_orders = []
        components = []
        coefs = []
        val = []
        scale_factor = []


        return new(numbers, node_orders, components, dofs, coefs, val, scale_factor)

    end


end

lin_constr_container = LinConstrArray()