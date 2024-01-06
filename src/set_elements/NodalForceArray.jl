"""
NodalForceArray

        Data structure for the "Main.SetElements.nodal_force_container" array. Contains the following data
        for each nodal force element of the element set.

            numbers::Vector{Int}                                   number of the element
            node_order::Vector{Int}                                order of the element node in th structural node set
            axe::Array{Vec3}                                       direction of the nodal force (dead load)
            time_function::Vector{String}                          description of the force time dependence
            params::Vector{Vector{Float64}}                        parameters of the time function      
    
            creation sequence:
    
                global force_container = NodalForceArray() 
                
"""
mutable struct NodalForceArray

    """

    structure defining a container for nodal forces

    """

    number::Vector{Int}
    node_order::Vector{Int}
    axe::Array{Vec3}
    time_function::Vector{String}
    params::Vector{Vector{Float64}}

    function NodalForceArray()

        number = []
        node_order = []
        axe = []
        time_function = []
        params = []

        return new(number, node_order, axe, time_function, params)

    end

end

global force_container = NodalForceArray()

