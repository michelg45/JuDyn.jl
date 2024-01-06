"""
    InequalityArray

        Data structure for the "Main.SetElements.inequality_container" array. Contains the following data
        for each inequality element of the element set.

            numbers::Vector{Int}                                    number of the element
            node::Vector{Vector{Int}}                               order of the element node in th structural node set
            points::Vector{Vector{Vec3}}                            set of point used to define the constraint
            direction::Vector{Vec3}                                 direction of the constraint
            type::Vector{String}                                    type of constraint: "linear" or "point_to_axis"
            bound_type::Vector{String}                              "upper" or "lower"
            time_function::Vector{String}                           function defining the time evolution of the constraint.
            bound_params::Vector{Vector{Float64}}                   parameters of the time function.
            
            creation sequence:
    
                global inequality_container = InequalityArray() 
"""
mutable struct InequalityArray


    number::Vector{Int}
    node::Vector{Int}
    points::Vector{Vector{Vec3}}
    direction::Vector{Vec3}
    type::Vector{String}
    bound_type::Vector{String}
    time_function::Vector{String}
    bound_params::Vector{Vector{Float64}}

    function InequalityArray()

        number = []
        node = []
        points = []
        direction = []
        type = []
        bound_type = []
        time_function = []
        bound_params = []

        return new(number, node, points, direction, type, bound_type, time_function, bound_params)

    end


end


inequality_container = InequalityArray()

