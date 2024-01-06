
"""
    NodalTorqueArray

        Data structure defining a container for nodal torques

            * number::Vector{Int} element number
            * node::Vector{Int} node of application
            * axe::Array{Vec3}  torque direction 
            * frame::Vector{String}   "absolute" or "material" type
            * time_function::Vector{String} time function 
            * params::Vector{Vector{Float64}} parameters of the time function

"""
mutable struct NodalTorqueArray



    number::Vector{Int}
    node::Vector{Int}
    axe::Array{Vec3}
    frame::Vector{String}
    time_function::Vector{String}
    params::Vector{Vector{Float64}}

    function NodalTorqueArray()

        number = []
        node = []
        axe = []
        time_function = []
        frame = []
        params = []

        return new(number, node, axe, frame, time_function,  params)

    end


end

# NodalForceArray()= NodalForceArray(Vector{Int}[],Vector{Int}[],Vector{Vec3}[],Vector{String}[],Vector{Vector{Float64}}[])

torque_container = NodalTorqueArray()