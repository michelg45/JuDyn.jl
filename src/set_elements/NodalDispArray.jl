"""
NodalDispArray

        data structure defining a container for nodal displacements

"""
mutable struct NodalDispArray

    number::Vector{Int}
    node_order::Vector{Int}
    axe::Array{Vec3}
    scale_factor::Vector{Float64}
    time_function::Vector{String}
    params::Vector{Vector{Float64}}

    function NodalDispArray()

        number = []
        node_order = []
        axe = []
        scale_factor = []
        time_function = []
        params = []

        return new(number, node_order, axe, scale_factor, time_function, params)

    end


end

# NodalDispArray()= NodalDispArray(Vector{Int}[],Vector{Int}[],Vector{Vec3}[],Vector{String}[],Vector{Vector{Float64}}[])

disp_container = NodalDispArray()
