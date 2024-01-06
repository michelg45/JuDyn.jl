"""
    PrismaticJointArray

    Data structure defining a container for prismatic joints

"""
mutable struct PrismaticJointArray

    number::Vector{Int}
    connecting_nodes::Vector{Vector{Int}}
    l_0::Vector{Float64}
    axis::Vector{Vec3}
    scale_factor::Vector{Float64}
    time_function::Vector{String}
    params::Vector{Vector{Float64}}

    function PrismaticJointArray()

        number = []
        connecting_nodes = []
        l_0 = []
        axis = []
        scale_factor = []
        time_function = []
        params = []

        return new(number, connecting_nodes,l_0,axis,scale_factor, time_function, params)

    end


end

# PrismaticJointArray()= PrismaticJointArray(Vector{Int}[],Vector{Int}[],Vector{Vec3}[],Vector{String}[],Vector{Vector{Float64}}[])

prism_container = PrismaticJointArray()