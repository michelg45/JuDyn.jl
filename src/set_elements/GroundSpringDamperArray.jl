"""
    GroundSpringDamperArray

    data structure for ground spring element

"""
mutable struct GroundSpringDamperArray



    number::Vector{Int}
    node_order::Vector{Int}
    position::Vector{Vec3}
    params::Vector{Vector{Float64}}
#    uniform_rotation::Vector{Bool}
#    rotation_speed::Vector{Vec3}

    function GroundSpringDamperArray()

        number = Vector{Int}[]
        node_order = Vector{Int}[]
        position = Vector{Vec3}[]
        params = Vector{Vector{Float64}}[]
        return new(number, node_order, position,params)

    #    uniform_rotation = Vector{Bool}[]
    #    rotation_speed = Vector{Vec3}[]

    #    return new(number, node, position,params,uniform_rotation,rotation_speed)
    end


end