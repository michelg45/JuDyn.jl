"""
    RigidMassArray
    
        structure defining a container for rigid masses

"""
mutable struct RigidMassArray

    numbers::Vector{Int}
    node_orders::Vector{Int}
    masses::Vector{Float64}
    # gravity::Array{Vec3,1}
    # uniform_rotation::Vector{Bool}
    # rotation_speed::Vector{Vec3}

    function RigidMassArray()

        numbers = Vector{Int}[]
        node_orders = Vector{Int}[]
        masses = Vector{Float64}[]

        return new(numbers,node_orders,masses)

        # gravity = Vector{Vec3}[]
        # uniform_rotation = Vector{Bool}[]
        # rotation_speed = Vector{Vec3}[]

        # return new(numbers,node_orders,masses,gravity,uniform_rotation,rotation_speed)
    end

end