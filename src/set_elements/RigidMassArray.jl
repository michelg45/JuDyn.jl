"""
    RigidMassArray
    
Data structure for the `SetElements.rigid_mass_container` array, created by the `set_rigid_mass` function. It contains the following data for each rigid body element of the element set:

|  |  |
|:-----------------------|:----------------------------  |
| numbers::Vector{Int} | number of the element |
| node_orders::Vector{Int} | order of the element node in the structural node set |
| masses::Vector{Float64} | mass of the rigid body |
                                                      
Creation sequence:
        
````{verbatim}
    global rigid_mass_container = RigidMassArray() 
```` 
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