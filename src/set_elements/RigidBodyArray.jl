"""
    RigidBodyArray

Data structure for the `SetElements.rigid_body_container` array, created by the `set_rigid_body` function. It contains the following data for each rigid body element of the element set:

|  |  |
|:-----------------------|:----------------------------  |
| numbers::Vector{Int} | number of the element |
| node_orders::Vector{Int} | order of the element node in the structural node set |
| masses::Vector{Float64} | mass of the rigid body |
|  inertia_tensors::Array{Mat3,1}  | inertia tensor of the rigid body   |
                                              
Creation sequence:

````{verbatim}
    global rigid_body_container = RigidBodyArray() 
````    
                
"""
mutable struct RigidBodyArray

    numbers::Vector{Int}
    node_orders::Vector{Int}
    masses::Vector{Float64}
    inertia_tensors::Array{Mat3,1}
    # gravity::Array{Vec3,1}
    # uniform_rotation::Vector{Bool}
    # rotation_speed::Vector{Vec3}

    function RigidBodyArray()

        numbers = Vector{Int}[]
        node_orders = Vector{Int}[]
        masses = Vector{Float64}[]
        inertia_tensors = Vector{Mat3}[]
        # gravity = Vector{Vec3}[]

        return new(numbers,node_orders,masses,inertia_tensors)

        # uniform_rotation = Vector{Bool}[]
        # rotation_speed = Vector{Vec3}[]

        # return new(numbers,node_orders,masses,inertia_tensors,gravity,uniform_rotation,rotation_speed)
    end
end