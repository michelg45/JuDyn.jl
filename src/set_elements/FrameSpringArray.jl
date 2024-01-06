"""
    FrameSpringArray

    Data structure for frame_spring_element
        * elemnt number
        * position of nodes in "node_container"
        * extension stiffness kernel : Vec3 vector
        * rotation stiffness kernel : Vec3 vector
        * initial position of nodes
        * initial orientation of nodes

"""
mutable struct FrameSpringArray

    numbers::Vector{Int}
    node_orders::Array{Vector{Int},1}
    extension_stiffness::Array{Vec3,1}
    rotation_stiffness::Array{Vec3,1}
    initial_displacements::Array{Vec3,1}
    initial_rotations::Array{RV3,1}

    function FrameSpringArray()

        numbers = []
        node-orders = []
        extension_stiffness = []
        rotation_stiffness = []
        initial_displacements = []
        initial_rotations = []

        return new(numbers, node_orders, extension_stiffness, rotation_stiffness,
         initial_displacements, initial_rotations)

    end
end


# FrameSpringArray()=FrameSpringArray(Vector{Int}[],Array{Vector{Int},1}[],Array{Vec3,1}[],Array{Vec3,1}[],Array{Vec3,1}[],Array{RV3,1}[])





