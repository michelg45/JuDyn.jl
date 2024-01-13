"""
    SuperBeamArray

        Data structure for the "Main.SetElements.super_beam_container" array. Contains the following data
        for each super_beam element of the element set.

            numbers::Vector{Int}                             element numbers
            ref_node_orders::Vector{Int}                     order of the reference nodes in the structural node set
            boundary_node_orders::Vector{Vector{Int}}        order of the connection nodes in the structural node set
            local_node_orientations::Array{Vector{RV3},1}    local reference frames at boundary nodes
            length::Vector{Float64}                          superelement lengths
            Jrot::Vector{Vector{Float64}}                    superelement tensors of inertia
            K_elast::Vector{Array{Float64,2}}                superelement stiffness matrices                     
            M_elast::Vector{Array{Float64,2}}                superelement mass matrices
            S::Vector{Vector{}}                              superelement gyroscopic matrices
            nl_correction::Vector{Bool}                      nonlinear correction parameters
"""
mutable struct SuperBeamArray

    numbers::Vector{Int}
    ref_node_orders::Vector{Int}
    boundary_node_orders::Vector{Vector{Int}}
    local_node_orientations::Array{Vector{RV3},1}
    length::Vector{Float64}
    mass::Vector{Float64}
    Jrot::Vector{Vector{Float64}}
    K_elast::Vector{Array{Float64,2}}
    M_elast::Vector{Array{Float64,2}}
    S::Vector{Vector{}}
    nl_correction::Vector{Bool}



end


SuperBeamArray()=SuperBeamArray(Vector{Int}[],Vector{Int}[],Vector{Vector{Int}}[],Array{Vector{RV3},1}[],Vector{Float64}[],Vector{Float64}[],Vector{Vector{Float64}}[],Vector{Array{Float64,2}}[],Vector{Array{Float64,2}}[],Vector{Vector{Any}}[],Vector{Bool}[])
