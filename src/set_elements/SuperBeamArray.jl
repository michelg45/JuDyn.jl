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
