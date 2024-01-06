mutable struct SuperElementArray

    numbers::Vector{Int}
    ref_node_orders::Vector{Int}
    boundary_node_numbers::Vector{Int}
    internal_mode_numbers::Vector{Int}
    boundary_node_orders::Vector{Vector{Int}}
    boundary_node_components::Vector{Vector{Int}}
    local_node_coordinates::Vector{Vector{Vec3}}
    local_node_orientations::Vector{Vector{RV3}}
    superelement_types::Vector{String}
    matrix_sets::Vector{Int}
    superelement_names::Vector{String}
    nl_correction::Vector{Bool}

end

SuperElementArray()=SuperElementArray(Vector{Int}[],Vector{Int}[],Vector{Int}[],Vector{Int}[],Vector{Vector{Int}}[],Vector{Vector{Int}}[],Vector{Vector{Vec3}}[],Vector{Vector{RV3}}[],Vector{String}[],Vector{Int}[],Vector{String}[],Vector{Bool}[])
