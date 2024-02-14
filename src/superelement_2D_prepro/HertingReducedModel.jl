mutable struct  HertingReducedModel
    
  Gyro::Bool
    Nrig::Int
    nbr_boundary_dofs::Int
    nbr_vibration_modes::Int
    nbr_boundary_nodes::Int
    boundary_node_coordinates::Array{Float64,2}
    boundary_nodes::Vector{Int}
    boundary_node_components::Vector{Int}
    Mrig::Array{Float64,2}
    K::Array{Float64,2}
    M::Array{Float64,2}
    S_rot::Vector{Array{Float64,2}}
    A_rot::Vector{Array{Float64,2}}
  
end
