mutable struct  LinearStructureModel
   Nodes::Int
   Ndofs::Int
   node_coordinates::Array{Float64,2}
   node_components::Vector{Int}
   loc_vector::Vector{Int}
   ip_nodes::Vector{Int}
   stiffness_matrix::Array{Float64,2}
   mass_matrix::Array{Float64,2}
   S_matrix::Vector{Array{Float64,2}}

   function LinearStructureModel(Nodes,Ndofs)
      node_coordinates = zeros(3,Nodes)
      node_components = zeros(Int,Nodes)
      loc_vector = zeros(Int,Ndofs)
      ip_nodes = zeros(Int,Ndofs)
      stiffness_matrix = zeros(Ndofs,Ndofs)
      mass_matrix = zeros(Ndofs,Ndofs)
      S_matrix = [zeros(Ndofs,Ndofs),zeros(Ndofs,Ndofs),zeros(Ndofs,Ndofs)]

      return new(Nodes, Ndofs, node_coordinates,
      node_components, loc_vector, ip_nodes,stiffness_matrix,
      mass_matrix,S_matrix)
   end
 end
