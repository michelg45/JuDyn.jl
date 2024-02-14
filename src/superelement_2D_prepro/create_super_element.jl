function create_super_element(JLD_file::String,SE_name::String,nbr_vibration_modes::Int,Ox::Bool, Oy::Bool, Oz::Bool, Phix::Bool, boundary_nodes::Vector{Int},boundary_node_components::Vector{Int})

println(pwd())

println("Phix ", Phix)

Nrig, Nelast,node_coordinates_CM_frame, Mrig, U, omega2, eigenmodes = vibration_properties(JLD_file)

omega2, eigenmodes = vibration_mode_selection(nbr_vibration_modes,Ox,Oy,Oz,Phix,Nelast,omega2,eigenmodes)
println("size(omega2) ", size(omega2))
println("size(eigenmodes) ", size(eigenmodes))

global vibration_model = VibrationModel(Nrig,nbr_vibration_modes,node_coordinates_CM_frame, Mrig,
U, omega2, eigenmodes)

file = jldopen(JLD_file,"r+")
"vibration_model" in keys(file) && (delete!(file,"vibration_model"))
write(file,"vibration_model",vibration_model)
close(file)

Nodes = size(node_coordinates_CM_frame,2)

println("Nrig ", Nrig)


herting_reduced_model(JLD_file,SE_name,boundary_nodes,boundary_node_components)

end

function create_super_element(JLD_file::String,SE_name::String,nbr_vibration_modes::Int,Ox::Bool, Oy::Bool, Oz::Bool, Phix::Bool)
    boundary_nodes = Int[1; Nodes]
    boundary_node_components = Int[6; 6]
    create_super_element(JLD_file,SE_name,nbr_vibration_modes,Ox, Oy, Oz, Phix, boundary_nodes, boundary_node_components)

end